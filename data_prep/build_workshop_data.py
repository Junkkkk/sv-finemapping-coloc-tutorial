"""
Purpose: Build the small locus-level data set used in the workshop notebooks from the
         full SPHERE Twin release, the AD/PD GWAS files and 1000 Genomes allele frequencies.
Inputs:  --twin-dir   Twin release folder (HF jpark01/sv-qtl-coloc-data, sv-qtl-coloc/input,
                      commit 02821371c26603cc9460396582eb4ebfbc3257df)
         --chain      UCSC hg19ToHg38.over.chain.gz (PD GWAS is GRCh37)
         --kg-dir     1000 Genomes GRCh38 (20170504 liftover release) VCF region slices, one per
                      colocalization window (GWAS lead +-500 kb), named <anything>.reference.tsv;
                      see data_prep/README.md
         --covariates folder with <Layer>.covariates.tsv.gz written by notebook 00
                      (results/covariates_recomputed)
         --out        output folder (data)
Outputs: data/ (loci.tsv, genotypes/, phenotypes/, covariates/, gwas/, reference/)
Depends: numpy, pandas, h5py; chain_liftover.py (same folder)
"""
import argparse
import glob
import gzip
import os
import shutil

import h5py
import numpy as np
import pandas as pd

from chain_liftover import lift_points, read_chain

ONE_MB = 1_000_000
COLOC_HALF_WINDOW = 500_000
GWAS_EAF_MIN, GWAS_EAF_MAX = 0.02, 0.98   # authors' GWAS formatting
FLOAT_FORMAT = "%.7g"                     # Twin values are float32

# One row per analysis: a molecular trait, its lead SV from the paper, and the GWAS leads
# the authors colocalized it with (Extended Data Table 2; Supplementary Figure 14).
ANALYSES = [
    # locus, layer, pheno_id, label, sv_id, gwas, gwas_lead_pos
    ("HP", "Plasma", "ENSG00000257017:HP.3054.3.2", "HP protein (plasma, aptamer 3054-3)", "Sniffles2.DEL.14B0MF", "", ""),
    ("HP", "Plasma", "ENSG00000257017:HP.7905.30.3", "HP protein (plasma, aptamer 7905-30)", "Sniffles2.DEL.14B0MF", "", ""),
    ("HP", "CSF", "HP.3054.3.2:ENSG00000257017:HP", "HP protein (CSF, aptamer 3054-3)", "Sniffles2.DEL.14B0MF", "", ""),
    ("HP", "CSF", "HP.7905.30.3:ENSG00000257017:HP", "HP protein (CSF, aptamer 7905-30)", "Sniffles2.DEL.14B0MF", "", ""),
    ("COG7", "Methylation", "chr16_10_113", "Methylation segment near GGA2/COG7", "Sniffles2.DEL.8FEMF", "AD", 23466676),
    ("MAPT", "RNA", "ENSG00000238083:LRRC37A2", "LRRC37A2 expression (PBMC)", "Sniffles2.INS.1142M10", "AD", 46613342),
    ("MAPT", "RNA", "ENSG00000238083:LRRC37A2", "LRRC37A2 expression (PBMC)", "Sniffles2.INS.1142M10", "AD", 46720553),
    ("MAPT", "RNA", "ENSG00000238083:LRRC37A2", "LRRC37A2 expression (PBMC)", "Sniffles2.INS.1142M10", "PD", 46276519),
    ("MAPT", "RNA", "ENSG00000238083:LRRC37A2", "LRRC37A2 expression (PBMC)", "Sniffles2.INS.1142M10", "PD", 46783637),
]
LOCUS_CHROM = {"HP": "chr16", "COG7": "chr16", "MAPT": "chr17"}


def read_bed_rows(path, ids):
    bed = pd.read_csv(path, sep="\t")
    key = bed.columns[3]
    out = bed[bed[key].isin(ids)]
    assert set(out[key]) == set(ids), (path, set(ids) - set(out[key]))
    return out.rename(columns={key: "pheno_id"})


def write_tsv_gz(df, path):
    with gzip.open(path, "wt", compresslevel=9) as f:
        df.to_csv(f, sep="\t", index=False, float_format=FLOAT_FORMAT)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--twin-dir", required=True)
    ap.add_argument("--chain", required=True)
    ap.add_argument("--kg-dir", required=True)
    ap.add_argument("--covariates", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = a.out
    for sub in ["genotypes", "phenotypes", "covariates", "gwas", "reference"]:
        os.makedirs(os.path.join(out, sub), exist_ok=True)

    analyses = pd.DataFrame(ANALYSES, columns=["locus", "layer", "pheno_id", "label", "sv_id", "gwas", "gwas_lead_pos"])
    analyses["chrom"] = analyses.locus.map(LOCUS_CHROM)

    # Phenotype coordinates fix the cis windows (phenotype start - 1 Mb, end + 1 Mb).
    pheno_rows = {}
    for layer, ids in analyses.groupby("layer").pheno_id:
        rows = read_bed_rows(os.path.join(a.twin_dir, layer, "phenotype.bed"), sorted(set(ids)))
        pheno_rows[layer] = rows
        rows.to_csv(os.path.join(out, "phenotypes", f"{layer}.phenotypes.tsv"), sep="\t", index=False, float_format=FLOAT_FORMAT)
    coords = pd.concat([r[["chrom", "start", "end", "pheno_id"]] for r in pheno_rows.values()])
    analyses = analyses.merge(coords.rename(columns={"chrom": "pheno_chrom", "start": "pheno_start", "end": "pheno_end"}), on="pheno_id")
    assert (analyses.pheno_chrom == analyses.chrom).all()

    # SV coordinates for the paper's lead SVs.
    sv_tables = {}
    for chrom in sorted(set(analyses.chrom)):
        sv_tables[chrom] = pd.read_csv(os.path.join(a.twin_dir, "genotype", "SV", f"SV.{chrom}.tsv.gz"), sep="\t")
    sv_pos = pd.concat([t[["#chr", "start", "end", "variant_id"]] for t in sv_tables.values()])
    analyses = analyses.merge(sv_pos.rename(columns={"variant_id": "sv_id", "start": "sv_start", "end": "sv_end"})[["sv_id", "sv_start", "sv_end"]], on="sv_id")

    # Region per locus = union of the phenotype cis windows (+-1 Mb), the lead SV +-1 Mb
    # (best-tagging SNV search) and the colocalization windows (GWAS lead +-500 kb).
    regions = []
    for locus, g in analyses.groupby("locus"):
        lo = min((g.pheno_start - ONE_MB).min(), (g.sv_start - ONE_MB).min())
        hi = max((g.pheno_end + ONE_MB).max(), (g.sv_end + ONE_MB).max())
        leads = pd.to_numeric(g.gwas_lead_pos, errors="coerce").dropna()
        if len(leads):
            lo, hi = min(lo, leads.min() - COLOC_HALF_WINDOW), max(hi, leads.max() + COLOC_HALF_WINDOW)
        regions.append((locus, g.chrom.iloc[0], int(max(lo, 0)), int(hi)))
    regions = pd.DataFrame(regions, columns=["locus", "chrom", "region_start", "region_end"])
    analyses = analyses.merge(regions, on=["locus", "chrom"])
    cols = ["locus", "chrom", "region_start", "region_end", "layer", "pheno_id", "label", "sv_id", "gwas", "gwas_lead_pos"]
    analyses[cols].to_csv(os.path.join(out, "loci.tsv"), sep="\t", index=False)
    print(regions.to_string(index=False))

    # Genotypes: SVs overlapping the region (all 544 people) and SNV dosages inside it (416 people).
    for _, r in regions.iterrows():
        sv = sv_tables[r.chrom]
        sv = sv[(sv.start <= r.region_end) & (sv.end >= r.region_start)]
        write_tsv_gz(sv, os.path.join(out, "genotypes", f"{r.locus}.SV.tsv.gz"))
        num = r.chrom.replace("chr", "")
        var = pd.read_csv(os.path.join(a.twin_dir, "genotype", "SNV", f"SNV.{r.chrom}.variants.tsv"), sep="\t")
        idx = np.flatnonzero((var.pos >= r.region_start) & (var.pos <= r.region_end))
        with h5py.File(os.path.join(a.twin_dir, "genotype", "SNV", f"SNV.{r.chrom}.h5"), "r") as h:
            dos = h["dosage"][idx[0]:idx[-1] + 1, :]
            samples = [s.decode() for s in h["dosage"].attrs["samples"]]
        assert len(idx) == idx[-1] - idx[0] + 1
        snv = var.iloc[idx].reset_index(drop=True)
        snv.insert(0, "chrom", "chr" + snv.pop("chr").astype(str))
        snv = snv[["chrom", "pos", "id", "counted_allele", "other_allele"]].rename(columns={"id": "variant_id"})
        snv = pd.concat([snv, pd.DataFrame(dos, columns=samples)], axis=1)
        write_tsv_gz(snv, os.path.join(out, "genotypes", f"{r.locus}.SNV.tsv.gz"))
        print(f"{r.locus}: {len(sv)} SVs, {len(snv)} SNVs ({num})")

    # Covariates (ancestry PCs + phenotype PCs) computed by notebook 00 on the full Twin.
    for f in sorted(glob.glob(os.path.join(a.covariates, "*.covariates.tsv.gz"))):
        shutil.copy(f, os.path.join(out, "covariates", os.path.basename(f)))

    # GWAS, formatted like the authors' files: EAF 0.02-0.98, GRCh38 positions.
    gwas_cols = ["rsid", "chr", "snp_pos", "pvalue", "effect_allele", "non_effect_allele", "beta", "se"]
    want = [(chrom.replace("chr", ""), lo, hi) for _, (_, chrom, lo, hi) in regions.iterrows()]
    want_chroms = sorted(set(w[0] for w in want))

    ad = []
    for c in pd.read_csv(os.path.join(a.twin_dir, "GWAS", "GCST90027158.tsv.gz"), sep="\t", chunksize=2_000_000,
                         dtype={"chromosome": str}):
        for chrom, lo, hi in want:
            ad.append(c[(c.chromosome == chrom) & c.base_pair_location.between(lo, hi)])
    ad = pd.concat(ad).drop_duplicates()
    ad = ad[ad.effect_allele_frequency.between(GWAS_EAF_MIN, GWAS_EAF_MAX)]
    ad = ad.rename(columns={"variant_id": "rsid", "chromosome": "chr", "base_pair_location": "snp_pos", "p_value": "pvalue",
                            "other_allele": "non_effect_allele", "standard_error": "se"})[gwas_cols]

    blocks = read_chain(a.chain)
    pd_rows = []
    for c in pd.read_csv(os.path.join(a.twin_dir, "GWAS", "GCST90275127.tsv.gz"), sep="\t", chunksize=2_000_000,
                         dtype={"chromosome": str}):
        c = c[c.chromosome.isin(want_chroms)]
        if not len(c):
            continue
        lifted = lift_points(blocks, ("chr" + c.chromosome).to_numpy(), c.base_pair_location.to_numpy())
        c = c.assign(chr38=lifted["chr"].to_numpy(), pos38=lifted["pos"].to_numpy(), status=lifted["status"].to_numpy())
        c = c[c.status != "unmapped"]
        c = c[c.chr38 == "chr" + c.chromosome]
        for chrom, lo, hi in want:
            pd_rows.append(c[(c.chromosome == chrom) & c.pos38.between(lo, hi)])
    pdg = pd.concat(pd_rows).drop_duplicates()
    pdg = pdg[pdg.effect_allele_frequency.between(GWAS_EAF_MIN, GWAS_EAF_MAX)]
    pdg = pdg.assign(snp_pos=pdg.pos38, chr=pdg.chromosome).rename(columns={
        "rs_id": "rsid", "p_value": "pvalue", "other_allele": "non_effect_allele", "standard_error": "se"})[gwas_cols]
    # Some PD rows carry a number instead of an rsID.
    pdg.loc[~pdg.rsid.astype(str).str.startswith("rs"), "rsid"] = np.nan

    for name, g in [("AD", ad), ("PD", pdg)]:
        for _, r in regions.iterrows():
            s = g[(g.chr == r.chrom.replace("chr", "")) & g.snp_pos.between(r.region_start, r.region_end)].sort_values("snp_pos")
            write_tsv_gz(s, os.path.join(out, "gwas", f"{name}.{r.locus}.tsv.gz"))
            print(f"GWAS {name} {r.locus}: {len(s)} variants")

    # 1000 Genomes reference alleles and frequencies for the colocalization windows.
    kg = []
    for f in sorted(glob.glob(os.path.join(a.kg_dir, "*.reference.tsv"))):
        v = pd.read_csv(f, sep="\t", dtype={"#CHROM": str})
        v["AF"] = v.INFO.str.extract(r"(?:^|;)AF=([^;]+)", expand=False)
        kg.append(v[["#CHROM", "POS", "ID", "REF", "ALT", "AF"]])
    kg = pd.concat(kg).drop_duplicates()
    for _, r in regions.iterrows():
        s = kg[(kg["#CHROM"] == r.chrom.replace("chr", "")) & kg.POS.between(r.region_start, r.region_end)].sort_values("POS")
        if len(s):
            write_tsv_gz(s, os.path.join(out, "reference", f"1000G.{r.locus}.tsv.gz"))
            print(f"1000G {r.locus}: {len(s)} records, {s.POS.min()}-{s.POS.max()}")


if __name__ == "__main__":
    main()
