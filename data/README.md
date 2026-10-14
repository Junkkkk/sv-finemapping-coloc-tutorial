# Workshop data

Locus-level files cut from the SPHERE Twin of the Stanford ADRC / SAMS long-read cohort, plus public GWAS
and 1000 Genomes data for the same regions. All coordinates are GRCh38. People are identified by Twin IDs
(`T0001`, ...), which do not correspond to real participants.

The Twin was generated with SPHERE
([He et al., bioRxiv 2026](https://www.biorxiv.org/content/10.64898/2026.09.01.748580v2)) and is released under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The files here are a subset of that release (three
regions) with covariate PCs computed by notebook 00.
Please acknowledge: "Stanford Alzheimer's Disease Research Center, NIH/NIA grant P30 AG066515."

The full Twin release is on Hugging Face:
[`jpark01/sv-qtl-coloc-data`](https://huggingface.co/datasets/jpark01/sv-qtl-coloc-data), folder
`sv-qtl-coloc/input`, commit `02821371c26603cc9460396582eb4ebfbc3257df`. `data_prep/build_workshop_data.py`
builds this folder from it.

## `loci.tsv`

One row per analysis. Columns: `locus`, `chrom`, `region_start`, `region_end` (the region cut for this
locus), `layer` (molecular layer), `pheno_id` (molecular trait), `label`, `sv_id` (SV highlighted in the paper),
`gwas` (AD or PD, empty for HP), `gwas_lead_pos` (GWAS lead SNV the authors colocalized with the trait;
Extended Data Table 2 and Supplementary Figure 14).

| Locus | Region | Traits |
|---|---|---|
| HP | chr16:71,054,505-73,058,850 | HP protein, two aptamers, plasma and CSF |
| COG7 | chr16:22,462,828-24,505,148 | methylation segment `chr16_10_113` near GGA2/COG7 |
| MAPT | chr17:44,915,490-47,511,509 | `LRRC37A2` expression in PBMC |

Each region covers the trait ±1 Mb, the highlighted SV ±1 Mb, and every GWAS lead ±500 kb.

## `genotypes/`

- `<locus>.SV.tsv.gz`: SVs overlapping the region ±1 Mb, 544 people. Columns `#chr, start, end, svlen, variant_id`,
  then one dosage column per person. `svlen` is negative for deletions. IDs are Sniffles2 IDs
  (`Sniffles2.<TYPE>.<id>`).
- `<locus>.SNV.tsv.gz`: SNVs and indels from short-read genomes inside the region, 416 people. Columns
  `chrom, pos, variant_id, counted_allele, other_allele`, then one dosage column per person. `variant_id` is
  `<chromosome number>:pos:counted_allele:other_allele`; dosages count `counted_allele`.

All Twin values, including genotype dosages, are continuous numbers.

## `phenotypes/`

`<Layer>.phenotypes.tsv`: all traits of each layer (`Plasma`, `CSF`, `RNA`, `Methylation`) that lie inside the
three regions.
Columns `chrom, start, end, pheno_id`, then one normalized value per person.

## `covariates/`

`<Layer>.covariates.tsv.gz`, one row per person in the layer: `ID, Sex, Age, AD_Dx, PD_Dx`, then SV-based
ancestry PCs (`Ancestry_PC1..10`; CSF `..5`), then molecular-phenotype PCs (`PC1..`; plasma 120, CSF 100,
RNA 120, methylation 80), in the column order of the authors' covariate files. Produced by
`notebooks/00_covariates_from_full_twin.ipynb`.

| Layer | People |
|---|---|
| Plasma | 490 |
| CSF | 225 |
| RNA | 268 |
| Methylation | 544 |

## `gwas/`

`<AD|PD>.<locus>.tsv.gz`, formatted like the authors' GWAS files: columns
`rsid, chr, snp_pos, pvalue, effect_allele, non_effect_allele, beta, se`; variants with effect-allele frequency
between 0.02 and 0.98.

- AD: Bellenguez *et al.* 2022, GWAS Catalog GCST90027158 (GRCh38).
- PD: Kim *et al.* 2024, GWAS Catalog GCST90275127 (random-effects p-values), lifted from GRCh37 to GRCh38
  with the UCSC `hg19ToHg38` chain (`data_prep/chain_liftover.py`).

## `reference/`

`1000G.<locus>.tsv.gz`: 1000 Genomes Phase 3 sites in the colocalization windows (COG7, MAPT), with columns
`#CHROM, POS, ID, REF, ALT, AF`, from the GRCh38 release used in the authors' colocalization configuration.
