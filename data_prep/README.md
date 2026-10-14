# How `data/` was built

Participants do not need this folder. It records how the workshop files were made from public sources.

1. Download the Twin release from Hugging Face
   ([`jpark01/sv-qtl-coloc-data`](https://huggingface.co/datasets/jpark01/sv-qtl-coloc-data), folder
   `sv-qtl-coloc/input`, commit `02821371c26603cc9460396582eb4ebfbc3257df`).
2. Run `notebooks/00_covariates_from_full_twin.ipynb` with `TWIN_FULL_DIR` pointing to that folder. It writes
   `results/covariates_recomputed/<Layer>.covariates.tsv.gz`.
3. Get the UCSC chain `hg19ToHg38.over.chain.gz`.
4. Make one 1000 Genomes slice per colocalization window (GWAS lead ±500 kb) from the GRCh38 liftover release
   `ALL.chr<N>_GRCh38.genotypes.20170504.vcf.gz`, keeping the header line and the first 8 VCF columns:

   ```bash
   tabix -h ALL.chr16_GRCh38.genotypes.20170504.vcf.gz 16:22966676-23966676 | grep -v '^##' | cut -f1-8 > COG7_AD_16_23466676.reference.tsv
   ```

5. Build:

   ```bash
   python build_workshop_data.py --twin-dir <twin>/sv-qtl-coloc/input --chain hg19ToHg38.over.chain.gz \
       --kg-dir <slices> --covariates ../results/covariates_recomputed --out ../data
   ```

Requires Python 3 with numpy, pandas and h5py.
