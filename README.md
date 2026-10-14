# Fine-mapping and colocalization with structural variants

Hands-on notebooks for the ASHG 2026 workshop **"Fine-mapping and Colocalization with Structural
Variants: Practical Workflows Using GWAS and Multi-omic Summary Statistics"**.

The notebooks walk through the structural-variant (SV) analyses of

> Jensen TD, Le Guen Y, Talozzi L, Yang S, *et al.* Long-read genome sequencing and multi-omics in aging
> and neurodegeneration. *medRxiv* (2025). doi:[10.1101/2025.10.10.25337775](https://doi.org/10.1101/2025.10.10.25337775)

using the authors' analysis code ([tjense25/ADRC-Flagship-lrGS-Code](https://github.com/tjense25/ADRC-Flagship-lrGS-Code))
on a **SPHERE Twin** of the Stanford ADRC / SAMS long-read cohort
([He et al., bioRxiv 2026](https://www.biorxiv.org/content/10.64898/2026.09.01.748580v2)). The Twin contains no
real person's record, so it can be shared openly, while it keeps the associations needed for QTL mapping,
fine-mapping and colocalization.

## Notebooks

| Notebook | What you do |
|---|---|
| [`01_sv_candidates_at_gwas_loci`](notebooks/01_sv_candidates_at_gwas_loci.ipynb) | Find SVs near AD/PD GWAS leads, their allele frequency and length, LD with the lead SNV, best-tagging SNVs; the HP deletion that SNVs tag poorly |
| [`02_finemapping_and_colocalization`](notebooks/02_finemapping_and_colocalization.ipynb) | SV-QTL mapping in four molecular layers, joint SV + SNV fine-mapping with SuSiE, colocalization with AD/PD GWAS (`coloc.abf`) with regional and LocusCompare plots, and the "SV in credible set + colocalization" call |
| [`00_covariates_from_full_twin`](notebooks/00_covariates_from_full_twin.ipynb) | Optional: download the full Twin (~550 MB) and recompute the ancestry and phenotype PCs used as covariates |

The notebooks are saved with their outputs, so you can read them without running anything.
Notebooks 01 and 02 use the small locus files in [`data/`](data/README.md) (about 20 MB); the covariates there
were produced by notebook 00, and you can rerun it to check that you get the same numbers.

## What the notebooks reproduce

| Locus | Finding in the paper |
|---|---|
| **HP** protein, plasma and CSF | A 1.7 kb deletion is the lead variant, alone in its credible set, and associates more strongly with CSF HP than its best-tagging SNV rs217181 |
| **GGA2 / COG7** methylation segment, AD | The segment colocalizes with AD, and a 208 bp tandem-repeat deletion near the COG7 promoter is in its credible set |
| **LRRC37A2** expression (MAPT region), AD and PD | Expression colocalizes with AD (PP.H4 = 0.79) and PD (PP.H4 = 0.98), with a 37 bp insertion as the lead variant of the credible set |

The notebooks reproduce these findings on the SPHERE Twin; exact numbers differ slightly from the paper.

## Running the notebooks

**Google Colab.** Open a notebook in Colab
([01](https://colab.research.google.com/github/Junkkkk/sv-finemapping-coloc-tutorial/blob/main/notebooks/01_sv_candidates_at_gwas_loci.ipynb),
[02](https://colab.research.google.com/github/Junkkkk/sv-finemapping-coloc-tutorial/blob/main/notebooks/02_finemapping_and_colocalization.ipynb),
[00](https://colab.research.google.com/github/Junkkkk/sv-finemapping-coloc-tutorial/blob/main/notebooks/00_covariates_from_full_twin.ipynb)),
choose *Runtime → Change runtime type → R*, and run all cells. The first cell downloads this repository and
installs missing R packages (a few minutes the first time). Run the first cell of notebook 02 while you work
on notebook 01, so its packages are ready.

**Your own computer.** R ≥ 4.3 and Jupyter with the R kernel:

```bash
pip install jupyterlab
Rscript -e 'install.packages("IRkernel", repos = "https://cloud.r-project.org"); IRkernel::installspec()'
git clone https://github.com/Junkkkk/sv-finemapping-coloc-tutorial.git
cd sv-finemapping-coloc-tutorial/notebooks
jupyter lab
```

R packages: `data.table`, `R.utils`, `ggplot2`, `RcppEigen`, `susieR`, `coloc` (installed by `R/setup.R` if
missing). The notebooks were run with R 4.3.1, susieR 0.14.2 and coloc 5.2.3.

## Contact

- Workshop: Tanner D. Jensen (tannerj@stanford.edu), Yann Le Guen (yleguen@stanford.edu)
- Notebooks and data: Junyoung Park (jpark01@stanford.edu)
- Greicius Lab, Stanford University: Michael D. Greicius (PI), greicius@stanford.edu, https://greiciuslab.stanford.edu

## Acknowledgments

- The notebooks and the SPHERE Twin data were prepared by Junyoung Park for the ASHG 2026 workshop organized and
  presented by Yann Le Guen and Tanner D. Jensen.
- SPHERE was developed in Zihuai He's group at Stanford University.
- The analysis code is adapted from [tjense25/ADRC-Flagship-lrGS-Code](https://github.com/tjense25/ADRC-Flagship-lrGS-Code)
  by Tanner D. Jensen. Colocalization inputs are prepared as in
  [mikegloudemans/alzheimers-parkinsons-colocalization](https://github.com/mikegloudemans/alzheimers-parkinsons-colocalization).
- Stanford Alzheimer's Disease Research Center, NIH/NIA grant P30 AG066515.

## Citation

If you use these notebooks, please cite the original study and SPHERE:

- Jensen TD, Le Guen Y, Talozzi L, Yang S, *et al.* Long-read genome sequencing and multi-omics in aging and
  neurodegeneration. *medRxiv* (2025). doi:[10.1101/2025.10.10.25337775](https://doi.org/10.1101/2025.10.10.25337775)
- He *et al.* SPHERE. *bioRxiv* (2026). [10.64898/2026.09.01.748580](https://www.biorxiv.org/content/10.64898/2026.09.01.748580v2)

Data sources: Alzheimer's disease GWAS, Bellenguez *et al.*, *Nat Genet* 2022
(doi:[10.1038/s41588-022-01024-z](https://doi.org/10.1038/s41588-022-01024-z), GCST90027158); Parkinson's
disease GWAS, Kim *et al.*, *Nat Genet* 2024 (doi:[10.1038/s41588-023-01584-8](https://doi.org/10.1038/s41588-023-01584-8),
GCST90275127); allele frequencies, 1000 Genomes Project (GRCh38).

## License

The analysis code adapted from tjense25/ADRC-Flagship-lrGS-Code is under the MIT License (notice in
[`LICENSES/`](LICENSES/)). The data in `data/` are derived from the SPHERE Twin release and shared under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
