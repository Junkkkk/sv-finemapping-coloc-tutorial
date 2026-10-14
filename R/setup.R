# Shared setup for the workshop notebooks.
# Finds the repository folder, installs missing CRAN packages, and defines paths.

find_repo_dir <- function() {
    d <- normalizePath(getwd())
    for (i in 1:5) {
        if (file.exists(file.path(d, "R", "setup.R")) && dir.exists(file.path(d, "data"))) return(d)
        d <- dirname(d)
    }
    stop("Run the notebooks from inside the workshop repository (or its notebooks/ folder).")
}

REPO_DIR <- find_repo_dir()
DATA_DIR <- file.path(REPO_DIR, "data")
RESULTS_DIR <- file.path(REPO_DIR, "results")
dir.create(RESULTS_DIR, showWarnings = FALSE)

needed <- c("data.table", "R.utils", "ggplot2", "RcppEigen", "susieR", "coloc")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
    message("Installing: ", paste(missing, collapse = ", "), " (a few minutes the first time)")
    repos <- "https://cloud.r-project.org"
    # On Linux (e.g. Google Colab) use precompiled packages from the Posit Package Manager.
    if (Sys.info()[["sysname"]] == "Linux" && file.exists("/etc/os-release")) {
        os <- readLines("/etc/os-release")
        codename <- sub("^VERSION_CODENAME=", "", grep("^VERSION_CODENAME=", os, value = TRUE))
        if (length(codename) == 1 && nzchar(codename)) {
            repos <- sprintf("https://packagemanager.posit.co/cran/__linux__/%s/latest", codename)
            options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(),
                                            paste(getRversion(), R.version$platform, R.version$arch, R.version$os)))
        }
    }
    install.packages(missing, repos = repos, Ncpus = max(1, parallel::detectCores()))
}

Sys.setenv(LANGUAGE = "en")
suppressWarnings(suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(RcppEigen)
    library(susieR)
    library(coloc)
}))

options(repr.plot.width = 10, repr.plot.height = 4.5, width = 160)
theme_set(theme_minimal(base_size = 12))

# The analysis loci used in the workshop.
LOCI <- if (file.exists(file.path(DATA_DIR, "loci.tsv"))) fread(file.path(DATA_DIR, "loci.tsv")) else NULL

# Read a genotype region file: 5 annotation columns, then one column per person.
read_genotypes <- function(path) {
    g <- fread(path)
    list(info = g[, 1:5], mat = as.matrix(g[, 6:ncol(g)], rownames.value = g$variant_id))
}
