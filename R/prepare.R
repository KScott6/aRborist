# prepare.R 
# one-time dependency setup for aRborist

# 1) Set a CRAN mirror if none is set 
if (is.null(getOption("repos")) || getOption("repos")["CRAN"] == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

# 2) CRAN packages
cran_pkgs <- c(
  "rentrez","stringr","plyr","dplyr","withr","XML",
  "data.table","tidyr","phylotools","scales",
  "purrr","readr","phytools","RColorBrewer",
  "maps","ggplot2","tidygeocoder","ggrepel","taxize","yaml"
  # NOTE: Biostrings is from Bioconductor, handled below
)

missing_cran <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(missing_cran)) install.packages(missing_cran)

# 3) Bioconductor packages (Biostrings)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
bioc_pkgs <- c("Biostrings")
missing_bioc <- setdiff(bioc_pkgs, rownames(installed.packages()))
if (length(missing_bioc)) BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)

message("aRborist dependencies have been installed.")
