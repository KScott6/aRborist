## Overview

aRborist is an automated sequence and metadata harvester designed to simplify the process of gathering and organizing sequence data from the **NCBI nucleotide database**. It retrieves accessions for specified taxa, extracts and standardizes metadata, and prepares sequences and metadata for downstream analyses. After using aRborist to pull metadata/sequence data, you can use other aRborist functions to help you make a phylogenetic tree. Or assign host information to your taxa of interest.

- **aRborist Input:** a list of taxa (e.g. a list of species names) and a few simple options (loci of interest, NCBI API, etc.)  
- **aRborist Output:** curated metadata tables, phylogenetic trees, host percentage


## Requirements

- R (≥ 4.2) and RStudio (optional but recommended)
- Internet access to query NCBI
- (Recommended) an NCBI API key

> Disclaimer and Limitations: The accuracy and completeness of your results depend on the quality of metadata available in NCBI. While aRborist applies consistent naming, standardization, and error-handling routines, it cannot correct for missing, inconsistent, or ambiguous source data. I encourage users to review and, if necessary, manually refine the curation rules for your own use.

---

## aRborist first time setup

### 1) Get the code

In R: 
```r
install.packages("usethis") # if needed
usethis::create_from_github("KScott6/aRborist")
```

Or, on GitHub, click Code → Download ZIP, unzip somewhere convenient.

In RStudio: File → Open Project… and select the unzipped folder (or set the working directory there).

### 2) (Recommended) Set your NCBI API key

This increases allowed request rates (10 requests/second with an API key, 3 requests/second without). Get your API key by logging into your NCBI account, open Account settings, and scroll down to "API Key Management". You can then copy your API key.

Create (or edit) a file named ~/.Renviron and add:

```bash
NCBI_API_KEY=YOUR_KEY_HERE
```

Save this file and restart your R instance. Then in R:

```R
Sys.getenv("NCBI_API_KEY")  # should show your key (or at least not be empty)
```

### 3) Install required packages

```R
source("prepare.R")   # installs any missing packages required by aRborist
```

### 4) Load aRborist helpers

```R
source(file.path("R","arborist_helpers.R"))
```

### 4) Set options

