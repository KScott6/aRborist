## Overview

aRborist is an automated sequence and metadata harvester designed to simplify the process of gathering and organizing sequence data from the **NCBI nucleotide database**. It retrieves accessions for specified taxa, extracts and standardizes metadata, and prepares sequences and metadata for downstream analyses. After using aRborist to pull metadata/sequence data, you can use other aRborist functions to help you make a phylogenetic tree. Or assign host information to your taxa of interest.

- **aRborist Input:** a list of taxa (a list of genus/species names) and a few simple options (loci of interest, your NCBI API, etc.)  
- **aRborist Output:** curated metadata sheets (can be used to create phylogenetic trees or assign host to taxa, using downstream aRborist pipelines)


### Requirements

- R (≥ 4.2)
- RStudio (optional but recommended)
- Internet access to query NCBI
- (Recommended) an NCBI API key

### Disclaimer and Limitations

The accuracy and completeness of your results depend on the quality of metadata available in NCBI. While aRborist applies consistent naming, standardization, and error-handling routines, it cannot correct for missing, inconsistent, or ambiguous source data. I encourage users to review and, if necessary, manually refine the curation rules for your own use.

---

## aRborist first time setup

### 1) (Recommended) Set your NCBI API key 

This increases the NCBI allowed request rates. Get your API key by logging into your NCBI account, open Account settings, and scroll down to "API Key Management". You can then copy your API key.

Create (or edit) a file named ~/.Renviron and add:

```bash
NCBI_API_KEY=YOUR_KEY_HERE
```

Save this file and restart your R instance. Then in R:

```R
Sys.getenv("NCBI_API_KEY")  # should show your key (or at least not be empty)
```

### 2) Get the code

In R: 
```r
install.packages("usethis") # if needed
usethis::create_from_github("KScott6/aRborist")
```

This will automatically create an R project file "aRborist.Rproj". You can working in this project file.

### 3) Install/load the required packages

```R
source(file.path("R", "prepare.R"))  # installs any missing packages required by aRborist
```

The first time setup is now complete.

---

## aRborist walkthrough

### 1) Load packages and helper scripts (run once per new project)

```R
source(file.path("R", "arborist_helpers.R"))
load_required_packages()
```

### 2) Create a new project workspace (run once per new project)

Set the location you prefer to have your aRborist project folder in; ("~") is default. Then, select a unique name for your project. You can revisit individual projects later on by referencing their project name.

```R
projects_dir <- file.path(path.expand("~"), "aRborist_Projects")
project_name <- "Blackwellomyces_tree_2025_10_16"
```

Then run:

```R
start_project(projects_dir, project_name) 
```

This will automatically create the necessary project folders and subfolders. 

### 3) Set options for your project

Exaplaination of options:

`taxa_of_interest` Provide a list of genus or species names to search on NCBI. 

`organism_scope` Change to your target taxa's correct kingdom. You can use either NCBI taxonomy codes or simple names ("txid2[Organism:exp]" = "Bacteria[Organism:exp]"; "txid33090[Organism:exp]" = "Viridiplantae"). Fungi (txid4751) is the default option. Using this option helps reduce the number of incorrect hits during your search. Leave blank ("") to not have any organism contraints in your search (not recommended).

`search_options` Change to fit your search needs. Use the terms for the [NCBI advanced search](https://www.ncbi.nlm.nih.gov/nuccore/advanced). The default search string is shown.

`max_acc_per_taxa` Provive integer value to specify the maximum number of accessions to obtain for each taxon name. Use the option "max" to retrieve **all** the matching NCBI hits -- but be warned that for taxa with many accessions (Fusarium, Alternaria, etc.) this can make the metadata retreival step take **<u>a really long time</u>** (days). 

`ncbi_api_key` Your NCBI API key. If you didn't set up your R environment with your API key, you can specify it here in quotes.

`my_lab_sequences` Optional. If you want to include your own lab sequences, provide a 5-column csv with Accession, strain, sequence, organism, and gene columns (see [example_lab_seq_input.csv](example_data/example_lab_seq_input.csv) for an example).

```R
# Set options for this project
taxa_of_interest   <- c("Blackwellomyces", "Flavocillium")
organism_scope <- "txid4751[Organism:exp]"
search_options <- "(biomol_genomic[PROP] AND (100[SLEN]:5000[SLEN])) NOT Contig[All Fields] NOT scaffold[All Fields] NOT genome[All Fields]"
max_acc_per_taxa   <- 1000 # "max" for all matching hits
ncbi_api_key <- Sys.getenv("NCBI_API_KEY")
my_lab_sequences   <- "" 

# Save the exact options you used in your project folder (for reproducibility)
save_project_config(
  project_dir            = getwd(),
  project_name           = project_name,
  taxa_of_interest       = taxa_of_interest,          
  my_lab_sequences       = my_lab_sequences,
  organism_scope         = organism_scope,
  search_options         = search_options,
  max_acc_per_taxa       = max_acc_per_taxa    
)

```

### 4) Run aRborist to collect metadata

```R
ncbi_data_fetch(project_name, max_acc_per_taxa, taxa_of_interest)
```

<br> 

# aRborist Phylogenetic tree pipeline

### set options for filtering

`regions_to_include` Provide a list of loci of interest. Try to use the spelling/abbreviations commonly used on NCBI (e.g. if you search "TEF1" you'll get many more hits than if you searched "ef-1"). **should these match the output of the curation terms??**

`min_region_requirement` The number of user-specified loci an isolate must have in order to progress to be included in final downstream analyses. 

```R
regions_to_include <- c("RPB2", "TEF", "ITS")
min_region_requirement <- length(regions_to_include)
```

<br>

# aRborist host assessment pipeline

