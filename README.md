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

This will automatically create an R project file "aRborist.Rproj". Continue working in this project file.

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

Name each of your projects something unique. You can revisit individual projects later on by referencing their project name.

```R
projects_dir <- file.path(path.expand("~"), "aRborist_Projects") # this creates a folder which will store all of your arborist projects
project_name <- "Blackwellomyces_tree_2025_10_16" # this is an examle project name. Set your project name to whatever you want.
start_project(projects_dir, project_name) # this will create a project folder in "aRborist_Projects", based on your provided project name
```

### 3) Set options for your project

Exaplaination of options:

`taxa_of_interest` Provide a list of genus or species names to search on NCBI. 

`organism_scope` Change to your organisms' correct kingdom. You can use NCBI taxonomy codes or simple names ("txid2[Organism:exp]" = "Bacteria[Organism:exp]"; "txid33090[Organism:exp]" = "Viridiplantae"). Fungi (txid4751) is the default option. Using this option helps reduce the number of incorrect hits during your search, but you could also leave this option blank ("") to not have any organism contraints in your search. 

`search_options` Change to fit your search needs. Use the terms for the [NCBI advanced search](https://www.ncbi.nlm.nih.gov/nuccore/advanced). 

`max_acc_per_taxa` Maximum number of accessions to obtain for each taxa/region search. Use the option "max" to retrieve **all** the matching NCBI hits -- but be warned that for taxa with many accessions (Fusarium, Alternaria, etc.) this can make the metadata retreival step take <u>**a really long time**<\u>. 

`ncbi_api_key` Your NCBI API key. If you didn't set up your R environment with your API key, you can do it here.

`my_lab_sequences` If you want to include your own lab sequences, provide a 5-column csv with Accession, strain, sequence, organism, and gene columns (see [example_lab_seq_input.csv](example_data/example_lab_seq_input.csv) for an example).

```R
# Set options for this project
taxa_of_interest   <- c("Blackwellomyces", "Flavocillium") # provide list of taxa
organism_scope <- "txid4751[Organism:exp]" # Set organism constraint. Provide NCBI taxonomy ID or common name. Fungi (txid4751) is the default option. Leave blank ("") to disregard organism constraint entirely (not recommended).
search_options <- "(biomol_genomic[PROP] AND (100[SLEN]:5000[SLEN])) NOT Contig[All Fields] NOT scaffold[All Fields] NOT genome[All Fields]" # default shown. Use NCBI advanced search terms. 
max_acc_per_taxa   <- 1000 # put "max" if you want all the matching hits per search/region search, otherwise indicate an integer.
ncbi_api_key <- Sys.getenv("NCBI_API_KEY")  # recommended to provide an API key. Provide your key here in quotes, if you didn't set it up in your environment previously.
my_lab_sequences   <- ""  # Optional. Leave blank or provide a path to a CSV with lab sequence info

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

`regions_to_include` Provide a list of loci of interest. Try to use the spelling/abbreviations commonly used on NCBI (e.g. if you search "TEF1" you'll get many more hits than if you searched "ef-1").

`min_region_requirement` The number of user-specified loci an isolate must have in order to progress to be included in final downstream analyses. 

```R
regions_to_include <- c("RPB2", "TEF", "ITS") # leave blank ("") if you don't want to filter your search by loci.
min_region_requirement <- length(regions_to_include) # leave as-is if you require all isolates to have all specified regions to be included in the final downstream analyses.
```

<br>

# aRborist host assessment pipeline

