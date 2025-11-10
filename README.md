## Overview

aRborist is an automated sequence and metadata harvester designed to simplify the process of gathering and organizing sequence data and metadata from the **NCBI nucleotide database**. It retrieves accessions for specified taxa, extracts and standardizes metadata, and prepares sequences and metadata for downstream analyses. After using aRborist to pull metadata/sequence data, you can use other aRborist functions to help you make a phylogenetic tree. Or assign host information to your taxa of interest.

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
max_acc_per_taxa   <- 1000 # specify "max" to obtain all matching hits
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

### 4) Collect metadata

Another warning - if you set the options to collect a lot of accesssions, this step can take quite a while. 

A note - when you search a term on NCBI, it will sometimes return accessions you are not interested in. For example, if you search "Pandora[organism]", NCBI will return all accessions explictedly labeled as "Pandora" in the "organism" field, as well as any accessions that have "Pandora" located anywhere in the metadata (such as in the "notes" or "Title" field). It will also include any accession that was historically named "Pandora" as well, I believe.  This is frustrating, as it will slow down your search by including accessions you don't care about. I haven't found a way around this yet. I've tried a few workarounds (e.g. "Pandora"[Organism:noexp]), but either they don't work or are too strict and result in too few hits. Later on in the curation steps, there is a step that automatically filtes out any accession whose organism name doesn't match to your list of target taxa. This means that you will probably have more accessions listed in your various intermediate files than you do in your final metadata file.

If you've already set the values for "project_name", "max_acc_per_taxa", and "taxa_of_interest", you don't need to change the following command - just copy/paste/run as-is.

```R
ncbi_data_fetch(project_name, max_acc_per_taxa, taxa_of_interest)
```

<br>

### 5) Curation of metadata

Now that you have your metadata, it's time to do some basic curation. 

NOTE:  Public metadata is highly inconsistent and often incomplete. Its quality depends entirely on the original submitter. aRborist attempts to standardize common fields and naming patterns, but you should expect irregularities like missing values, inconsistent strain naming, or unusual formatting. Review your curated data before downstream analyses and keep these limitations in mind.

These are the curation steps that are peformed:

1) Assign a universal strain name
   Each accession receives a unified strain identifier (strain.standard) drawn from the following metadata fields, in order of priority:
specimen_voucher → strain → isolate → Accession. If there is no voucher, strain, or isolate data recorded, the accession number itself is used.

2) Standardize strain names
   All spaces and special characters are stripped. The standardized strain name is called "strain.standard".
   Example:  Both "ARSEF 1234" and "ARSEF-1234" become "ARSEF1234". 

3) Flag accessions from type material. 
   If an accession’s metadata indicates type material (e.g., holotype, isotype, ex-type, etc.), "TYPE" is appended to the standardized strain name.
   Example: "ARSEF1234.TYPE" in the column strain.standard.type.

4) Optional taxon filtering. 
   By default, any accession whose "organism" name does not match your specified taxa_of_interest is removed. You can disable this by passing taxa_of_interest = NULL.


Running the basic curation:

Perform basic curation with taxon filtering (recommended):

```R
curate_metadata_basic(project_name)
```

Perform basic curation without taxon filtering:

```R
curate_metadata_basic(project_name, taxa_of_interest = NULL)
```

After this step completes, a new file will be created in your project directory:

./metadata_files/all_accessions_pulled_metadata_<project_name>_curated_basic.csv

<br> 

### 6) Standardizing region names

Now that your basic metadata has been standardized, the next step is to curate gene region information. The goal is to assign a consistent set of region identifiers for each accession, even when the original records use messy or compound descriptions.

1) Before running this step, make sure you have run the basic curation step and have this file: ./metadata_files/all_accessions_pulled_metadata_<project_name>_curated_basic.csv
   
This step uses a user-editable "replacement patterns" file to detect and standardize region names (e.g., ITS, TEF, RPB2, LSU, SSU). You can add as many fragments as you want to catch multi-region descriptions. Each hit appends to the component list for that field. If you forget to include a pattern, aRborist will still flag common regions (ITS, LSU, SSU) automatically as a fallback. Any accessions with unmatched gene, product, AND acc_title categories will be logged in a separate file.

1) The region replacement patterns is file located here: ./aRborist/example_data/region_replacement_patterns.csv

Each "pattern" is a regular expression (regex) that will be searched (case-insensitive) in the "gene", "product", and "acc_title" metadata text fields.
Each "standard" is the region name or label you want assigned when that pattern is found.

You can add as many lines as you want, the file acts as a flexible compound detector.

For example, if a gene description contains:

> "internal transcribed spacer 1; 5.8S ribosomal RNA; large subunit ribosomal RNA"

and your mapping file includes those three patterns, the resulting cell will record gene.components as:

>ITS;5.8S;LSU

aRborist will then perform the same pattern search for the "product" and "acc_title" categories for that accession. 

Then, aRborist combines the component fields to assign a final "region.standard" column using the priority: gene.region.components > product.region.components > acc_title.region.components

If you notice accessions in the new unmatched regions file (./metadata_files/unmatched_regions_Blackwellomyces_tree.csv), you can simply add the necessary pattern information to the replacement patterns file, save, and re-run the region curation step. It's not essential that every accession be assinged a region - you can stop whenever you feel you have captured all the useful information from the accessions you care about. 

To run the region curation step:

```R
curate_metadata_regions(project_name)
```

<br>

### End of basic aRborist pipeline

At this point, you should have a massive metadata file with curated information. Hopefully, this is in a format that is useful to you! 

I have several downstream pipelines that directly build off the output from these inital steps. These include:

1) Phylogenetic tree pipeline
   
   Semi-automated pipeline to build phylogeneies from one or more genes.

2) Host assessment pipeline

   Automatically parses through massive amounts of public data to assign host percentage at different taxonomic levels.

I create these pipelines primarily for myself, as needed for different projects, so I am always adding new offshoots of the aRborist pipeline, so this set of pipelines may expand in the future.

<br>

# aRborist Phylogenetic tree pipeline (includes fasta generation)

### set options for filtering

`regions_to_include` Provide a list of loci of interest. Try to use the spelling/abbreviations commonly used on NCBI (e.g. if you search "TEF1" you'll get many more hits than if you searched "ef-1"). **should these match the output of the curation terms??**

`min_region_requirement` The number of user-specified loci an isolate must have in order to progress to be included in final downstream analyses. 

```R
regions_to_include <- c("RPB2", "TEF", "ITS")
min_region_requirement <- length(regions_to_include)
```

<br>

# aRborist host assessment pipeline