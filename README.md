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
<br>

### End of basic aRborist pipeline

At this point, you should have a massive metadata file with curated information. Hopefully, this is in a format that is useful to you! 

I have several downstream pipelines that directly build off the output from these inital steps. These include:

1) Phylogenetic tree pipeline
   
   Semi-automated pipeline to build phylogeneies from one or more genes.

2) Host assessment pipeline

   Semi-automatically parses through massive amounts of public data to assign host percentage at different taxonomic levels.

I create these pipelines primarily for myself, making them as needed for different projects. I am always adding new offshoots of the aRborist pipeline, so this set of pipelines may expand in the future.

<br>
<br>

# aRborist Phylogenetic tree pipeline (includes fasta generation)

## Setup and software download

This pipeline is optional and is fully independent of the host assessment pipeline. You can run either or both of these pipelines, in any order, after completing the basic metadata curation step of the basic arborist pipeline.

You will need to download some external software before you proceed:  MAFFT (for alignment), TrimAl (for sequence trimming), and IQ-TREE (to actually create the phylogenetic trees). You will also need a software to view the phylogenies, such as [FigTree](https://github.com/rambaut/figtree/releases) or [TreeViewer](https://treeviewer.org/).

If you have conda installed on your computer, you can easily install the software:

> conda install -c bioconda trimal mafft iqtree

Or, you can manually install at their respective websites:  [TrimAl](https://vicfero.git)  [MAFFT](https://mafft.cbrc.jp/alignment/software/source.html) [IQ-TREE](https://iqtree.github.io/)

Take note of the full paths of your downloaded software. 

Open your .Renviron file:

```R
usethis::edit_r_environ()
```

Edit this file to include the lines:

```R
>MAFFT_PATH=<path_to_mafft_install>
# my install path was: /Users/scott/miniconda3/bin/mafft
>TRIMAL_PATH=<path_to_trimal_install>
# my install path was: /Users/scott/miniconda3/bin/trimal
>IQTREE_PATH=<path_to_iqtree_install>
# my install path was: /Users/scott/miniconda3/bin/iqtree
```

Save, and restart your R instance. 

If you don't want to set the paths permanently in the .Renviron file, you can just set the paths each time you open R.

To test if you have properly installed the software and R can find the binaries, run:

```R
Sys.getenv("MAFFT_PATH")
Sys.getenv("TRIMAL_PATH")
Sys.getenv("IQTREE_PATH")
```

If you see the full paths you just set, you are good to go. 

<br>

## 1) Standardizing region names

You need to further curate your metadata so that the gene region information is useable. The goal is to assign a consistent set of region identifiers for each accession, even when the original records use messy or compound descriptions.

a) Before running this step, make sure you have run the basic curation step and have this file: ./metadata_files/all_accessions_pulled_metadata_<project_name>_curated.csv
   
This step uses a user-editable "replacement patterns" file to detect and standardize region names (e.g., ITS, TEF, RPB2, LSU, SSU). You can add as many fragments as you want to catch multi-region descriptions. Each hit appends to the component list for that field. If you forget to include a pattern, aRborist will still flag common regions (ITS, LSU, SSU) automatically as a fallback. Any accessions with unmatched gene, product, AND acc_title categories will be logged in a separate file.

b) The region replacement patterns is file located here: ./aRborist/example_data/region_replacement_patterns.csv

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

## 2) Filter metadata to desired regions

Now that you have a more curated metadata file for the project, the next step is to narrow it down to just the accessions that will be used to build a tree. In aRborist, we do this by telling the pipeline which marker(s) we want to use (e.g. ITS, TEF, RPB2), and the script will pull out only the accessions that match those markers. This produces a clean, region-specific dataset that the alignment/trim steps can use later.

Remember - you can restart or jump between projects at any time by running the start_project command with the desired project name. If you just got done restarting your R instance to install the alignment and trimming software and you wanted to restart the test project, run these commands:

```R
# load up arborist packages
arborist_repo <- normalizePath("~/github/aRborist")
source(file.path(arborist_repo,"R", "arborist_helpers.R"))
load_required_packages()

# restart with the name of your project
start_project(project_name = "Blackwellomyces_tree")
```

Then, to filter the metadata, 

```R
select_regions(project_name,
               regions_to_include = c("ITS", "TEF"),
               acc_to_exclude = character(0),
               min_region_requirement = 2)
```

`regions_to_include` Provide a list of loci of interest. These terms will match the "standard" region names you specified in the region curation step. Otherwise, use standard NCBI region names such as ITS, TEF, RPB1, etc.

`min_region_requirement` The number of user-specified loci an isolate must have in order to progress to be included in final downstream analyses. 

What this step does:

1) Filters the metadata to only consider accessions of regions that you are interested in.

2) Filters the data to only include strain with X number of regions (you set the cutoff).

If you set

> min_region_requirement <- length(regions_to_include)

then only strains that have ALL of the regions you asked for will be kept (strict mode).

And if you set

> min_region_requirement <- 1

then any strain that has at least one of the regions will be kept.

Important note about duplicates:  Public metadata is messy, and it’s common to have more than one accession for the same strain and the same region (for example, two ITS sequences submitted at different times). In this step, the script keeps only one accession per strain × region combination. 

So if you are trying to include a particular accession, but you find that a duplicate entry or entires keeps being used in place of your desired accession, you can specify to remove those particular accessions with the "acc_to_exclude" option, like so:

> acc_to_exclude = "PP46469,PP464690"

<br>

## 3) Create multifastas

To create multifastas that contain the RAW sequence data from NCBI, run this command:

```R
create_multifastas(project_name, regions_to_include)
```

You will see a multifasta appear for each of the regions you specified. For example: "./Blackwellomyces_tree/phylogenies/ITS.RPB2.TEF/prep/ITS/Blackwellomyces_tree.ITS.RPB2.TEF.ITS.raw.fasta"


## 4) Align each region

Once the raw multifasta files are created for each region, the next step is to align the sequences. This step is carried out separately for each region you specified. By default, aRborist uses the alignment software MAFFT (although I may add more alignment software options in the future).

This step is run with:

```R
align_regions_mafft(project_name,
                    regions_to_include,
                    threads = max(1, parallel::detectCores() - 1),
                    mafft_args = c("--auto", "--reorder"),
                    force = TRUE)
```

Important parameters:

`threads` : how many CPUs MAFFT will use. Defaults to all but one available core. 

`extra_args` : allows you to pass different MAFFT parameters. For example:
   - "auto" : automatically selects the best algorithm based on the number and length of sequences
   - "--reorder" : lets MAFFT rearrange sequences internally to speed up the alignment
   - check out the MAFFT manual for more options

`force` : if TRUE, will overwrite preexisting alignment files in the project folder

After this step, each region folder will contain the aligned multifastas (<region>.aligned.fasta) as well as the log file from the mafft run (<region>.mafft.log). The aligned files are ready for trimming. 

Note:  Before proceeding further, I recommened checking the alignments with a alignment GUI just to make sure there isn't any rouge sequneces messing up the alignment. If someone uploaded a TEF sequence but labeled it as TEF, this could really mess up your alignment and any downstream processes. 

<br>

## 5) Trim each region

After alignment, many columns in the alignment may contain mostly gaps or poorly aligned positions. We also need to ensure that all the sequences for a particular region are the same length. aRborist uses TrimAl to perform these steps. 

The step is run with:

```R
trim_regions_trimal(project_name,
                    regions_to_include,
                    trimal_args = c("-automated1"))
```

Important parameters:

`trimal_args` : allows you to pass different TrimAl parameters. For example:
   - "-automated1" : trimAl automatically select thresholds for maximum allowed gap percentage per column, minimum overlap between sequences, and conservation scores.
   - check the TrimAl manual for more options

After this step, you will see a multifasta for the trimmed and aligned files (<region>.trimmed.fasta) as well as the log file from the mafft run (<region>.trimal.log).

<br>

## 6) Create single-gene trees

Once you have trimmed alignments for each gene, the next step is to generate individual maximum-likelihood trees for each of your specified regions. If you are only interested in making a phylogeny from a single region, you can stop after this step as you will have your final tree (./single_gene_trees/<region>/<project>.<region>.modeltest.contree). 

If you are going to make a multi-gene tree, this step is still essential to identify the best substitution model for region region, as well as helping you find problematic loci, identify outliers, and confirm that sequences are behaving as expected before concatenation. 

This is how you create the single-gene trees with your trimmed alignments:

```R
iqtree_modelfinder_per_region(
  project_name,
  regions_to_include,
  threads = 8, # or whatever you like
  single_gene_bootstraps = 1000, # default bootstrap #          
  iqtree_args = c("-m", "MFP+MERGE"),   # MFP+MERGE necessary for model ID; you can add more IQ-TREE options here if needed
  force = TRUE
)
```

<br>

## 7) Create files necessary for multi-gene tree creation

The next step is to create the necessary files for the multi-gene tree in IQ-TREE. This step will create a concatenated supermatrix from the trimmed and aligned sequences, as well as a nexus (.nex) file that will store the sequence length and best substition model for each region. 

```R
concatenate_and_write_partitions(project_name, regions_to_include)
```

<br>

## 8)  Create multi-gene tree with partitioned analysis

When creating phylogenies from multiple genes, I prefer to run a [partitioned analysis](https://iqtree.github.io/doc/Advanced-Tutorial) rather than use a single substituion model with the concatenated supermatrix. 

Each gene evolves under its own substitution dynamics- this means that the rates of evolution, base composition, among-site rate heterogeneity, and patterns of selective constraint can differ widely across loci. If you force a single substitution model onto one giant concatenated alignment, you assume all sites evolve exactly the same, an assumption that is almost always unrealistic in multigene datasets. In most cases, I find that running a paritioned analysis results in a tree with a structure that makes more sense and has greatly improved bootstrap support values.

This command will run a partitioned analysis in IQ-TREE:

```R
iqtree_multigene_partitioned(
  project_name,
  regions_to_include,
  threads = 8,
  multigene_bootstraps = 1000,
  iqtree_args = c("-redo"),
  force = TRUE
)
```

After this step, the multi-gene phylogeny pipeline is complete. You can find your final consensus tree file here: (./multi_gene_trees/iqtree_<project_name>.<regions>.contree). 

<br>

## Software citations

Make sure you cite all the software and methods used wrapped into this pipeline. 

MAFFT:

TrimAL:

IQ-TREE:

Partitioned analysis: 

<br>
<br>

# aRborist Host Assessment Pipeline

Taxonomic curation and summary of host associations for each species included in your metadata.

This host assessment pipeline takes your curated metadata from the basic arborist pipeline, cleans and standardizes it, looks up the host taxonomy via NCBI, and provides a helpful summary. This pipeline is optional and is fully independent of the phylogenetic pipeline. You can run either or both of these pipelines, in any order, after completing the basic metadata curation step of the basic arborist pipeline. 

All the output from this pipeline will be stored inside your project folder in a folder called "host_assessment".

<br> 

## 1) Initial host term extraction and taxonomy lookup

This step will create a new column in your metadata (host.standardized) and use the NCBI taxonomy database to look up the full taxonomy for each unique term.

At the end of the lookup process, you will be told how many host names failed the search. If by some miracle you have zero failed names, or if you don't care about using as much of the metadata as possible, you may proceed directly to step 3. 

```R
run_host_assessment_initial_pass(
  project_name,
  use_isolation_source = FALSE, 
  overwrite_host_standardized = TRUE
)
```

Explanation of options:

`use_isolation_source` : If the "host" metadata field is empty for an accession, will instead use the entry for "isolation_source".

Sometimes, when looking at metadata, it's really obvious that someone put down host infomation in the "isolation_source" category rather than the correct "host" category. I made this option in case I wanted to wring every bit of somewhat applicable information out of the metadata. Turning this option on will drastically increase the number of terms you need to search and edit, plus, chances are some of the isolation_source data truly is inappropriate to be considered as host data. Overall, I would recommend against using this option. 

`overwrite_host_standardized` : if TRUE, will overwrite the host_standarized column in your metadata file. Turn this on if you want to start the host assessment pipeline from scratch and need to re-do this step.

<br> 

## 2) Curation and re-attempt to lookup host terms

You probably had at least a few terms fail the NCBI taxonomy lookup. Metadata will often contain messy, misspelled, or ambiguous terms - this does not play well with automated searching. *If* you want to rescue as much metadata as possible, you'll need to do some manual editing.

The previous step has created a file listing all the failed terms:  ./host_assessment/failed_host_terms_<project>.csv

It contains two columns: "original_term" and "replacement_term". Open the file and change the values in the "replacement_term" column to a more approprate term. If there is a term you know you don't care about or want to skip (e.g. "soil", "culture from", nonsense) leave it as "NA". 

Some examples:

* "on insect cocoon buried in soil" ->  "Insecta"

* "Crinipellis pernikiosa" -> "Crinipellis perniciosa"

* "lepidopteran larva" -> "Lepidoptera"


Some notes:  I recommend being as conservative as possible when changing terms. Double check commonly mispelled names, or if an organism has more than one name. Also, if a name is not in the NCBI taxonomy database, it will not return the taxonomy (I have run into this problem a lot with esoteric plant taxa).


Once you have made the necessary edits, save the file. Then, run the following to re-attempt the taxonomy lookup with just the failed terms:

```R
run_host_assessment_refinement_pass(project_name)
```

You can re-run the refinement step as many times as needed. Once you feel good about the state of your data, move onto the next step.

<br>

## 3) Summarize host data

Now it's time to add your host taxonomy data back to your master datasheet (./metadata_files/all_accessions_pulled_metadata_<project_name>_curated.csv) and summarize the info so it's in a useable form. 

```R
run_host_assessment_summary(
  project_name,
  host_rank = "phylum",   # or "order", "class", "family", "genus", etc.
  keep_NAs  = FALSE       # include/exclude NA hosts from percentages
)
```

Explanation of options: 

`host_rank` :  specify the taxonomy rank you want to investigate. 

Many accessions only have high-level host information available, so it is best to begin with broader ranks (phylum, order) and only move to finer levels if the dataset supports it.

`keep_NAs` : control how missing or unusable host information affects percentage calculations.

If TRUE, aRborist will summarize host usage only among accessions with known host information. If FALSE, aRborist will incorporate NAs into the calculations; percentages become more conservative and may be strongly diluted by missing data. Not including the NAs will better reflect the overall data completeness, but will probably weaken the ecological signal because unknown hosts will probably dominate the totals.

<br>

With this step complete, the host assessment pipeline is done. You will have a final report of the host assessment of your metadata here: ./host_assessment/host_usage_by_<taxon_level>_<project_name>.csv

This file has one row for each species in your dataset, with the host breakdown at the specified taxon level. The "top_host_category" column reports the host taxa with the highest percetage per target speices. The "host_profile" column summarizes the host breakdown. 