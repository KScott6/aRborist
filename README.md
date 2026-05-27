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

<br>

### 2) Get the code

In R: 
```r
install.packages("usethis") # if needed
usethis::create_from_github("KScott6/aRborist")
```

This will automatically create an R project file "aRborist.Rproj". You can working in this project file.

<br>

### 3) Install/load the required packages

In R, change your working directory to where you have downloaded the aRborist scripts. 

For me, it was: `/Users/scott/github/aRborist`

```R
setwd("/Users/scott/github/aRborist") # change to your download location.

source(file.path("R", "prepare.R"))  # installs any missing packages required by aRborist
```

The first time setup is now complete.

---

<br>
<br>

## aRborist walkthrough

<br>

### 1) Load packages and helper scripts (run once per new project)

Start by moving into your local aRborist GitHub folder, sourcing the helper script, and loading the required R packages.

```R
setwd("~/github/aRborist") # change to your arborist download location
source(file.path("R", "arborist_helpers.R"))
load_required_packages()
```

<br>

### 2) Create or reopen a project workspace

Choose a unique project name. This name will be used to create the project folder and to name many of the output files.

```R
project_name <- "Blackwellomyces_tree_2025_10_16"
```

Then create the project workspace:

```R
start_project(project_name = project_name)
```

By default, this creates a project folder inside the aRborist project directory, for example: ~/github/aRborist/projects/Blackwellomyces_tree_2025_10_16

If the project already exists, start_project() can be run again to reload the project settings and continue working in the same project.

<br>

### 3) Set options for your project

Exaplaination of options:

`taxa_of_interest` Provide one or more genus/species names to search on NCBI, or provide a file that contains a list.

`organism_scope` Change to your target taxa's correct kingdom. This help reduce incorrect hits. Leave as "" to remove this restriction, though that is usually not recommended. Examples:

* Fungi: "txid4751[Organism:exp]"
* Bacteria: "txid2[Organism:exp]"
* Plants: "txid33090[Organism:exp]"

`search_include` A character vector of NCBI search filters that must be included in the search.

`search_exclude` A character vector of NCBI search filters that should be excluded from the search.

`raw_entrez_terms` Provide a string to be directly searched on NCBI. This option bypasses the automatic query construction and instead uses a manually supplied Entrez query exactly as written. If you specify anything for this option, all terms in `search_include` and `search_exclude` are ignored, so do not include these variables in your ncbi_data_fetch command.

   Normally, aRborist assumes each entry in taxa_of_interest is a plain organism name (e.g. "Fusarium"), and automatically builds an Entrez query using the helper filters and organism scope settings. This works well for most projects, but doesn't work so well when very specific NCBI search behavior is needed. You will need to pass a placeholder term in taxa_of_interest which will be used as the internal label used by aRborist for filenames and downstream grouping. If you use this option, you cannot use the optional taxa filtering in the inital curation step of the metadata curation; you must pass taxa_of_interest = NULL .

`max_acc_per_taxa` Provive integer value to specify the maximum number of accessions to obtain for each taxon name. Use the option "max" to retrieve **all** the matching NCBI hits -- but be warned that for taxa with many accessions (Fusarium, Alternaria, etc.) this can make the metadata retreival step take **<u>a really long time</u>** (days). 

`ncbi_api_key` Your NCBI API key. If you didn't set up your R environment with your API key, you can specify it here in quotes.

`my_lab_sequences` Optional. If you want to include your own lab sequences, provide a 5-column csv with Accession, strain, sequence, organism, and gene columns (see [example_lab_seq_input.csv](example_data/example_lab_seq_input.csv) for an example).

`literature_accessions` Optional, used in multi-gene tree pipeline. A table of literature-derived accessions that should be tracked and optionally prioritized during region selection. Table must contain columns named paper_id, region, and accession. Optional additional columns such as strain or organism names may also be included. (see [example_literature_acc.tsv](example_data/example_literature_acc.tsv) for an example).

   The provided literature accessions are flagged in the curated metadata, tracked in the attendance sheets, (optionally) prioritized during region selection, and any literature accessions not found in your metadata are recorded to the missing_literature_accessions_<project>.csv report.

```R
# Set options for this project
taxa_of_interest <- c("Blackwellomyces", "Flavocillium")

# or, read in a file like this:
# taxa_of_interest <- read_lines("/Users/$USER/Desktop/genera.txt")

organism_scope <- "txid4751[Organism:exp]"

search_include <- c(
  "biomol_genomic[PROP]",
  "(100[SLEN]:5000[SLEN])"
)

search_exclude <- c(
  "Contig[All Fields]",
  "scaffold[All Fields]",
  "genome[All Fields]"
)

# or, for raw string searching:

#raw_entrez_terms <- list(
#  Nectriaceae_unclassified = 'Nectriaceae sp.[porgn:__txid1756110] NOT uncultured[All Fields]'
#)

max_acc_per_taxa <- 1000   # use "max" to retrieve all matching hits
ncbi_api_key <- Sys.getenv("NCBI_API_KEY")
my_lab_sequences <- "/path/to/my_lab_sequences.tsv" # put "" if you have no lab sequecnes to add
literature_accessions <- "/path/to/literature_accessions.tsv" # put "" if you have no literature accesssions to flag

# Save the exact options you used in your project folder
save_project_config(
  project_name = project_name,
  taxa_of_interest = taxa_of_interest,
  my_lab_sequences = my_lab_sequences,
  organism_scope = organism_scope,
  max_acc_per_taxa = max_acc_per_taxa
)

```

<br>

### 4) Collect metadata

**Important:** This can be VERY time-intensive for large datasets.

(With my default parameters, I retrieved ~410,000 accessions and it took ~4 days to get all the metadata)

**Tip:** For very large runs, consider testing your pipeline on a small subset first (e.g., max_acc_per_taxa = 50) to confirm that your search parameters behave as expected before scaling up.

About NCBI search behavior:  standard NCBI searches are not perfectly constrained to the "organism" field. For example, if you search "Pandora[organism]", NCBI will return all accessions explictedly labeled as "Pandora" in the "organism" field, as well as any accessions that have "Pandora" located anywhere in the metadata (such as in the "notes" or "Title" field). It will also include any accession that was historically named "Pandora" as well, I believe.  This is frustrating, as it will slow down your search by including accessions you don't care about. I haven't found a foolproof way around this yet. I've tried a few workarounds (like constraining the search with "Pandora"[Organism:noexp]"), but this appears to still let a few unwanted accessions appear in the search. I've addressed this later on in the curation steps - there is a step that by default filters out any accession whose organism name doesn't match to your list of target taxa. As a result, you will probably have more accessions listed in your various intermediate files than you do in your final metadata file; this is normal and not a cause for concern.

Metadata retrieval is checkpointed by taxon (e.g., genus). Each taxon is written to its own file during the run (./metadata_files/metadata_checkpoints/metadata_<taxon>.csv). If the run is interrupted (e.g., laptop sleeps, internet drops, R crashes), progress is not lost! Re-running the same command with "resume = TRUE" will automatically skip metadata retreival for taxa that have already completed and will continue from where the run left off.

<br>

Running metadata retrieval:

```R
ncbi_data_fetch(
  taxa_list = taxa_of_interest,
  max_acc_per_taxa = max_acc_per_taxa,
  organism_scope = organism_scope,
  include_filters = search_include,
  exclude_filters = search_exclude,
  project_name = project_name
)
```

Resuming metadata collection (in case of metadata retreival interruption):

The ncbi_data_fetch() function runs both accession retrieval and metadata collection. If your metadata run is interrupted, you do not need to rerun everything. Instead, you can resume metadata collection directly:

```R
retrieve_ncbi_metadata(project_name, resume = TRUE)
```

<br>

### 5) Curation of metadata

Now that you have your metadata, it's time to do some basic curation. 

NOTE:  Public metadata is highly inconsistent and often incomplete. Its quality depends entirely on the original submitter. aRborist attempts to standardize common fields and naming patterns, but you should expect irregularities like missing values, inconsistent strain naming, or unusual formatting. Review your curated data before downstream analyses and keep these limitations in mind.

These are the curation steps that are peformed:

1) Integrate custom sequences (optional)
   If you provided a file via my_lab_sequences, your custom sequences and metadata are merged into the NCBI metadata before any curation steps. This allows your data to be treated identically to public accessions throughout the pipeline.

2) Assign a universal strain name
   Each accession receives a unified strain identifier (strain.standard) drawn from the following metadata fields, in order of priority:
specimen_voucher → strain → isolate → Accession. If none of these fields are available, the accession number is used.

1) Standardize strain names
   All spaces and special characters are removed to ensure compatibility with downstream analyses and FASTA headers. The standardized strain name is called "strain.standard".
   Example:  Both "ARSEF 1234" and "ARSEF-1234" become "ARSEF1234"

2) Flag accessions from type material. 
   If an accession is associated with type material (e.g., holotype, isotype, ex-type, etc.), "TYPE" is appended to the standardized strain name. Stored in: "strain.standard.type".
   Example: "ARSEF1234.TYPE"

3) Optional taxon filtering. 
   By default, any accession whose "organism" name does not match your specified taxa_of_interest is removed. You can disable this with taxa_of_interest = NULL.
   
4) Add strain/accession taxonomy
   aRborist automatically extracts full fungal taxonomy for each accession from NCBI metadata category **GBSeq_taxonomy**. THe following columns are added: Strain.taxonomy, Strain.phylum, Strain.class, Strain.order, Strain.family, Strain.genus, Strain.species

<br>

Running the basic curation:

Perform basic curation with taxon filtering (recommended):

```R
data_curate(project_name)
```


After this step completes, a new file will be created in your project directory:

./metadata_files/all_accessions_pulled_metadata_<project_name>_curated.csv


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

## Restarting a project

Remember - you can restart or jump between projects at any time by running the start_project command with the desired project name. If you just got done restarting your R instance to install the alignment and trimming software and you wanted to restart the test project, run these commands:

```R
# load up arborist packages
setwd("~/github/aRborist") # change to your arborist download location
source(file.path("R", "arborist_helpers.R"))
load_required_packages()

# specify the name of your project
arborist_home <- "/Users/scott/aRborist_Projects"
project_name <- "Blackwellomyces_tree_2025_10_16"

# then start your project again
start_project(project_name = project_name)
```

<br>

## 0) Optional: adding lab sequences and/or flagging accessions

If you are using custom lab sequences or literature accession flags, run those steps before region curation:

```R
merge_metadata_with_custom_file(project_name)

flag_literature_accessions(
  project_name = project_name,
  literature_accessions = literature_accessions
)
```

## 1) Standardizing region names

You need to further curate your metadata so that the gene region information is useable. The goal is to assign a consistent set of region identifiers for each accession, even when the original records use messy or compound descriptions.

(!) Before running this step, make sure you have run the basic curation step and have this file: ./metadata_files/all_accessions_pulled_metadata_<project_name>_curated.csv
   
This step uses a user-editable "replacement patterns" file (example_data/region_replacement_patterns.csv) to detect and standardize region names (e.g., ITS, TEF, RPB2, LSU, SSU). 

Each "pattern" is a regular expression (regex) that will be searched (case-insensitive) in the "gene", "product", and "accession_title" metadata text fields.

Each "standard" is the region name or label you want assigned when that pattern is found.

You can add as many fragments as you want to catch multi-region descriptions. Make sure to adjust for ambiguous fragments (e.g. specify matching to any "\bACT\b" instead of just "ACT", so phrases like "D**act**ylonectria beta-tubulin" aren't mis-labeled)

Each hit appends to the component list for that field. Any accessions with unmatched gene, product, AND acc_title categories will be logged in a separate file.

You can add as many lines as you want, the file acts as a flexible compound detector.

For example, if a gene description contains:

> "internal transcribed spacer 1; 5.8S ribosomal RNA; large subunit ribosomal RNA"

and your mapping file includes those three patterns, the resulting cell will record gene.components as:

>ITS;5.8S;LSU

aRborist will then perform the same pattern search for the "product" and "acc_title" categories for that accession. 

Then, aRborist combines the component fields to assign a final "region.standard" column using the priority: gene.region.components > product.region.components > acc_title.region.components

If you notice accessions in the new unmatched regions file (./metadata_files/unmatched_regions_Blackwellomyces_tree.csv), you can simply add the necessary pattern information to the replacement patterns file, save, and re-run the region curation step. It is not necessary for every accession to receive a region assignment. The practical goal is to capture the regions needed for the phylogeny.

<br>

To run the region curation step:

```R
curate_metadata_regions(project_name)
```

 <br> 

## 2) Filter metadata to desired regions

After region standardization, the next step is to select the loci that will be used for phylogenetic analysis and generate a filtered region attendance sheet.

This step creates a subfolder in /phyogenies named after the specified regions of interest. All downstream analyses will be stored in this subfolder, unless the regions of interest are changed. 

```R
regions_to_include <- c("ITS", "TEF")

select_regions(project_name,
               regions_to_include = regions_to_include,
               acc_to_exclude = character(0),
               min_region_requirement = 2)
```

`regions_to_include` Provide a vector of the standardized region names you want to include in the phylogeny. These terms will match the "standard" region names you specified in the region curation step. Otherwise, use standard NCBI region names such as ITS, TEF, RPB1, etc.

`min_region_requirement` Controls how many of the requested loci a strain must possess to remain in the downstream dataset.

`acc_to_exclude` Optional vector of specific accessions to remove before region selection. This is useful when duplicate or low-quality accessions exist for the same strain and region.

   example: acc_to_exclude <- c("PP464689", "PP464690")

`allow_compound_regions_for` Allows accessions to be assigned to "compound" regions.

   Some NCBI accessions contain multiple loci in a single sequence, for example: ITS;LSU;SSU. By default, aRborist allows compound-region accessions to satisfy searches for ITS (since it is so often submitted with partial LSU and SSU seqences). This means that an accession annotated as ITS;LSU;SSU will by default be selected when requesting ITS sequences. For all other loci, exact region matching is required unless explicitly added to allow_compound_regions_for. For example: allow_compound_regions_for = c("ITS", "LSU") would allow compound annotations containing LSU to be included during LSU selection.

`prefer_literature_accessions` Optional setting that prioritizes literature-derived accessions during duplicate resolution. When TRUE, accessions flagged in flag_literature_accessions() are preferentially selected if multiple accessions exist for the same strain and region. Only one accession per strain × region combination is retained in the attendance sheet. Manual edits to the attendance sheet always override automatic accession selection.

<br>

This step produces:

1) selected_accessions_metadata_<project>.<region_set>.csv : A long-format metadata file containing only the selected regions.

2) Region_attendance_sheet_<project>.<region_set>.csv : A wide-format attendance sheet showing one row per strain and one column per region.

3) region_selection_policy_<project>.<region_set>.txt : a record of the exact filtering settings used during selection.


(!) Important note about duplicates:  Public metadata is messy, and it’s common to have more than one accession for the same strain and the same region (for example, two ITS sequences submitted at different times). In this step, the script keeps only one accession per strain × region combination. 

So if you are trying to include a particular accession, but you find that a duplicate entry or entires keeps being used in place of your desired accession, you can specify to remove those particular accessions with the "acc_to_exclude" option, like so:

> acc_to_exclude = "PP46469,PP464690"

<br>

## 3) Creating a phylogeny run directory

After selecting the desired loci and generating the initial attendance sheet, the next step is to create a dedicated phylogeny run directory. This step creates a self-contained working directory for a specific phylogenetic analysis within the region set subfolder in /phylogenies.

Many phylogenetic projects involve multiple rounds of filtering, alignment trimming, manual editing, or exploratory analyses. Rather than overwriting earlier outputs, aRborist stores each phylogeny attempt in its own run folder.

You can and should modfiy the run_label and re-run this command and whenever you are modifying the input accessions for tree creation. This way, you can compare files/trees across your analyses.

```R
run_dir <- start_phylogeny_run(
  project_name,
  regions_to_include,
  run_label = "standard1"
)
```

`run_label` A short descriptive label for the phylogeny run. 

<br>

## 4) Optional (BUT HIGHLY RECOMMENDED): Manual editing of strains/accessions

Before generating multifastas and alignments, it is often useful to manually review and edit the input data. Chances are, there are strains or accessions that you don't especially want in any downstream analysis connected to this project. You can manually edit the accessions/strains that are included in future analyses by editing the primary attenance sheet. 

Any change made to this primary attendance sheet will be propagated to the subsequent run folders for this project. (i.e. if you remove an accession from the primary attendance sheet, any analyses created after this change will not have this accession.)

The primary attendance sheet is located here:

```text
~/github/aRborist/projects/<project_name>/phylogenies/<region_set>/Region_attendance_sheet_<project_name>.<region_set>.csv
```

This file controls which strains and accessions are included in the downstream phylogeny workflow.

Each row represents a strain, and each region column contains the accession selected for that locus.

You can control the sequences present in the downstream analyses by:

* disregard certain strains from being included (recommend action: see include_in_tree section)
* permanently removing problematic strains (delete entire row)
* permanently removing accessions known to be of poor quality or contaminant origin (replace select accession with NA)
* replacing accessions (edit accession)


#### include_in_tree : Temporary inclusion/exclusion

Sometimes you may want to temporarily remove strains from a particular analysis without permanently deleting them from the attendance sheet.

To do this, add a new column named: include_in_tree

Then, for any strain you DO want in your analysis, specify TRUE in that column.

If you created this column, when running the next step (create_multifastas) you need to specify use_tree_filter = TRUE

<br>

## 5) Create multifastas

To create multifastas that contain the RAW sequence data from NCBI, run this command:

```R
create_multifastas(
  project_name,
  regions_to_include,
  run_dir = run_dir,
  use_tree_filter = TRUE
)
```

Only specify use_tree_filter = TRUE if you actually created the include_in_tree column and want to use those selected strains. 

You will see a multifasta appear for each of the regions you specified.

<br>

## 6) Align each region

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

Note:  Before proceeding further, I recommened checking the alignments with a alignment GUI just to make sure there isn't any rouge sequneces messing up the alignment. If someone uploaded a TEF sequence but labeled it as ITS, this could really mess up your alignment and any downstream processes. 

<br>

## 7) Trim each region

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
   - Some of my favorite options: 
     - "-gt", "0.9", 
     - "-cons", "60", 
     - "-resoverlap", "0.8", 
     - "-seqoverlap", "75"

After this step, you will see a multifasta for the trimmed and aligned files (<region>.trimmed.fasta) as well as the log file from the mafft run (<region>.trimal.log).

<br>

## 8) Generate the final region attendance sheet 

After alignment and trimming, some sequences may be automatically removed during the trimming process. For example, highly incomplete, poorly aligned, or problematic sequences may no longer be present in the final trimmed FASTA files depending on which trimal parameters you used. So, a new attendance sheet must be produced to reflect the final set of accessions used.

To generate a final attendance sheet reflecting only the sequences that survived trimming, run:

```R
write_final_region_attendance_sheet(
  project_name,
  regions_to_include,
  run_dir = run_dir
)
```

This step compares:

* the intended attendance sheet,
* the aligned FASTA files,
* the trimmed FASTA files,

and determines which accessions successfully made it into the final trimmed alignments.

<br>

## 9) Create single-gene trees

Once you have trimmed alignments for each gene, the next step is to generate individual maximum-likelihood trees for each of your specified regions. If you are only interested in making a phylogeny from a single region, you can stop after this step as you will have your final tree.

If you are going to make a multi-gene tree, this step is still essential to identify the best substitution model for region region, as well as helping you find problematic loci, identify outliers, and confirm that sequences are behaving as expected before concatenation. 

This is how you create the single-gene trees with your trimmed alignments:

```R
iqtree_modelfinder_per_region(
  project_name,
  regions_to_include,
  threads = 8, # or whatever you like
  single_gene_bootstraps = 1000, # default bootstrap #          
  iqtree_args = c("-m", "MFP+MERGE"),   # MFP+MERGE necessary for model ID; you can add more IQ-TREE options here if needed
  run_dir = run_dir,
  force = TRUE
)
```

<br>

## 10) Create files necessary for multi-gene tree creation

The next step is to create the necessary files for the multi-gene tree in IQ-TREE. This step will create a concatenated supermatrix from the trimmed and aligned sequences, as well as a nexus (.nex) file that will store the sequence length and best substition model for each region. 

```R
concatenate_and_write_partitions(
  project_name,
  regions_to_include,
  run_dir = run_dir
)
```

<br>

## 11)  Create multi-gene tree with partitioned analysis

When creating phylogenies from multiple genes, I prefer to run a [partitioned analysis](https://iqtree.github.io/doc/Advanced-Tutorial) rather than use a single substituion model with the concatenated supermatrix. 

Each gene evolves under its own substitution dynamics- this means that the rates of evolution, base composition, among-site rate heterogeneity, and patterns of selective constraint can differ widely across loci. If you force a single substitution model onto one giant concatenated alignment, you assume all sites evolve exactly the same, an assumption that is almost always unrealistic in multigene datasets. In most cases, I find that running a paritioned analysis results in a tree with a structure that makes more sense and has greatly improved bootstrap support values.

This command will run a partitioned analysis in IQ-TREE:

```R
iqtree_multigene_partitioned(
  project_name,
  regions_to_include,
  run_dir = run_dir,
  threads = 8,
  multigene_bootstraps = 1000,
  iqtree_args = c("-redo"),
  force = TRUE
)
```

After this step, the multi-gene phylogeny pipeline is complete. You can find your final consensus tree file here: (./multi_gene_trees/iqtree_<project_name>.<regions>.contree). 

<br>
<br>
<br>

## Software citations

Make sure you cite all the software and methods used wrapped into this pipeline. 

MAFFT: Katoh, K. and Standley, D.M., 2013. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Molecular biology and evolution, 30(4), pp.772-780.

TrimAL: Capella-Gutiérrez, S., Silla-Martínez, J.M. and Gabaldón, T., 2009. trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. Bioinformatics, 25(15), pp.1972-1973.

IQ-TREE: Wong, T.K., Ly-Trong, N., Ren, H., Baños, H., Roger, A.J., Susko, E., Bielow, C., De Maio, N., Goldman, N., Hahn, M.W. and Huttley, G., 2025. IQ-TREE 3: Phylogenomic Inference Software using Complex Evolutionary Models.

Partitioned analysis: Chernomor, O., Von Haeseler, A. and Minh, B.Q., 2016. Terrace aware data structure for phylogenomic inference from supermatrices. Systematic biology, 65(6), pp.997-1008.

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