# Base aRborist pipeline - collect and curate NCBI nucleotide sequences and metadata

## LOAD HELPER FUNCTIONS
arborist_repo <- normalizePath("~/github/aRborist")
source(file.path(arborist_repo,"R", "arborist_helpers.R"))

## USER OPTIONS (edit these for your project needs)
base_dir <- path.expand("~")
project_name <- "Blackwellomyces_tree"
taxa_of_interest <- c("Blackwellomyces", "Flavocillium")
regions_to_include <- c("RPB2", "TEF", "ITS")
max_acc_per_taxa <- "max"
my_lab_sequences <- ""  # path to input file, put "" if none
ncbi_api_key <- Sys.getenv("NCBI_API_KEY")  # or just set once in ~/.Renviron
acc_to_exclude <- "" # include any accessions you know you don't want to include
min_region_requirement <- length(regions_to_include) # default

## aRborist setup
load_required_packages()
start_project(project_name)

## NCBI data fetch
ncbi_data_fetch(
  taxa_list        = taxa_of_interest,
  max_acc_per_taxa = max_acc_per_taxa,
  organism_scope   = "txid4751[Organism:exp]",
  include_filters  = c("is_nuccore[filter]", "(00000000010[SLEN] : 00000020000[SLEN])"),
  exclude_filters  = c("mitochondrion[filter]")
)

## Basic metadata curation
curate_metadata_basic(project_name)

## Region curation
curate_metadata_regions(project_name)


# Phlyogeny pipeline

## Filter metadata to desired regions
select_regions(project_name, regions_to_include, acc_to_exclude, min_region_requirement)

## create raw multifastas
create_multifastas(project_name, regions_to_include)

## automatic alignment with mafft
align_regions_mafft(project_name,
                    regions_to_include,
                    threads = max(1, parallel::detectCores() - 1),
                    mafft_args = c("--auto", "--reorder"),
                    force = TRUE)

## automatic trimming with trimal
trim_regions_trimal(project_name,
                    regions_to_include,
                    trimal_args = c("-automated1")) 

## single-gene tree(s) and prep for multi-gene tree, if running
iqtree_modelfinder_per_region(
  project_name,
  regions_to_include,
  threads = 8, # or whatever you like
  single_gene_bootstraps = 1000, # default bootstrap #          
  iqtree_args = c("-m", "MFP+MERGE"),   # you can add more IQ-TREE options here if needed
  force = TRUE
)

## create files for multigene tree
concatenate_and_write_partitions(project_name, regions_to_include)

## multigene tree
iqtree_multigene_partitioned(
  project_name,
  regions_to_include,
  threads = 8,
  multigene_bootstraps = 1000,
  iqtree_args = c("-redo"),
  force = TRUE
)
