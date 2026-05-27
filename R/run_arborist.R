############################
# Base aRborist pipeline - collect and curate NCBI nucleotide sequences and metadata
############################


setwd("~/github/aRborist") # change to your arborist download location
source(file.path("R", "arborist_helpers.R"))
load_required_packages()

arborist_home <- "/Users/scott/aRborist_Projects"
project_name <- "Dactylonectria_tree"

start_project(project_name = project_name)

# Set options for this project
taxa_of_interest <- c("Nectriaceae")

organism_scope <- "txid4751[Organism:exp]"


raw_entrez_terms <- list(
  Nectriaceae_unclassified = 'Nectriaceae sp.[porgn:__txid1756110] NOT uncultured[All Fields]'
)

taxa_of_interest <- c("Nectriaceae_unclassified")

# not using these in this run
search_include <- c(
  "biomol_genomic[PROP]",
  "(100[SLEN]:5000[SLEN])",
  "txid1756110[Organism:exp]"
)

# not using these in this run
search_exclude <- c(
  "Contig[All Fields]",
  "scaffold[All Fields]",
  "genome[All Fields]",
  "uncultured[All Fields]"
)

max_acc_per_taxa <- "max"   # use "max" to retrieve all matching hits
ncbi_api_key <- Sys.getenv("NCBI_API_KEY")
my_lab_sequences <- "/Users/scott/github/aRborist/example_data/dactylonectria_extracted_sequences.tsv"
literature_accessions <- "/Users/scott/github/aRborist/example_data/Dactylonectria_literature_acc.tsv"

# Save the exact options you used in your project folder
save_project_config(
  project_name = project_name,
  taxa_of_interest = taxa_of_interest,
  my_lab_sequences = my_lab_sequences,
  organism_scope = organism_scope,
  max_acc_per_taxa = max_acc_per_taxa
)

ncbi_data_fetch(
  taxa_list = taxa_of_interest,
  max_acc_per_taxa = max_acc_per_taxa,
  project_name = project_name
)

data_curate(project_name, taxa_of_interest = NULL, my_lab_sequences = my_lab_sequences)


############################
# Multi-gene tree aRborist pipeline
############################

merge_metadata_with_custom_file(project_name)

flag_literature_accessions(
  project_name = project_name,
  literature_accessions = literature_accessions
)

curate_metadata_regions(project_name)

regions_to_include <- c("His3", "ITS", "TEF", "BTUB","LSU")

select_regions(project_name,
               regions_to_include,
               acc_to_exclude = character(0),
               min_region_requirement = 3, 
               prefer_literature_accessions = TRUE,
               allow_compound_regions_for = c("ITS"))


## then manual editing of the attendance file at this point

run_dir <- start_phylogeny_run(
  project_name,
  regions_to_include,
  run_label = "filtered1"
)

create_multifastas(
  project_name,
  regions_to_include,
  run_dir = run_dir,
  use_tree_filter = TRUE
)

align_regions_mafft(
  project_name,
  regions_to_include,
  run_dir = run_dir,
  force = TRUE
)

trim_regions_trimal(
  project_name,
  regions_to_include,
  run_dir = run_dir,
  trimal_args = c(
    "-gt", "0.9",
  "-cons", "60",
  "-resoverlap", "0.8",
  "-seqoverlap", "75"),
  force = TRUE
)

  
write_final_region_attendance_sheet(
  project_name,
  regions_to_include,
  run_dir = run_dir
)

iqtree_modelfinder_per_region(
  project_name,
  regions_to_include,
  run_dir = run_dir,
  force = TRUE
)

concatenate_and_write_partitions(
  project_name,
  regions_to_include,
  run_dir = run_dir
)

iqtree_multigene_partitioned(
  project_name,
  regions_to_include,
  run_dir = run_dir,
  force = TRUE
)