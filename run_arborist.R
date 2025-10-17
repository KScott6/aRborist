# LOAD HELPER FUNCTIONS
source(file.path("R", "arborist_helpers.R"))

# USER OPTIONS (edit these for your project needs)
base_dir <- path.expand("~")
project_name <- "Blackwellomyces_tree"
taxa_of_interest <- c("Blackwellomyces", "Flavocillium")
regions_to_include <- c("RPB2", "TEF", "ITS")
max_acc_per_taxa <- "max"
min_region_requirement <- length(regions_to_include)
my_lab_sequences <- ""  # "" if none
ncbi_api_key <- Sys.getenv("NCBI_API_KEY")  # set once in ~/.Renviron

# RUN SETUP (run these in order)
setup_project_structure(base_dir)   # creates folders and setwd()
load_required_packages()

# NCBI data fetch
ncbi_data_fetch(project_name, max_acc_per_taxa, taxa_of_interest)

# Curation + FASTAs
data_curate(project_name, taxa_of_interest, NULL, min_region_requirement)