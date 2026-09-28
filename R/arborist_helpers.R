# ============================================================
# aRborist
# ============================================================

# ============================================================
# aRborist helper functions and defaults
# ============================================================

# Package setup
required_packages <- c(
  "rentrez","stringr","plyr","dplyr","withr","XML",
  "data.table","tidyr","phylotools","scales",
  "purrr","readr","phytools","RColorBrewer",
  "maps","ggplot2","tidygeocoder",
  "ggrepel","taxize","Biostrings", "yaml","readr"
)

installed_packages <- required_packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(required_packages[!installed_packages])
}

load_required_packages <- function() {
  invisible(lapply(required_packages, library, character.only = TRUE))
}


# ============================================================
# Build local NCBI taxonomy databases for host assessment
# ============================================================

setup_ncbi_taxonomy_database <- function(
    taxonomy_dir,
    output_dir = NULL,
    overwrite = FALSE
) {
  
  # ==========================================================
  # Check dependencies
  # ==========================================================
  
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop(
      "Package 'data.table' is required to build the local ",
      "NCBI taxonomy database.\n",
      "Install it with:\n",
      "  install.packages(\"data.table\")"
    )
  }
  
  
  # ==========================================================
  # Validate input directory
  # ==========================================================
  
  taxonomy_dir <- normalizePath(
    taxonomy_dir,
    mustWork = FALSE
  )
  
  if (!dir.exists(taxonomy_dir)) {
    stop(
      "NCBI taxonomy directory does not exist:\n  ",
      taxonomy_dir
    )
  }
  
  
  # ==========================================================
  # Determine output directory
  # ==========================================================
  
  if (is.null(output_dir)) {
    
    if (exists(
      ".host_example_data_dir",
      mode = "function"
    )) {
      
      output_dir <- .host_example_data_dir()
      
    } else {
      
      arborist_root <- get0(
        "arborist_repo",
        envir = .GlobalEnv,
        ifnotfound = normalizePath(
          "~/github/aRborist",
          mustWork = FALSE
        )
      )
      
      output_dir <- file.path(
        arborist_root,
        "example_data"
      )
    }
  }
  
  output_dir <- normalizePath(
    output_dir,
    mustWork = FALSE
  )
  
  if (!dir.exists(output_dir)) {
    dir.create(
      output_dir,
      recursive = TRUE
    )
  }
  
  
  # ==========================================================
  # Locate required NCBI taxonomy files
  # ==========================================================
  
  find_taxonomy_file <- function(filename) {
    
    hits <- list.files(
      taxonomy_dir,
      pattern = paste0(
        "^",
        gsub(
          "\\.",
          "\\\\.",
          filename
        ),
        "$"
      ),
      recursive = TRUE,
      full.names = TRUE
    )
    
    if (length(hits) == 0) {
      stop(
        "Required NCBI taxonomy file not found:\n  ",
        filename,
        "\n\nSearched within:\n  ",
        taxonomy_dir
      )
    }
    
    if (length(hits) > 1) {
      warning(
        "Multiple copies of ",
        filename,
        " were found. Using:\n  ",
        hits[[1]]
      )
    }
    
    hits[[1]]
  }
  
  
  names_file <- find_taxonomy_file(
    "names.dmp"
  )
  
  nodes_file <- find_taxonomy_file(
    "nodes.dmp"
  )
  
  rankedlineage_file <- find_taxonomy_file(
    "rankedlineage.dmp"
  )
  
  
  message("")
  message("NCBI taxonomy source files:")
  message("  names.dmp:         ", names_file)
  message("  nodes.dmp:         ", nodes_file)
  message("  rankedlineage.dmp: ", rankedlineage_file)
  
  
  # ==========================================================
  # Define output files
  # ==========================================================
  
  name_lookup_output <- file.path(
    output_dir,
    "ncbi_name_lookup.rds"
  )
  
  ranked_taxonomy_output <- file.path(
    output_dir,
    "ncbi_ranked_taxonomy.rds"
  )
  
  
  existing_outputs <- c(
    name_lookup_output,
    ranked_taxonomy_output
  )
  
  if (
    any(file.exists(existing_outputs)) &&
    !overwrite
  ) {
    
    existing_outputs <- existing_outputs[
      file.exists(existing_outputs)
    ]
    
    stop(
      "One or more local NCBI taxonomy database files already exist:\n  ",
      paste(
        existing_outputs,
        collapse = "\n  "
      ),
      "\n\nUse overwrite = TRUE to rebuild them."
    )
  }
  
  
# Helper for cleaning NCBI dump fields
clean_ncbi_field <- function(x) {
    
    x <- as.character(x)
    
    # Remove whitespace introduced by the NCBI .dmp format
    x <- trimws(x)
    
    # Treat empty values as NA
    x[
      is.na(x) |
        x == ""
    ] <- NA_character_
    
    x
  }
  
  # ==========================================================
  # Build NCBI name lookup database
  #
  # names.dmp fields:
  #
  # TaxID | name | unique_name | name_class |
  #
  # Keep ALL name classes. This allows aRborist to recognize
  # scientific names, synonyms, common names, misspellings,
  # equivalent names, etc. Ambiguous names are handled later
  # by the host-assessment pipeline.
  # ==========================================================
  
  message("")
  message("Reading NCBI names.dmp...")
  
  ncbi_names <- data.table::fread(
    names_file,
    sep = "|",
    header = FALSE,
    quote = "",
    fill = TRUE,
    data.table = FALSE,
    showProgress = TRUE
  )
  
  if (ncol(ncbi_names) < 4) {
    stop(
      "names.dmp did not contain the expected four fields."
    )
  }
  
  ncbi_names <- ncbi_names[
    ,
    1:4,
    drop = FALSE
  ]
  
  names(ncbi_names) <- c(
    "TaxID",
    "name",
    "unique_name",
    "name_class"
  )
  
  ncbi_names$TaxID <- clean_ncbi_field(
    ncbi_names$TaxID
  )
  
  ncbi_names$name <- clean_ncbi_field(
    ncbi_names$name
  )
  
  ncbi_names$unique_name <- clean_ncbi_field(
    ncbi_names$unique_name
  )
  
  ncbi_names$name_class <- clean_ncbi_field(
    ncbi_names$name_class
  )
  
  
  # ----------------------------------------------------------
  # Normalized version used for fast host-name matching
  # ----------------------------------------------------------
  
  if (exists(
    ".host_normalize_name",
    mode = "function"
  )) {
    
    ncbi_names$normalized_name <-
      .host_normalize_name(
        ncbi_names$name
      )
    
  } else {
    
    normalized_name <- trimws(
      as.character(
        ncbi_names$name
      )
    )
    
    normalized_name <- gsub(
      "[[:space:]]+",
      " ",
      normalized_name
    )
    
    ncbi_names$normalized_name <-
      tolower(
        normalized_name
      )
  }
  
  
  # Remove malformed rows, but retain all legitimate
  # NCBI name classes.
  
  ncbi_names <- ncbi_names[
    !is.na(ncbi_names$TaxID) &
      !is.na(ncbi_names$name) &
      !is.na(ncbi_names$name_class) &
      !is.na(ncbi_names$normalized_name) &
      nzchar(ncbi_names$normalized_name),
    ,
    drop = FALSE
  ]
  
  
  ncbi_names <- unique(
    ncbi_names
  )
  
  
  message(
    "NCBI name records retained: ",
    format(
      nrow(ncbi_names),
      big.mark = ","
    )
  )
  
  
  # ==========================================================
  # Read rankedlineage.dmp
  #
  # rankedlineage.dmp fields:
  #
  # TaxID
  # Scientific.name
  # Species
  # Genus
  # Family
  # Order
  # Class
  # Phylum
  # Kingdom
  # Superkingdom
  # ==========================================================
  
  message("")
  message("Reading NCBI rankedlineage.dmp...")
  
  ranked_taxonomy <- data.table::fread(
    rankedlineage_file,
    sep = "|",
    header = FALSE,
    quote = "",
    fill = TRUE,
    data.table = FALSE,
    showProgress = TRUE
  )
  
  if (ncol(ranked_taxonomy) < 10) {
    stop(
      "rankedlineage.dmp did not contain the expected ",
      "ten taxonomy fields."
    )
  }
  
  ranked_taxonomy <- ranked_taxonomy[
    ,
    1:10,
    drop = FALSE
  ]
  
  names(ranked_taxonomy) <- c(
    "TaxID",
    "Scientific.name",
    "Species",
    "Genus",
    "Family",
    "Order",
    "Class",
    "Phylum",
    "Kingdom",
    "Superkingdom"
  )
  
  
  for (column_name in names(ranked_taxonomy)) {
    
    ranked_taxonomy[[column_name]] <-
      clean_ncbi_field(
        ranked_taxonomy[[column_name]]
      )
  }
  
  
  # ==========================================================
  # Read nodes.dmp
  #
  # I use nodes.dmp to retain the ACTUAL rank of each TaxID.
  # rankedlineage.dmp provides the standard lineage columns,
  # but does not itself tell us whether the queried TaxID is
  # a species, subspecies, strain, genus, etc.
  # ==========================================================
  
  message("")
  message("Reading NCBI nodes.dmp...")
  
  ncbi_nodes <- data.table::fread(
    nodes_file,
    sep = "|",
    header = FALSE,
    quote = "",
    fill = TRUE,
    data.table = FALSE,
    showProgress = TRUE,
    select = 1:3
  )
  
  if (ncol(ncbi_nodes) < 3) {
    stop(
      "nodes.dmp did not contain the expected first ",
      "three fields."
    )
  }
  
  names(ncbi_nodes) <- c(
    "TaxID",
    "parent_TaxID",
    "taxon_rank"
  )
  
  ncbi_nodes$TaxID <- clean_ncbi_field(
    ncbi_nodes$TaxID
  )
  
  ncbi_nodes$parent_TaxID <- clean_ncbi_field(
    ncbi_nodes$parent_TaxID
  )
  
  ncbi_nodes$taxon_rank <- clean_ncbi_field(
    ncbi_nodes$taxon_rank
  )
  
  
  # Only TaxID and rank are required in the final lookup.
  
  rank_lookup <- ncbi_nodes[
    ,
    c(
      "TaxID",
      "taxon_rank"
    ),
    drop = FALSE
  ]
  
  
  rank_lookup <- rank_lookup[
    !duplicated(rank_lookup$TaxID),
    ,
    drop = FALSE
  ]
  

  # Add actual NCBI rank to ranked taxonomy
  message("")
  message("Adding NCBI taxon ranks...")
  
  ranked_taxonomy <- merge(
    ranked_taxonomy,
    rank_lookup,
    by = "TaxID",
    all.x = TRUE,
    sort = FALSE
  )
  
  # Final cleanup
  ranked_taxonomy <- ranked_taxonomy[
    !is.na(ranked_taxonomy$TaxID) &
      !is.na(ranked_taxonomy$Scientific.name),
    ,
    drop = FALSE
  ]
  
  
  ranked_taxonomy <- ranked_taxonomy[
    !duplicated(ranked_taxonomy$TaxID),
    ,
    drop = FALSE
  ]
  
  
  message(
    "NCBI ranked taxonomy records retained: ",
    format(
      nrow(ranked_taxonomy),
      big.mark = ","
    )
  )
  
  # Save databases
  message("")
  message("Saving local NCBI taxonomy databases...")
  
  saveRDS(
    ncbi_names,
    name_lookup_output
  )
  
  saveRDS(
    ranked_taxonomy,
    ranked_taxonomy_output
  )
  
  
  # Validate output against aRborist readers
  message("")
  message("Validating local taxonomy databases...")
  
  if (exists(
    ".read_ncbi_name_lookup",
    mode = "function"
  )) {
    
    name_test <- .read_ncbi_name_lookup(
      name_lookup_output
    )
    
    if (nrow(name_test) == 0) {
      stop(
        "The NCBI name lookup database was created ",
        "but contains zero records."
      )
    }
  }
  
  if (exists(
    ".read_ncbi_ranked_taxonomy",
    mode = "function"
  )) {
    
    taxonomy_test <- .read_ncbi_ranked_taxonomy(
      ranked_taxonomy_output
    )
    
    if (nrow(taxonomy_test) == 0) {
      stop(
        "The NCBI ranked taxonomy database was created ",
        "but contains zero records."
      )
    }
  }

  # Report creation
  message("")
  message("============================================================")
  message("LOCAL NCBI TAXONOMY DATABASE SETUP COMPLETE")
  message("============================================================")
  message("")
  message(
    "Name lookup database:\n  ",
    name_lookup_output
  )
  message(
    "  Records: ",
    format(
      nrow(ncbi_names),
      big.mark = ","
    )
  )
  message("")
  message(
    "Ranked taxonomy database:\n  ",
    ranked_taxonomy_output
  )
  message(
    "  Records: ",
    format(
      nrow(ranked_taxonomy),
      big.mark = ","
    )
  )
  message("")
  message(
    "The local NCBI taxonomy databases are ready for ",
    "aRborist host assessment."
  )
  
  invisible(
    list(
      name_lookup_file =
        name_lookup_output,
      
      ranked_taxonomy_file =
        ranked_taxonomy_output,
      
      name_records =
        nrow(ncbi_names),
      
      taxonomy_records =
        nrow(ranked_taxonomy)
    )
  )
}

# region sorting
# makes sure there is a consistent, sorted set of regions at each step
sort_regions <- function(regions_to_include) {
  sorted <- sort(unique(regions_to_include))
  if (!identical(regions_to_include, sorted)) {
    message(
      "Reordering regions_to_include to alphabetical order:\n  ",
      paste(regions_to_include, collapse = ", "),
      "  -->  ",
      paste(sorted, collapse = ", ")
    )
  }
  sorted
}

# Defaults for entrez search term
# My preferred search term defaults (can be overridden by user)
default_organism_scope <- "txid4751[Organism:exp]"  # fungi

default_search_include <- c(
  "biomol_genomic[PROP]",
  "is_nuccore[filter]",
  "(00000000100[SLEN] : 00000005000[SLEN])"
)
default_search_exclude <- c(
  "mitochondrion[filter]"
)

# ============================================================
# NCBI database selection
# ============================================================

default_ncbi_database <- "nucleotide"

normalize_ncbi_database <- function(
    ncbi_database = default_ncbi_database
) {
  
  if (
    is.null(ncbi_database) ||
    length(ncbi_database) == 0 ||
    is.na(ncbi_database[1])
  ) {
    ncbi_database <- default_ncbi_database
  }
  
  x <- tolower(
    trimws(
      as.character(ncbi_database[1])
    )
  )
  
  if (x %in% c("nucleotide", "nuccore")) {
    return("nucleotide")
  }
  
  if (x %in% c(
    "biosample",
    "bio_sample",
    "bio-sample"
  )) {
    return("biosample")
  }
  
  stop(
    "Unsupported NCBI database: ",
    ncbi_database,
    "\nCurrently supported options are:",
    "\n  nucleotide",
    "\n  biosample"
  )
}


get_entrez_database_name <- function(
    ncbi_database = default_ncbi_database
) {
  
  ncbi_database <- normalize_ncbi_database(
    ncbi_database
  )
  
  switch(
    ncbi_database,
    nucleotide = "nuccore",
    biosample = "biosample"
  )
}


# Nucleotide defaults remain unchanged
default_search_include <- c(
  "biomol_genomic[PROP]",
  "is_nuccore[filter]",
  "(00000000100[SLEN] : 00000005000[SLEN])"
)

default_search_exclude <- c(
  "mitochondrion[filter]"
)


# BioSample has completely different searchable properties,
# so do NOT inherit the Nucleotide filters by default.
default_biosample_search_include <- character(0)
default_biosample_search_exclude <- character(0)


# Initialize .GlobalEnv overrides if they don't exist
if (!exists("organism_scope", .GlobalEnv) ||
    is.null(get("organism_scope", .GlobalEnv)) ||
    !nzchar(get("organism_scope", .GlobalEnv))) {
  assign("organism_scope", default_organism_scope, .GlobalEnv)
}

if (!exists("search_include", .GlobalEnv) ||
    is.null(get("search_include", .GlobalEnv))) {
  assign("search_include", default_search_include, .GlobalEnv)
}

if (!exists("search_exclude", .GlobalEnv) ||
    is.null(get("search_exclude", .GlobalEnv))) {
  assign("search_exclude", default_search_exclude, .GlobalEnv)
}

# making full search term for entrez
compose_entrez_term <- function(
    taxon,
    organism_scope = NULL,
    include_filters = NULL,
    exclude_filters = NULL,
    ncbi_database = get0(
      "ncbi_database",
      envir = .GlobalEnv,
      ifnotfound = default_ncbi_database
    )
) {
  
  ncbi_database <- normalize_ncbi_database(
    ncbi_database
  )
  
  # ============================================================
  # 1. Resolve organism scope
  # ============================================================
  
  if (is.null(organism_scope)) {
    
    organism_scope <- get0(
      "organism_scope",
      envir = .GlobalEnv,
      ifnotfound = default_organism_scope
    )
  }
  
  
  # ============================================================
  # 2. Resolve database-specific include/exclude defaults
  # ============================================================
  
  if (is.null(include_filters)) {
    
    if (ncbi_database == "nucleotide") {
      
      include_filters <- get0(
        "search_include",
        envir = .GlobalEnv,
        ifnotfound = default_search_include
      )
      
    } else {
      
      include_filters <- get0(
        "biosample_search_include",
        envir = .GlobalEnv,
        ifnotfound = default_biosample_search_include
      )
    }
  }
  
  
  if (is.null(exclude_filters)) {
    
    if (ncbi_database == "nucleotide") {
      
      exclude_filters <- get0(
        "search_exclude",
        envir = .GlobalEnv,
        ifnotfound = default_search_exclude
      )
      
    } else {
      
      exclude_filters <- get0(
        "biosample_search_exclude",
        envir = .GlobalEnv,
        ifnotfound = default_biosample_search_exclude
      )
    }
  }
  
  
  # ============================================================
  # 3. Base organism term
  # ============================================================
  
  q <- sprintf(
    '"%s"[Organism]',
    taxon
  )
  
  
  # ============================================================
  # 4. Add organism scope
  # ============================================================
  
  if (
    !is.null(organism_scope) &&
    length(organism_scope) > 0 &&
    !is.na(organism_scope[1]) &&
    nzchar(organism_scope[1])
  ) {
    
    q <- paste(
      q,
      organism_scope,
      sep = " AND "
    )
  }
  
  
  # ============================================================
  # 5. Positive filters
  # ============================================================
  
  include_filters <- as.character(
    include_filters
  )
  
  include_filters <- include_filters[
    !is.na(include_filters) &
      nzchar(include_filters)
  ]
  
  if (length(include_filters) > 0) {
    
    include_str <- paste(
      include_filters,
      collapse = " AND "
    )
    
    q <- paste(
      q,
      include_str,
      sep = " AND "
    )
  }
  
  
  # ============================================================
  # 6. Negative filters
  # ============================================================
  
  exclude_filters <- as.character(
    exclude_filters
  )
  
  exclude_filters <- exclude_filters[
    !is.na(exclude_filters) &
      nzchar(exclude_filters)
  ]
  
  if (length(exclude_filters) > 0) {
    
    not_str <- paste(
      paste(
        "NOT",
        exclude_filters
      ),
      collapse = " "
    )
    
    q <- paste(
      q,
      not_str
    )
  }
  
  
  q
}



# default metadata categories to keep in search
if (!exists("metadata_categories_keep", .GlobalEnv)) {
  metadata_categories_keep <- c(
    "GBSeq_locus","GBSeq_length","GBSeq_strandedness","GBSeq_moltype",
    "GBSeq_update.date","GBSeq_create.date","GBSeq_definition",
    "GBSeq_accession.version","GBSeq_project","GBSeq_organism","GBSeq_taxonomy",
    "GBSeq_sequence","GBSeq_feature.table","_title","_journal","ref_id","pubmed"
  )
}

# setup oroject structure
setup_project_structure <- function(project_dir,
                                    subdirs = c("intermediate_files", "metadata_files")) {
  if (!dir.exists(project_dir)) dir.create(project_dir)
  for (dir in subdirs) {
    full_path <- file.path(project_dir, dir)
    if (!dir.exists(full_path)) dir.create(full_path, recursive = TRUE)
  }
  setwd(project_dir)
}


configure_entrez_key <- function(quiet = FALSE) {
  
  key <- ""
  
  if (exists("ncbi_api_key", envir = .GlobalEnv, inherits = FALSE)) {
    key <- get("ncbi_api_key", envir = .GlobalEnv)
  } else if (exists("api_key", envir = .GlobalEnv, inherits = FALSE)) {
    key <- get("api_key", envir = .GlobalEnv)
  }
  
  if (is.null(key) || length(key) == 0 || is.na(key[1])) {
    key <- ""
  } else {
    key <- trimws(as.character(key[1]))
  }
  
  if (!nzchar(key)) {
    key <- trimws(Sys.getenv("NCBI_API_KEY", unset = ""))
  }
  
  if (!nzchar(key)) {
    key <- trimws(Sys.getenv("ENTREZ_KEY", unset = ""))
  }
  
  if (nzchar(key)) {
    rentrez::set_entrez_key(key)
    
    if (!quiet) {
      message("NCBI API key configured for rentrez.")
    }
    
    return(invisible(TRUE))
  }
  
  if (!quiet) {
    message("No NCBI API key found; using the slower no-key request rate.")
  }
  
  invisible(FALSE)
}

# Sleep helper (uses ncbi_api_key from global env)
# this is not using the "10 requests/sec with API, 3 requests/sec without" timing because I noticed the requests were being bunched up and sent in groups, resulting in a noticable percentage of my requests getting denied. 
# you can mess with the timings if you want, but watch out for denied requests
# Find and configure an NCBI API key.
#
# Priority:
#   1. Global ncbi_api_key object
#   2. Global api_key object
#   3. NCBI_API_KEY environment variable
#   4. ENTREZ_KEY environment variable
get_sleep_duration <- function() {
  
  key_present <- configure_entrez_key(quiet = TRUE)
  
  if (key_present) {
    0.2
  } else {
    0.5
  }
}

# a probably not very accurate countdown timer when downloading metadata
format_duration <- function(seconds) {
  
  if (is.na(seconds) || !is.finite(seconds)) {
    return("unknown")
  }
  
  seconds <- max(0, round(seconds))
  
  days <- seconds %/% 86400
  hours <- (seconds %% 86400) %/% 3600
  minutes <- (seconds %% 3600) %/% 60
  secs <- seconds %% 60
  
  parts <- character(0)
  
  if (days > 0) {
    parts <- c(parts, paste0(days, "d"))
  }
  
  if (hours > 0 || days > 0) {
    parts <- c(parts, paste0(hours, "h"))
  }
  
  if (minutes > 0 || hours > 0 || days > 0) {
    parts <- c(parts, paste0(minutes, "m"))
  }
  
  if (length(parts) == 0) {
    parts <- paste0(secs, "s")
  }
  
  paste(parts, collapse = " ")
}

format_progress <- function(
    current,
    total,
    start_time,
    min_items_for_eta = 5
) {
  
  elapsed_seconds <- as.numeric(
    difftime(
      Sys.time(),
      start_time,
      units = "secs"
    )
  )
  
  percent_complete <- if (total > 0) {
    current / total * 100
  } else {
    0
  }
  
  if (
    current >= min_items_for_eta &&
    current > 0 &&
    total > current
  ) {
    
    average_seconds_per_item <- elapsed_seconds / current
    
    remaining_seconds <- average_seconds_per_item *
      (total - current)
    
    remaining_text <- format_duration(
      remaining_seconds
    )
    
  } else if (current >= total && total > 0) {
    
    remaining_text <- "0s"
    
  } else {
    
    remaining_text <- "calculating"
  }
  
  paste0(
    current,
    " / ",
    total,
    " (",
    sprintf("%.1f", percent_complete),
    "%)",
    " | elapsed: ",
    format_duration(elapsed_seconds),
    " | remaining: ",
    remaining_text
  )
}


# Save the run options for this project so it's reproducible later
# still need to implement this in a meaningful way
save_project_config <- function(project_dir = getwd(),
                                project_name,
                                taxa_of_interest,
                                regions_to_include = NULL,
                                max_acc_per_taxa = NULL,
                                min_region_requirement = NULL,
                                my_lab_sequences = "",
                                acc_to_exclude = NULL,
                                organism_scope = if (exists("organism_scope", .GlobalEnv)) get("organism_scope", .GlobalEnv) else NULL,
                                search_options = if (exists("search_options", .GlobalEnv)) get("search_options", .GlobalEnv) else NULL,
                                ncbi_api_key_present = {
                                  if (exists("ncbi_api_key", .GlobalEnv)) nzchar(get("ncbi_api_key", .GlobalEnv))
                                  else nzchar(Sys.getenv("NCBI_API_KEY"))
                                }) {
  
  cfg <- list(
    project_name = project_name,
    saved_at     = as.character(Sys.time()),
    options = list(
      taxa_of_interest       = taxa_of_interest,
      regions_to_include     = regions_to_include,
      max_acc_per_taxa       = max_acc_per_taxa,
      min_region_requirement = min_region_requirement,
      my_lab_sequences       = my_lab_sequences,
      acc_to_exclude         = acc_to_exclude,
      organism_scope         = organism_scope,
      search_options         = search_options,
      ncbi_api_key_present   = ncbi_api_key_present
    )
  )
  
  if (!dir.exists(project_dir)) dir.create(project_dir, recursive = TRUE)
  yaml::write_yaml(cfg, file.path(project_dir, "config.yml"))
  message("Saved project config: ", file.path(project_dir, "config.yml"))
  invisible(cfg)
}

# Load a project's config and return a named list (does not auto-assign)
load_project_config <- function(projects_dir, project_name) {
  cfg_path <- file.path(projects_dir, project_name, "config.yml")
  if (!file.exists(cfg_path)) stop("No config.yml found for project: ", cfg_path)
  yaml::read_yaml(cfg_path)
}

# ============================================================



# ============================================================
# arborist main functions
# ============================================================

# cleaning taxon names, prevents weirdness from names with square brackets
clean_taxon_name <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  
  # Remove square brackets around taxon names, e.g. [Neocosmospora] -> Neocosmospora
  x <- gsub("^\\[(.*)\\]$", "\\1", x)
  
  # Clean extra whitespace again
  x <- trimws(x)
  
  x[x == ""] <- NA_character_
  
  x
}


.biosample_esummary_to_accessions <- function(
    summaries,
    search_group
) {
  
  if (
    is.null(summaries) ||
    length(summaries) == 0
  ) {
    
    return(
      data.frame(
        Accession = character(0),
        EntrezUID = character(0),
        search_group = character(0),
        ncbi_database = character(0),
        stringsAsFactors = FALSE
      )
    )
  }
  
  
  get_summary_value <- function(
    x,
    field
  ) {
    
    nms <- names(x)
    
    if (
      is.null(nms) ||
      length(nms) == 0
    ) {
      return(NA_character_)
    }
    
    hit <- which(
      tolower(nms) ==
        tolower(field)
    )
    
    if (length(hit) == 0) {
      return(NA_character_)
    }
    
    value <- x[[hit[1]]]
    
    if (
      is.null(value) ||
      length(value) == 0
    ) {
      return(NA_character_)
    }
    
    as.character(value[1])
  }
  
  
  summary_names <- names(
    summaries
  )
  
  
  rows <- lapply(
    seq_along(summaries),
    function(i) {
      
      this_summary <- summaries[[i]]
      
      accession <- get_summary_value(
        this_summary,
        "accession"
      )
      
      uid <- get_summary_value(
        this_summary,
        "uid"
      )
      
      if (
        (is.na(uid) || !nzchar(uid)) &&
        !is.null(summary_names) &&
        length(summary_names) >= i
      ) {
        uid <- summary_names[i]
      }
      
      data.frame(
        Accession = accession,
        EntrezUID = uid,
        search_group = search_group,
        ncbi_database = "biosample",
        stringsAsFactors = FALSE
      )
    }
  )
  
  
  out <- dplyr::bind_rows(
    rows
  )
  
  
  out <- out[
    !is.na(out$Accession) &
      nzchar(out$Accession),
    ,
    drop = FALSE
  ]
  
  
  out <- dplyr::distinct(
    out,
    Accession,
    .keep_all = TRUE
  )
  
  
  out
}

.write_accession_checkpoint <- function(
    x,
    checkpoint_file
) {
  
  checkpoint_dir <- dirname(
    checkpoint_file
  )
  
  
  if (!dir.exists(checkpoint_dir)) {
    dir.create(
      checkpoint_dir,
      recursive = TRUE
    )
  }
  
  
  temp_file <- tempfile(
    pattern = "accession_checkpoint_",
    tmpdir = checkpoint_dir,
    fileext = ".csv"
  )
  
  
  write.csv(
    x,
    temp_file,
    row.names = FALSE,
    quote = FALSE
  )
  
  
  renamed <- file.rename(
    temp_file,
    checkpoint_file
  )
  
  
  if (!renamed) {
    
    copied <- file.copy(
      temp_file,
      checkpoint_file,
      overwrite = TRUE
    )
    
    
    unlink(
      temp_file
    )
    
    
    if (!copied) {
      
      stop(
        "Could not safely write accession checkpoint: ",
        checkpoint_file
      )
    }
  }
  
  
  invisible(
    checkpoint_file
  )
}


# Fetch accessions from NCBI - searching using entrez
fetch_accessions_for_taxon <- function(
    taxon,
    max_acc = max_acc_per_taxa,
    organism_scope = NULL,
    include_filters = NULL,
    exclude_filters = NULL,
    checkpoint_file = NULL,
    checkpoint_every = 500,
    resume = TRUE,
    overwrite = FALSE,
    accession_fetch_batch_size = 500,
    page_max_retries = 3,
    retry_wait = 5,
    max_history_refreshes = 3,
    ncbi_database = get0(
      "ncbi_database",
      envir = .GlobalEnv,
      ifnotfound = default_ncbi_database
    )
) {
  
  ncbi_database <- normalize_ncbi_database(
    ncbi_database
  )
  
  entrez_db <- get_entrez_database_name(
    ncbi_database
  )
  
  cat(
    "Searching term:",
    taxon,
    "\n"
  )
  
  cat(
    "NCBI database:",
    ncbi_database,
    " (",
    entrez_db,
    ")\n",
    sep = ""
  )
  
  
  # ============================================================
  # Existing checkpoint
  # ============================================================
  
  existing_accessions <- data.frame(
    Accession = character(0),
    EntrezUID = character(0),
    search_group = character(0),
    ncbi_database = character(0),
    stringsAsFactors = FALSE
  )
  
  
  if (
    !is.null(checkpoint_file) &&
    file.exists(checkpoint_file)
  ) {
    
    if (overwrite) {
      
      message(
        "Removing existing accession checkpoint: ",
        checkpoint_file
      )
      
      file.remove(
        checkpoint_file
      )
      
    } else if (resume) {
      
      existing_accessions <- read.csv(
        checkpoint_file,
        stringsAsFactors = FALSE,
        colClasses = "character"
      )
      
      
      if (
        !"Accession" %in%
        names(existing_accessions)
      ) {
        
        stop(
          "Existing accession checkpoint is missing the Accession column: ",
          checkpoint_file
        )
      }
      
      
      if (
        !"search_group" %in%
        names(existing_accessions)
      ) {
        
        existing_accessions$search_group <- rep(
          taxon,
          nrow(existing_accessions)
        )
      }
      
      
      if (
        !"EntrezUID" %in%
        names(existing_accessions)
      ) {
        
        existing_accessions$EntrezUID <- NA_character_
      }
      
      
      if (
        !"ncbi_database" %in%
        names(existing_accessions)
      ) {
        
        existing_accessions$ncbi_database <-
          ncbi_database
      }
      
      
      stored_databases <- unique(
        existing_accessions$ncbi_database[
          !is.na(existing_accessions$ncbi_database) &
            nzchar(existing_accessions$ncbi_database)
        ]
      )
      
      
      if (
        length(stored_databases) > 0 &&
        any(stored_databases != ncbi_database)
      ) {
        
        stop(
          "Accession checkpoint belongs to a different NCBI database:\n  ",
          checkpoint_file,
          "\nStored database: ",
          paste(
            stored_databases,
            collapse = ", "
          ),
          "\nRequested database: ",
          ncbi_database
        )
      }
      
      
      existing_accessions <-
        dplyr::distinct(
          existing_accessions,
          Accession,
          .keep_all = TRUE
        )
      
      
      message(
        "Resuming from checkpoint with ",
        nrow(existing_accessions),
        " accession(s): ",
        checkpoint_file
      )
    }
  }
  
  
  # ============================================================
  # Build Entrez search
  # ============================================================
  
  if (
    exists(
      "raw_entrez_terms",
      envir = .GlobalEnv
    ) &&
    taxon %in%
    names(
      get(
        "raw_entrez_terms",
        envir = .GlobalEnv
      )
    )
  ) {
    
    filters <- get(
      "raw_entrez_terms",
      envir = .GlobalEnv
    )[[taxon]]
    
    message(
      "Using raw Entrez query for ",
      taxon,
      ": ",
      filters
    )
    
  } else {
    
    filters <- compose_entrez_term(
      taxon = taxon,
      organism_scope = organism_scope,
      include_filters = include_filters,
      exclude_filters = exclude_filters,
      ncbi_database = ncbi_database
    )
  }
  
  
  message(
    "\nFINAL ENTREZ QUERY:\n",
    filters,
    "\n"
  )
  
  
  # ============================================================
  # Helper: create a fresh Web History
  # ============================================================
  
  create_history_search <- function() {
    
    last_error <- NULL
    
    
    for (attempt in seq_len(page_max_retries)) {
      
      search_result <- tryCatch(
        {
          
          rentrez::entrez_search(
            db = entrez_db,
            term = filters,
            use_history = TRUE,
            retmax = 0
          )
          
        },
        error = function(e) {
          
          last_error <<-
            conditionMessage(e)
          
          NULL
        }
      )
      
      
      if (
        !is.null(search_result) &&
        !is.null(search_result$web_history)
      ) {
        
        return(
          search_result
        )
      }
      
      
      if (
        attempt < page_max_retries
      ) {
        
        wait_seconds <-
          retry_wait * attempt
        
        
        message(
          "Entrez search / Web History creation failed ",
          "(attempt ",
          attempt,
          " / ",
          page_max_retries,
          "). Waiting ",
          wait_seconds,
          " seconds before retry..."
        )
        
        
        Sys.sleep(
          wait_seconds
        )
      }
    }
    
    
    stop(
      "Unable to create NCBI Web History after ",
      page_max_retries,
      " attempt(s).\n",
      "Last error: ",
      last_error
    )
  }
  
  
  # ============================================================
  # Initial Web History
  # ============================================================
  
  search <- create_history_search()
  
  
  total_accession_count <- as.integer(
    search$count
  )
  
  
  if (
    is.na(total_accession_count) ||
    total_accession_count == 0
  ) {
    
    cat(
      "No accessions found for:",
      taxon,
      "\n"
    )
    
    return(
      data.frame(
        Accession = character(0),
        EntrezUID = character(0),
        search_group = character(0),
        ncbi_database = character(0),
        stringsAsFactors = FALSE
      )
    )
  }
  
  
  max_n <- if (
    identical(
      max_acc,
      "max"
    )
  ) {
    Inf
  } else {
    as.numeric(max_acc)
  }
  
  
  pull_n <- min(
    total_accession_count,
    max_n
  )
  
  
  cat(
    total_accession_count,
    " accessions available for ",
    taxon,
    " - pulling a maximum of ",
    pull_n,
    "\n",
    sep = ""
  )
  
  
  already_retrieved <- nrow(
    existing_accessions
  )
  
  
  if (
    already_retrieved >=
    pull_n
  ) {
    
    message(
      "Checkpoint already contains ",
      already_retrieved,
      " accession(s), which meets or exceeds the requested total of ",
      pull_n,
      "."
    )
    
    
    df <- existing_accessions[
      seq_len(
        min(
          nrow(existing_accessions),
          pull_n
        )
      ),
      ,
      drop = FALSE
    ]
    
    
    return(df)
  }
  
  
  # ============================================================
  # Helper: retrieve one page with retries and Web History refresh
  # ============================================================
  
  fetch_accession_page <- function(
    seq_start,
    this_retmax
  ) {
    
    history_refresh_count <- 0L
    last_error <- NULL
    
    
    repeat {
      
      # --------------------------------------------------------
      # Try current Web History several times
      # --------------------------------------------------------
      
      for (attempt in seq_len(page_max_retries)) {
        
        new_df <- tryCatch(
          {
            
            # --------------------------------------------------
            # Nucleotide
            # --------------------------------------------------
            
            if (
              ncbi_database ==
              "nucleotide"
            ) {
              
              recs <- rentrez::entrez_fetch(
                db = "nuccore",
                web_history = search$web_history,
                rettype = "acc",
                retmax = this_retmax,
                retstart = seq_start
              )
              
              
              new_acc <- unlist(
                strsplit(
                  recs,
                  "\\s+"
                )
              )
              
              
              new_acc <- new_acc[
                nzchar(new_acc)
              ]
              
              
              if (
                length(new_acc) !=
                this_retmax
              ) {
                
                stop(
                  "NCBI returned ",
                  length(new_acc),
                  " accession(s), but ",
                  this_retmax,
                  " were requested at retstart ",
                  seq_start,
                  "."
                )
              }
              
              
              data.frame(
                Accession = new_acc,
                EntrezUID = rep(
                  NA_character_,
                  length(new_acc)
                ),
                search_group = rep(
                  taxon,
                  length(new_acc)
                ),
                ncbi_database = rep(
                  "nucleotide",
                  length(new_acc)
                ),
                stringsAsFactors = FALSE
              )
              
              
            } else {
              
              # ------------------------------------------------
              # BioSample
              # ------------------------------------------------
              
              summaries <-
                rentrez::entrez_summary(
                  db = "biosample",
                  web_history =
                    search$web_history,
                  retmax =
                    this_retmax,
                  retstart =
                    seq_start,
                  always_return_list =
                    TRUE
                )
              
              
              new_df <-
                .biosample_esummary_to_accessions(
                  summaries =
                    summaries,
                  search_group =
                    taxon
                )
              
              
              if (
                nrow(new_df) !=
                this_retmax
              ) {
                
                stop(
                  "NCBI returned ",
                  nrow(new_df),
                  " BioSample accession(s), but ",
                  this_retmax,
                  " were requested at retstart ",
                  seq_start,
                  "."
                )
              }
              
              
              new_df
            }
            
          },
          error = function(e) {
            
            last_error <<-
              conditionMessage(e)
            
            NULL
          }
        )
        
        
        if (!is.null(new_df)) {
          return(
            new_df
          )
        }
        
        
        if (
          attempt <
          page_max_retries
        ) {
          
          wait_seconds <-
            retry_wait * attempt
          
          
          message(
            "  Accession page failed at retstart ",
            seq_start,
            " (attempt ",
            attempt,
            " / ",
            page_max_retries,
            ")."
          )
          
          
          message(
            "  Error: ",
            last_error
          )
          
          
          message(
            "  Waiting ",
            wait_seconds,
            " seconds before retry..."
          )
          
          
          Sys.sleep(
            wait_seconds
          )
        }
      }
      
      
      # ========================================================
      # Current Web History repeatedly failed
      # ========================================================
      
      if (
        history_refresh_count >=
        max_history_refreshes
      ) {
        
        stop(
          "Accession retrieval failed at retstart ",
          seq_start,
          " after ",
          page_max_retries,
          " retry attempt(s) per Web History and ",
          history_refresh_count,
          " Web History refresh(es).\n",
          "The existing accession checkpoint will be preserved.\n",
          "Last NCBI error: ",
          last_error
        )
      }
      
      
      history_refresh_count <-
        history_refresh_count + 1L
      
      
      message("")
      
      message(
        "  Refreshing NCBI Web History after repeated failure..."
      )
      
      
      message(
        "  Web History refresh ",
        history_refresh_count,
        " / ",
        max_history_refreshes
      )
      
      
      refreshed_search <- create_history_search()
      
      
      refreshed_count <- as.integer(
        refreshed_search$count
      )
      
      
      # --------------------------------------------------------
      # Safety check
      #
      # If the search count changed, do not continue blindly
      # using the same retstart positions.
      # --------------------------------------------------------
      
      if (
        is.na(refreshed_count) ||
        refreshed_count !=
        total_accession_count
      ) {
        
        stop(
          "NCBI search result count changed while refreshing Web History.\n",
          "Original count: ",
          total_accession_count,
          "\nRefreshed count: ",
          refreshed_count,
          "\nStopping rather than risk skipping or duplicating accessions.\n",
          "The existing accession checkpoint will be preserved."
        )
      }
      
      
      search <<-
        refreshed_search
      
      
      message(
        "  Fresh Web History created successfully."
      )
      
      
      message(
        "  Retrying the same accession page at retstart ",
        seq_start,
        "."
      )
      
      
      Sys.sleep(
        retry_wait
      )
    }
  }
  
  
  # ============================================================
  # Retrieve accession pages
  #
  # Keep current page size = 50 for now.
  # ============================================================
  
  record_chunks <- list()
  
  new_since_checkpoint <- 0L
  
  
  for (
    seq_start in seq(
      already_retrieved,
      pull_n - 1,
      by = accession_fetch_batch_size
    )
  ) {
    
    this_retmax <- min(
      accession_fetch_batch_size,
      pull_n - seq_start
    )
    
    
    new_df <- fetch_accession_page(
      seq_start = seq_start,
      this_retmax = this_retmax
    )
    
    
    if (
      !is.null(new_df) &&
      nrow(new_df) > 0
    ) {
      
      record_chunks[[
        length(record_chunks) + 1L
      ]] <- new_df
      
      
      new_since_checkpoint <-
        new_since_checkpoint +
        nrow(new_df)
    }
    
    
    # ==========================================================
    # Accession checkpoint
    # ==========================================================
    
    if (
      !is.null(checkpoint_file) &&
      new_since_checkpoint >=
      checkpoint_every
    ) {
      
      newly_retrieved <-
        dplyr::bind_rows(
          record_chunks
        )
      
      
      checkpoint_df <-
        dplyr::bind_rows(
          existing_accessions,
          newly_retrieved
        )
      
      
      checkpoint_df <-
        dplyr::distinct(
          checkpoint_df,
          Accession,
          .keep_all = TRUE
        )
      
      
      .write_accession_checkpoint(
        checkpoint_df,
        checkpoint_file
      )
      
      
      message(
        "Accession checkpoint written: ",
        checkpoint_file,
        " | ",
        nrow(checkpoint_df),
        " unique accession(s)"
      )
      
      
      existing_accessions <-
        checkpoint_df
      
      
      record_chunks <- list()
      new_since_checkpoint <- 0L
    }
    
    
    Sys.sleep(
      get_sleep_duration()
    )
  }
  
  
  # ============================================================
  # Combine final accession set
  # ============================================================
  
  remaining_df <- if (
    length(record_chunks) > 0
  ) {
    
    dplyr::bind_rows(
      record_chunks
    )
    
  } else {
    
    data.frame(
      Accession = character(0),
      EntrezUID = character(0),
      search_group = character(0),
      ncbi_database = character(0),
      stringsAsFactors = FALSE
    )
  }
  
  
  df <- dplyr::bind_rows(
    existing_accessions,
    remaining_df
  )
  
  
  df <- dplyr::distinct(
    df,
    Accession,
    .keep_all = TRUE
  )
  
  
  # ============================================================
  # Validate final accession count
  # ============================================================
  
  if (
    nrow(df) !=
    pull_n
  ) {
    
    stop(
      "Accession retrieval ended with ",
      nrow(df),
      " unique accession(s), but ",
      pull_n,
      " were expected.\n",
      "The checkpoint will be preserved and this search group ",
      "will NOT be treated as successfully completed."
    )
  }
  
  
  # ============================================================
  # Final accession checkpoint
  # ============================================================
  
  if (
    !is.null(checkpoint_file)
  ) {
    
    .write_accession_checkpoint(
      df,
      checkpoint_file
    )
    
    
    message(
      "Final accession checkpoint written: ",
      checkpoint_file,
      " | ",
      nrow(df),
      " unique accession(s)"
    )
  }
  
  
  cat(
    "Accession retrieval for ",
    taxon,
    " successful: ",
    nrow(df),
    " accessions\n\n",
    sep = ""
  )
  
  
  df
}


# pulling accessions using provided taxa names
get_accessions_for_all_taxa <- function(
    taxa_list,
    max_acc_per_taxa,
    organism_scope = NULL,
    include_filters = NULL,
    exclude_filters = NULL,
    timing_file = "./intermediate_files/fetch_times_accessions.csv",
    checkpoint_every = 500,
    resume = TRUE,
    overwrite = FALSE,
    accession_fetch_batch_size = 500,
    page_max_retries = 3,
    retry_wait = 5,
    max_history_refreshes = 3,
    ncbi_database = get0(
      "ncbi_database",
      envir = .GlobalEnv,
      ifnotfound = default_ncbi_database
    )
) {
  
  ncbi_database <- normalize_ncbi_database(
    ncbi_database
  )
  
  
  taxa_frame_acc <- vector(
    "list",
    length(taxa_list)
  )
  
  
  timing_log <- data.frame(
    Taxon = character(),
    NCBI_database = character(),
    Num_accessions = integer(),
    Start_time = character(),
    End_time = character(),
    Elapsed_minutes = numeric(),
    Status = character(),
    Error = character(),
    stringsAsFactors = FALSE
  )
  
  
  overall_start <- Sys.time()
  
  
  for (
    i in seq_along(taxa_list)
  ) {
    
    term <- taxa_list[i]
    
    
    cat(
      "\n=== Starting ",
      term,
      " (",
      i,
      " of ",
      length(taxa_list),
      ") ===\n",
      sep = ""
    )
    
    
    safe_term <- gsub(
      "[^A-Za-z0-9_.-]+",
      "_",
      term
    )
    
    
    checkpoint_prefix <- if (
      ncbi_database ==
      "nucleotide"
    ) {
      "Accessions_for_"
    } else {
      "BioSample_accessions_for_"
    }
    
    
    outfile_name <- file.path(
      "./intermediate_files",
      paste0(
        checkpoint_prefix,
        safe_term,
        ".csv"
      )
    )
    
    
    start_time <- Sys.time()
    
    
    retrieval_error <- NULL
    
    
    tempdf <- tryCatch(
      {
        
        fetch_accessions_for_taxon(
          taxon = term,
          max_acc =
            max_acc_per_taxa,
          organism_scope =
            organism_scope,
          include_filters =
            include_filters,
          exclude_filters =
            exclude_filters,
          checkpoint_file =
            outfile_name,
          checkpoint_every =
            checkpoint_every,
          resume =
            resume,
          overwrite =
            overwrite,
          accession_fetch_batch_size =
            accession_fetch_batch_size,
          page_max_retries =
            page_max_retries,
          retry_wait =
            retry_wait,
          max_history_refreshes =
            max_history_refreshes,
          ncbi_database =
            ncbi_database
        )
      },
      error = function(e) {
        
        retrieval_error <<-
          conditionMessage(e)
        
        NULL
      }
    )
    
    
    end_time <- Sys.time()
    
    
    elapsed <- as.numeric(
      difftime(
        end_time,
        start_time,
        units = "mins"
      )
    )
    
    
    # ==========================================================
    # Retrieval failed
    # ==========================================================
    
    if (is.null(tempdf)) {
      
      cat(
        "\nERROR while retrieving ",
        term,
        ":\n",
        retrieval_error,
        "\n",
        sep = ""
      )
      
      
      checkpoint_rows <- NA_integer_
      
      
      if (file.exists(outfile_name)) {
        
        checkpoint_rows <- tryCatch(
          {
            
            checkpoint_check <- read.csv(
              outfile_name,
              stringsAsFactors = FALSE,
              colClasses = "character"
            )
            
            
            nrow(
              checkpoint_check
            )
          },
          error = function(e) {
            NA_integer_
          }
        )
        
        
        message(
          "Existing accession checkpoint preserved: ",
          outfile_name,
          if (!is.na(checkpoint_rows)) {
            paste0(
              " | ",
              checkpoint_rows,
              " accession(s)"
            )
          } else {
            ""
          }
        )
      }
      
      
      timing_log <- rbind(
        timing_log,
        data.frame(
          Taxon = term,
          NCBI_database =
            ncbi_database,
          Num_accessions =
            checkpoint_rows,
          Start_time = format(
            start_time,
            "%Y-%m-%d %H:%M:%S"
          ),
          End_time = format(
            end_time,
            "%Y-%m-%d %H:%M:%S"
          ),
          Elapsed_minutes =
            round(
              elapsed,
              2
            ),
          Status =
            "FAILED",
          Error =
            retrieval_error,
          stringsAsFactors = FALSE
        )
      )
      
      
      write.csv(
        timing_log,
        timing_file,
        row.names = FALSE
      )
      
      
      stop(
        "Accession retrieval stopped because search group '",
        term,
        "' did not complete.\n",
        "The existing accession checkpoint has been preserved.\n",
        "Rerun with resume = TRUE after the NCBI error has cleared.\n",
        "Original error:\n",
        retrieval_error,
        call. = FALSE
      )
    }
    
    
    # ==========================================================
    # Retrieval succeeded
    # ==========================================================
    
    cat(
      sprintf(
        "Finished %s in %.2f minutes\n",
        term,
        elapsed
      )
    )
    
    
    timing_log <- rbind(
      timing_log,
      data.frame(
        Taxon = term,
        NCBI_database =
          ncbi_database,
        Num_accessions =
          nrow(tempdf),
        Start_time = format(
          start_time,
          "%Y-%m-%d %H:%M:%S"
        ),
        End_time = format(
          end_time,
          "%Y-%m-%d %H:%M:%S"
        ),
        Elapsed_minutes =
          round(
            elapsed,
            2
          ),
        Status =
          "COMPLETE",
        Error =
          "",
        stringsAsFactors = FALSE
      )
    )
    
    
    taxa_frame_acc[[i]] <-
      tempdf
    
    
    write.csv(
      timing_log,
      timing_file,
      row.names = FALSE
    )
  }
  
  
  # ============================================================
  # All search groups completed successfully
  # ============================================================
  
  total_elapsed <- as.numeric(
    difftime(
      Sys.time(),
      overall_start,
      units = "mins"
    )
  )
  
  
  cat(
    "\nAll taxa completed in ",
    round(total_elapsed, 2),
    " minutes.\n",
    sep = ""
  )
  
  
  non_empty <- Filter(
    function(x) {
      !is.null(x) &&
        nrow(x) > 0
    },
    taxa_frame_acc
  )
  
  
  all_memberships <- if (
    length(non_empty)
  ) {
    
    dplyr::bind_rows(
      non_empty
    )
    
  } else {
    
    data.frame(
      Accession = character(0),
      EntrezUID = character(0),
      search_group = character(0),
      ncbi_database = character(0),
      stringsAsFactors = FALSE
    )
  }
  
  
  write.csv(
    all_memberships,
    "./intermediate_files/all_pulled_accession_memberships.csv",
    row.names = FALSE
  )
  
  
  accession_list <-
    dplyr::distinct(
      all_memberships,
      Accession,
      .keep_all = TRUE
    )
  
  
  write.csv(
    accession_list,
    "./intermediate_files/all_pulled_accessions.csv",
    row.names = FALSE
  )
  
  
  write.csv(
    timing_log,
    timing_file,
    row.names = FALSE
  )
  
  
  cat(
    "Timing log written to ",
    timing_file,
    "\n",
    sep = ""
  )
  
  
  accession_list
}


fetch_ncbi_metadata_batch_xml <- function(
    accessions,
    post_chunk_size = 500,
    biosample_fetch_chunk_size = 50,
    ncbi_database = get0(
      "ncbi_database",
      envir = .GlobalEnv,
      ifnotfound = default_ncbi_database
    )
) {
  
  ncbi_database <- normalize_ncbi_database(
    ncbi_database
  )
  
  
  accessions <- unique(
    trimws(
      as.character(
        accessions
      )
    )
  )
  
  
  accessions <- accessions[
    !is.na(accessions) &
      nzchar(accessions)
  ]
  
  
  if (
    length(accessions) == 0
  ) {
    stop(
      "No valid accessions were supplied."
    )
  }
  
  
  # ============================================================
  # NUCLEOTIDE
  #
  # Preserve existing EPost + History Server behavior.
  # ============================================================
  
  if (
    ncbi_database ==
    "nucleotide"
  ) {
    
    post_chunks <- split(
      accessions,
      ceiling(
        seq_along(accessions) /
          post_chunk_size
      )
    )
    
    
    message(
      "Posting ",
      length(accessions),
      " accession(s) to NCBI History Server in ",
      length(post_chunks),
      " chunk(s)..."
    )
    
    
    web_history <- NULL
    
    
    for (
      i in seq_along(
        post_chunks
      )
    ) {
      
      this_chunk <- unname(
        post_chunks[[i]]
      )
      
      
      message(
        "  Posting chunk ",
        i,
        " / ",
        length(post_chunks),
        " (",
        length(this_chunk),
        " accession(s))"
      )
      
      
      if (is.null(web_history)) {
        
        web_history <-
          rentrez::entrez_post(
            db = "nuccore",
            id = this_chunk
          )
        
      } else {
        
        web_history <-
          rentrez::entrez_post(
            db = "nuccore",
            id = this_chunk,
            web_history =
              web_history
          )
      }
      
      
      Sys.sleep(
        get_sleep_duration()
      )
    }
    
    
    message(
      "Fetching ",
      length(accessions),
      " Nucleotide record(s) from NCBI History Server..."
    )
    
    
    xml_text <-
      rentrez::entrez_fetch(
        db = "nuccore",
        web_history =
          web_history,
        rettype = "gb",
        retmode = "xml",
        retmax =
          length(accessions)
      )
    
    
    if (
      is.null(xml_text) ||
      !nzchar(xml_text)
    ) {
      
      stop(
        "NCBI returned an empty Nucleotide response for a batch of ",
        length(accessions),
        " accession(s)."
      )
    }
    
    
    doc <- XML::xmlParse(
      xml_text
    )
    
    
    record_nodes <- XML::getNodeSet(
      doc,
      "//GBSeq"
    )
    
    
    if (
      length(record_nodes) == 0
    ) {
      
      stop(
        "No GBSeq records were found in the NCBI XML response."
      )
    }
    
    
    message(
      "Requested ",
      length(accessions),
      " accession(s); NCBI returned ",
      length(record_nodes),
      " GBSeq record(s)."
    )
    
    
    return(
      record_nodes
    )
  }
  
  
  # ============================================================
  # BIOSAMPLE
  #
  # BioSample full records are returned directly as XML.
  # Keep requests small enough that accession lists remain
  # manageable.
  # ============================================================
  
  fetch_chunks <- split(
    accessions,
    ceiling(
      seq_along(accessions) /
        biosample_fetch_chunk_size
    )
  )
  
  
  message(
    "Fetching ",
    length(accessions),
    " BioSample record(s) in ",
    length(fetch_chunks),
    " XML chunk(s)..."
  )
  
  
  record_nodes <- list()
  
  
  for (
    i in seq_along(
      fetch_chunks
    )
  ) {
    
    this_chunk <- unname(
      fetch_chunks[[i]]
    )
    
    
    message(
      "  Fetching BioSample XML chunk ",
      i,
      " / ",
      length(fetch_chunks),
      " (",
      length(this_chunk),
      " accession(s))"
    )
    
    
    xml_text <-
      rentrez::entrez_fetch(
        db = "biosample",
        id = this_chunk,
        rettype = "full",
        retmode = "xml"
      )
    
    
    if (
      is.null(xml_text) ||
      !nzchar(xml_text)
    ) {
      
      stop(
        "NCBI returned an empty BioSample response for ",
        length(this_chunk),
        " accession(s)."
      )
    }
    
    
    doc <- XML::xmlParse(
      xml_text
    )
    
    
    this_nodes <- XML::getNodeSet(
      doc,
      "//BioSample"
    )
    
    
    if (
      length(this_nodes) == 0
    ) {
      
      stop(
        "No BioSample records were found in the NCBI XML response."
      )
    }
    
    
    record_nodes <- c(
      record_nodes,
      this_nodes
    )
    
    
    Sys.sleep(
      get_sleep_duration()
    )
  }
  
  
  message(
    "Requested ",
    length(accessions),
    " BioSample accession(s); NCBI returned ",
    length(record_nodes),
    " BioSample record(s)."
  )
  
  
  record_nodes
}


fetch_metadata_for_accession_batch <- function(
    accessions,
    ncbi_database = get0(
      "ncbi_database",
      envir = .GlobalEnv,
      ifnotfound = default_ncbi_database
    )
) {
  
  ncbi_database <- normalize_ncbi_database(
    ncbi_database
  )
  
  
  record_nodes <-
    fetch_ncbi_metadata_batch_xml(
      accessions = accessions,
      ncbi_database =
        ncbi_database
    )
  
  
  parser <- if (
    ncbi_database ==
    "nucleotide"
  ) {
    parse_gbseq_node
  } else {
    parse_biosample_node
  }
  
  
  metadata_list <- lapply(
    record_nodes,
    parser
  )
  
  
  metadata_df <- dplyr::bind_rows(
    metadata_list
  )
  
  
  metadata_df
}


fetch_metadata_batch_resilient <- function(
    accessions,
    min_batch_size = 1,
    max_retries = 2,
    retry_wait = 5,
    ncbi_database = get0(
      "ncbi_database",
      envir = .GlobalEnv,
      ifnotfound = default_ncbi_database
    )
) {
  
  ncbi_database <- normalize_ncbi_database(
    ncbi_database
  )
  
  
  accessions <- unique(
    trimws(
      as.character(
        accessions
      )
    )
  )
  
  
  accessions <- accessions[
    !is.na(accessions) &
      nzchar(accessions)
  ]
  
  
  if (
    length(accessions) == 0
  ) {
    return(NULL)
  }
  
  
  # ============================================================
  # First try the entire requested batch
  # ============================================================
  
  last_error <- NULL
  
  
  for (
    attempt in seq_len(
      max_retries
    )
  ) {
    
    result <- tryCatch(
      {
        
        fetch_metadata_for_accession_batch(
          accessions,
          ncbi_database =
            ncbi_database
        )
      },
      error = function(e) {
        
        last_error <<-
          conditionMessage(e)
        
        NULL
      }
    )
    
    
    if (!is.null(result)) {
      return(result)
    }
    
    
    if (
      attempt <
      max_retries
    ) {
      
      message(
        "  Batch of ",
        length(accessions),
        " failed (attempt ",
        attempt,
        " / ",
        max_retries,
        "). Waiting ",
        retry_wait,
        " seconds before retry..."
      )
      
      
      Sys.sleep(
        retry_wait
      )
    }
  }
  
  
  # ============================================================
  # Split failed batches progressively
  # ============================================================
  
  if (
    length(accessions) >
    min_batch_size
  ) {
    
    split_point <- ceiling(
      length(accessions) / 2
    )
    
    
    first_half <- accessions[
      seq_len(split_point)
    ]
    
    
    second_half <- accessions[
      seq.int(
        split_point + 1L,
        length(accessions)
      )
    ]
    
    
    message(
      "  Batch of ",
      length(accessions),
      " failed after ",
      max_retries,
      " attempt(s). Splitting into ",
      length(first_half),
      " + ",
      length(second_half),
      "."
    )
    
    
    first_result <-
      fetch_metadata_batch_resilient(
        accessions =
          first_half,
        min_batch_size =
          min_batch_size,
        max_retries =
          max_retries,
        retry_wait =
          retry_wait,
        ncbi_database =
          ncbi_database
      )
    
    
    second_result <- if (
      length(second_half) > 0
    ) {
      
      fetch_metadata_batch_resilient(
        accessions =
          second_half,
        min_batch_size =
          min_batch_size,
        max_retries =
          max_retries,
        retry_wait =
          retry_wait,
        ncbi_database =
          ncbi_database
      )
      
    } else {
      NULL
    }
    
    
    successful_results <- Filter(
      Negate(is.null),
      list(
        first_result,
        second_result
      )
    )
    
    
    if (
      length(successful_results) == 0
    ) {
      return(NULL)
    }
    
    
    return(
      dplyr::bind_rows(
        successful_results
      )
    )
  }
  
  
  # ============================================================
  # Single accession still failed
  # ============================================================
  
  message(
    "  FAILED | ",
    accessions[1],
    " | ",
    last_error
  )
  
  
  NULL
}

accession_completion_file <- function() {
  "./intermediate_files/accession_retrieval_complete.rds"
}



write_accession_completion_marker <- function(
    taxa_list,
    accession_file = "./intermediate_files/all_pulled_accessions.csv",
    ncbi_database = get0(
      "ncbi_database",
      envir = .GlobalEnv,
      ifnotfound = default_ncbi_database
    )
) {
  
  ncbi_database <- normalize_ncbi_database(
    ncbi_database
  )
  
  
  # ============================================================
  # Accession manifest must exist
  # ============================================================
  
  if (!file.exists(accession_file)) {
    
    stop(
      "Cannot write accession completion marker because the ",
      "accession manifest does not exist:\n  ",
      accession_file
    )
  }
  
  
  # ============================================================
  # Read and validate accession manifest
  # ============================================================
  
  accession_manifest <- read.csv(
    accession_file,
    stringsAsFactors = FALSE,
    colClasses = "character"
  )
  
  
  if (
    !"Accession" %in%
    names(accession_manifest)
  ) {
    
    stop(
      "Cannot write accession completion marker because the ",
      "accession manifest is missing the Accession column:\n  ",
      accession_file
    )
  }
  
  
  valid_accessions <- accession_manifest$Accession[
    !is.na(accession_manifest$Accession) &
      nzchar(trimws(accession_manifest$Accession))
  ]
  
  
  n_accessions <- length(
    unique(valid_accessions)
  )
  
  
  # ============================================================
  # An empty accession manifest should never be marked complete
  # ============================================================
  
  if (n_accessions == 0) {
    
    stop(
      "Cannot write accession completion marker because the ",
      "accession manifest contains zero accessions:\n  ",
      accession_file
    )
  }
  
  
  # ============================================================
  # Build completion marker
  # ============================================================
  
  marker <- list(
    completed = TRUE,
    taxa_list =
      as.character(taxa_list),
    n_search_groups =
      length(taxa_list),
    n_accessions =
      n_accessions,
    ncbi_database =
      ncbi_database,
    completed_at =
      as.character(Sys.time()),
    accession_file =
      accession_file
  )
  
  
  # ============================================================
  # Save marker
  # ============================================================
  
  saveRDS(
    marker,
    accession_completion_file()
  )
  
  
  message(
    "Accession retrieval completion marker written: ",
    accession_completion_file(),
    " | ",
    n_accessions,
    " unique accession(s)"
  )
  
  
  invisible(marker)
}


accession_retrieval_is_complete <- function(
    taxa_list,
    ncbi_database = get0(
      "ncbi_database",
      envir = .GlobalEnv,
      ifnotfound = default_ncbi_database
    )
) {
  
  ncbi_database <- normalize_ncbi_database(
    ncbi_database
  )
  
  
  marker_path <- accession_completion_file()
  
  
  # ============================================================
  # No marker = definitely not complete
  # ============================================================
  
  if (!file.exists(marker_path)) {
    return(FALSE)
  }
  
  
  # ============================================================
  # Read completion marker safely
  # ============================================================
  
  marker <- tryCatch(
    readRDS(marker_path),
    error = function(e) NULL
  )
  
  
  if (is.null(marker)) {
    
    message(
      "Accession completion marker could not be read. ",
      "Accession retrieval will be checked again."
    )
    
    return(FALSE)
  }
  
  
  # ============================================================
  # Marker must explicitly say completed
  # ============================================================
  
  if (!isTRUE(marker$completed)) {
    return(FALSE)
  }
  
  
  # ============================================================
  # Determine accession manifest path
  # ============================================================
  
  accession_file <- if (
    !is.null(marker$accession_file) &&
    length(marker$accession_file) > 0 &&
    !is.na(marker$accession_file[1]) &&
    nzchar(marker$accession_file[1])
  ) {
    
    marker$accession_file[1]
    
  } else {
    
    "./intermediate_files/all_pulled_accessions.csv"
  }
  
  
  # ============================================================
  # Manifest must actually exist
  # ============================================================
  
  if (!file.exists(accession_file)) {
    
    message(
      "Accession completion marker exists, but the accession manifest ",
      "is missing:\n  ",
      accession_file,
      "\nAccession retrieval will be checked again."
    )
    
    return(FALSE)
  }
  
  
  # ============================================================
  # Database compatibility
  # ============================================================
  
  # Old completion markers predate database tracking.
  # Interpret those as Nucleotide runs.
  marker_database <- if (
    is.null(marker$ncbi_database)
  ) {
    
    "nucleotide"
    
  } else {
    
    normalize_ncbi_database(
      marker$ncbi_database
    )
  }
  
  
  same_database <- identical(
    marker_database,
    ncbi_database
  )
  
  
  if (!same_database) {
    
    message(
      "Existing accession completion marker belongs to NCBI database '",
      marker_database,
      "', not '",
      ncbi_database,
      "'. Accession retrieval will be checked again."
    )
    
    return(FALSE)
  }
  
  
  # ============================================================
  # Search-group compatibility
  # ============================================================
  
  same_taxa <- identical(
    as.character(marker$taxa_list),
    as.character(taxa_list)
  )
  
  
  if (!same_taxa) {
    
    message(
      "Existing accession completion marker belongs to a different ",
      "set of search groups. Accession retrieval will be checked again."
    )
    
    return(FALSE)
  }
  
  
  # ============================================================
  # Validate accession manifest
  #
  # This specifically prevents an empty all_pulled_accessions.csv
  # from being accepted as a completed run.
  # ============================================================
  
  manifest_check <- tryCatch(
    {
      
      manifest <- read.csv(
        accession_file,
        stringsAsFactors = FALSE,
        colClasses = "character"
      )
      
      
      if (
        !"Accession" %in%
        names(manifest)
      ) {
        
        stop(
          "Accession column is missing."
        )
      }
      
      
      valid_accessions <- manifest$Accession[
        !is.na(manifest$Accession) &
          nzchar(trimws(manifest$Accession))
      ]
      
      
      length(
        unique(valid_accessions)
      )
    },
    error = function(e) {
      
      message(
        "Could not validate accession manifest: ",
        conditionMessage(e)
      )
      
      NA_integer_
    }
  )
  
  
  # ============================================================
  # Unreadable manifest = not complete
  # ============================================================
  
  if (is.na(manifest_check)) {
    
    message(
      "Existing accession manifest could not be validated. ",
      "Accession retrieval will be checked again."
    )
    
    return(FALSE)
  }
  
  
  # ============================================================
  # Empty manifest = definitely not complete
  # ============================================================
  
  if (manifest_check == 0) {
    
    message(
      "Accession completion marker points to an empty accession manifest. ",
      "This run will NOT be treated as complete."
    )
    
    return(FALSE)
  }
  
  # ============================================================
  # For newer markers, verify that the manifest still contains
  # the same number of accessions that were present when the
  # completion marker was written.
  # ============================================================
  
  if (
    !is.null(marker$n_accessions)
  ) {
    
    expected_accessions <- as.integer(
      marker$n_accessions
    )
    
    
    if (
      is.na(expected_accessions) ||
      manifest_check != expected_accessions
    ) {
      
      message(
        "Accession completion marker does not match the current manifest."
      )
      
      
      message(
        "Marker count: ",
        expected_accessions,
        " | Current manifest count: ",
        manifest_check
      )
      
      
      message(
        "Accession retrieval will be checked again."
      )
      
      
      return(FALSE)
    }
  }
  
  # ============================================================
  # Existing marker and manifest appear internally consistent
  # ============================================================
  
  TRUE
}


parse_gbseq_node <- function(gbseq_node) {
  # keep list (top-level + qualifier-level)
  if (!exists("metadata_categories_keep", .GlobalEnv)) {
    metadata_categories_keep <- c(
      # top-level
      "GBSeq_locus","GBSeq_length","GBSeq_strandedness","GBSeq_moltype",
      "GBSeq_update-date","GBSeq_create-date","GBSeq_definition",
      "GBSeq_accession-version","GBSeq_project","GBSeq_organism","GBSeq_taxonomy",
      "GBSeq_sequence",
      # qualifiers
      "isolation_source","host","country","lat_lon","collection_date","geo_loc_name",
      "strain","isolate","culture_collection","specimen_voucher",
      "type_material","identified_by","note","gene","product","db_xref"
    )
  }
  
  doc <- gbseq_node
  
  # top-level fields
  top_locus     <- XML::xpathSApply(doc, "./GBSeq_locus", XML::xmlValue)
  top_len       <- XML::xpathSApply(doc, "./GBSeq_length", XML::xmlValue)
  top_strand    <- XML::xpathSApply(doc, "./GBSeq_strandedness", XML::xmlValue)
  top_moltype   <- XML::xpathSApply(doc, "./GBSeq_moltype", XML::xmlValue)
  top_upd       <- XML::xpathSApply(doc, "./GBSeq_update-date", XML::xmlValue)
  top_create    <- XML::xpathSApply(doc, "./GBSeq_create-date", XML::xmlValue)
  top_def       <- XML::xpathSApply(doc, "./GBSeq_definition", XML::xmlValue)
  top_accver    <- XML::xpathSApply(doc, "./GBSeq_accession-version", XML::xmlValue)
  top_proj      <- XML::xpathSApply(doc, "./GBSeq_project", XML::xmlValue)
  top_org       <- XML::xpathSApply(doc, "./GBSeq_organism", XML::xmlValue)
  top_tax       <- XML::xpathSApply(doc, "./GBSeq_taxonomy", XML::xmlValue)
  top_seq       <- XML::xpathSApply(doc, "./GBSeq_sequence", XML::xmlValue)
  
  # qualifiers
  # Parse each qualifier node individually so names and values stay aligned,
  # even when a qualifier has no GBQualifier_value element.
  qualifier_nodes <- XML::getNodeSet(
    doc,
    ".//GBQualifier"
  )
  
  if (length(qualifier_nodes) > 0) {
    
    qualifier_rows <- lapply(
      qualifier_nodes,
      function(node) {
        
        name_node <- XML::getNodeSet(
          node,
          "./GBQualifier_name"
        )
        
        value_node <- XML::getNodeSet(
          node,
          "./GBQualifier_value"
        )
        
        qualifier_name <- if (length(name_node) > 0) {
          XML::xmlValue(name_node[[1]])
        } else {
          NA_character_
        }
        
        qualifier_value <- if (length(value_node) > 0) {
          XML::xmlValue(value_node[[1]])
        } else {
          ""
        }
        
        data.frame(
          name = qualifier_name,
          value = qualifier_value,
          stringsAsFactors = FALSE
        )
      }
    )
    
    quals <- dplyr::bind_rows(qualifier_rows)
    
  } else {
    
    quals <- data.frame(
      name = character(0),
      value = character(0),
      stringsAsFactors = FALSE
    )
  }
  
  # make a named list for qualifiers we care about
  get_q <- function(nm) {
    v <- quals$value[quals$name == nm]
    if (length(v) == 0) "" else paste(unique(v), collapse = "; ")
  }
  
  out <- data.frame(
    Accession          = if (length(top_locus)) top_locus else accession,
    GBSeq_length       = if (length(top_len)) top_len else NA,
    GBSeq_strandedness = if (length(top_strand)) top_strand else NA,
    GBSeq_moltype      = if (length(top_moltype)) top_moltype else NA,
    GBSeq_update.date  = if (length(top_upd)) top_upd else NA,
    GBSeq_create.date  = if (length(top_create)) top_create else NA,
    accession_title    = if (length(top_def)) top_def else NA,
    GBSeq_accession.version = if (length(top_accver)) top_accver else accession,
    GBSeq_project      = if (length(top_proj)) top_proj else NA,
    organism           = if (length(top_org)) top_org else NA,
    GBSeq_taxonomy     = if (length(top_tax)) top_tax else NA,
    sequence           = if (length(top_seq)) top_seq else NA,
    isolation_source   = get_q("isolation_source"),
    host               = get_q("host"),
    geo_loc_name       = get_q("geo_loc_name"),
    country            = get_q("country"),
    lat_lon            = get_q("lat_lon"),
    collection_date    = get_q("collection_date"),
    strain             = get_q("strain"),
    isolate            = get_q("isolate"),
    culture_collection = get_q("culture_collection"),
    specimen_voucher   = get_q("specimen_voucher"),
    type_material      = get_q("type_material"),
    identified_by      = get_q("identified_by"),
    note               = get_q("note"),
    gene               = get_q("gene"),
    product            = get_q("product"),
    db_xref            = get_q("db_xref"),
    stringsAsFactors   = FALSE
  )
  
  out
}


parse_biosample_node <- function(
    biosample_node
) {
  
  # ============================================================
  # Small XML helpers
  # ============================================================
  
  xml_attr <- function(
    node,
    attribute,
    default = NA_character_
  ) {
    
    attrs <- XML::xmlAttrs(
      node
    )
    
    
    if (
      is.null(attrs) ||
      !attribute %in%
      names(attrs)
    ) {
      return(default)
    }
    
    
    value <- as.character(
      attrs[[attribute]]
    )
    
    
    if (
      length(value) == 0 ||
      is.na(value) ||
      !nzchar(value)
    ) {
      return(default)
    }
    
    
    value
  }
  
  
  xml_first_value <- function(
    node,
    path,
    default = NA_character_
  ) {
    
    hits <- XML::getNodeSet(
      node,
      path
    )
    
    
    if (length(hits) == 0) {
      return(default)
    }
    
    
    value <- XML::xmlValue(
      hits[[1]]
    )
    
    
    if (
      is.null(value) ||
      length(value) == 0 ||
      !nzchar(trimws(value))
    ) {
      return(default)
    }
    
    
    trimws(
      as.character(value)
    )
  }
  
  
  normalize_attribute_name <- function(
    x
  ) {
    
    x <- tolower(
      trimws(
        as.character(x)
      )
    )
    
    
    x <- gsub(
      "[^a-z0-9]+",
      "_",
      x
    )
    
    
    x <- gsub(
      "^_+|_+$",
      "",
      x
    )
    
    
    x
  }
  
  
  # ============================================================
  # Basic BioSample information
  # ============================================================
  
  biosample_accession <- xml_attr(
    biosample_node,
    "accession"
  )
  
  
  biosample_id <- xml_attr(
    biosample_node,
    "id"
  )
  
  
  biosample_access <- xml_attr(
    biosample_node,
    "access"
  )
  
  
  submission_date <- xml_attr(
    biosample_node,
    "submission_date"
  )
  
  
  publication_date <- xml_attr(
    biosample_node,
    "publication_date"
  )
  
  
  last_update <- xml_attr(
    biosample_node,
    "last_update"
  )
  
  
  title <- xml_first_value(
    biosample_node,
    "./Description/Title"
  )
  
  
  owner <- xml_first_value(
    biosample_node,
    "./Owner/Name"
  )
  
  
  package <- xml_first_value(
    biosample_node,
    "./Package"
  )
  
  
  description_nodes <- XML::getNodeSet(
    biosample_node,
    "./Description/Comment/Paragraph"
  )
  
  
  biosample_description <- if (
    length(description_nodes) > 0
  ) {
    
    description_values <- vapply(
      description_nodes,
      XML::xmlValue,
      character(1)
    )
    
    
    description_values <- unique(
      trimws(
        description_values[
          nzchar(
            trimws(
              description_values
            )
          )
        ]
      )
    )
    
    
    if (
      length(description_values) > 0
    ) {
      
      paste(
        description_values,
        collapse = "; "
      )
      
    } else {
      NA_character_
    }
    
  } else {
    NA_character_
  }
  
  
  # ============================================================
  # Organism
  # ============================================================
  
  organism_nodes <- XML::getNodeSet(
    biosample_node,
    "./Description/Organism"
  )
  
  
  if (
    length(organism_nodes) > 0
  ) {
    
    organism_node <-
      organism_nodes[[1]]
    
    
    organism <- xml_attr(
      organism_node,
      "taxonomy_name"
    )
    
    
    taxid <- xml_attr(
      organism_node,
      "taxonomy_id"
    )
    
  } else {
    
    organism <- NA_character_
    taxid <- NA_character_
  }
  
  
  # ============================================================
  # Models
  # ============================================================
  
  model_nodes <- XML::getNodeSet(
    biosample_node,
    "./Models/Model"
  )
  
  
  biosample_models <- if (
    length(model_nodes) > 0
  ) {
    
    model_values <- vapply(
      model_nodes,
      XML::xmlValue,
      character(1)
    )
    
    
    model_values <- unique(
      trimws(
        model_values[
          nzchar(
            trimws(
              model_values
            )
          )
        ]
      )
    )
    
    
    paste(
      model_values,
      collapse = "; "
    )
    
  } else {
    NA_character_
  }
  
  
  # ============================================================
  # Other identifiers
  # ============================================================
  
  id_nodes <- XML::getNodeSet(
    biosample_node,
    "./Ids/Id"
  )
  
  
  biosample_other_ids <- if (
    length(id_nodes) > 0
  ) {
    
    id_pairs <- vapply(
      id_nodes,
      function(node) {
        
        db_name <- xml_attr(
          node,
          "db",
          default = ""
        )
        
        
        id_value <- trimws(
          XML::xmlValue(
            node
          )
        )
        
        
        if (nzchar(db_name)) {
          
          paste0(
            db_name,
            "=",
            id_value
          )
          
        } else {
          id_value
        }
      },
      character(1)
    )
    
    
    id_pairs <- unique(
      id_pairs[
        nzchar(id_pairs)
      ]
    )
    
    
    paste(
      id_pairs,
      collapse = "; "
    )
    
  } else {
    NA_character_
  }
  
  
  # ============================================================
  # BioSample attributes
  # ============================================================
  
  attribute_nodes <- XML::getNodeSet(
    biosample_node,
    "./Attributes/Attribute"
  )
  
  
  if (
    length(attribute_nodes) > 0
  ) {
    
    attribute_rows <- lapply(
      attribute_nodes,
      function(node) {
        
        original_name <- xml_attr(
          node,
          "attribute_name",
          default = ""
        )
        
        
        harmonized_name <- xml_attr(
          node,
          "harmonized_name",
          default = ""
        )
        
        
        value <- trimws(
          XML::xmlValue(
            node
          )
        )
        
        
        preferred_name <- if (
          nzchar(harmonized_name)
        ) {
          harmonized_name
        } else {
          original_name
        }
        
        
        data.frame(
          original_name =
            original_name,
          harmonized_name =
            harmonized_name,
          preferred_name =
            preferred_name,
          normalized_original =
            normalize_attribute_name(
              original_name
            ),
          normalized_preferred =
            normalize_attribute_name(
              preferred_name
            ),
          value =
            value,
          stringsAsFactors = FALSE
        )
      }
    )
    
    
    attributes <- dplyr::bind_rows(
      attribute_rows
    )
    
    
  } else {
    
    attributes <- data.frame(
      original_name = character(0),
      harmonized_name = character(0),
      preferred_name = character(0),
      normalized_original = character(0),
      normalized_preferred = character(0),
      value = character(0),
      stringsAsFactors = FALSE
    )
  }
  
  
  # ============================================================
  # Get one or more BioSample attribute values
  # ============================================================
  
  get_biosample_attribute <- function(
    names_to_find
  ) {
    
    if (
      nrow(attributes) == 0
    ) {
      return("")
    }
    
    
    keys <- normalize_attribute_name(
      names_to_find
    )
    
    
    hit <- (
      attributes$normalized_preferred %in%
        keys
    ) |
      (
        attributes$normalized_original %in%
          keys
      )
    
    
    values <- attributes$value[
      hit
    ]
    
    
    values <- unique(
      trimws(
        values[
          !is.na(values) &
            nzchar(
              trimws(values)
            )
        ]
      )
    )
    
    # BioSample placeholder values should be treated as missing data
    missing_values <- c(
      "missing",
      "not provided",
      "not collected",
      "not applicable",
      "not available",
      "unknown",
      "na",
      "n/a"
    )
    
    values <- values[
      !tolower(values) %in% missing_values
    ]
    
    
    if (length(values) == 0) {
      return("")
    }
    
    
    paste(
      values,
      collapse = "; "
    )
  }
  
  
  # ============================================================
  # Preserve every attribute in one human-readable column
  # ============================================================
  
  all_attributes <- if (
    nrow(attributes) > 0
  ) {
    
    attribute_pairs <- paste0(
      attributes$preferred_name,
      "=",
      attributes$value
    )
    
    
    attribute_pairs <- unique(
      attribute_pairs[
        nzchar(
          attributes$preferred_name
        )
      ]
    )
    
    
    paste(
      attribute_pairs,
      collapse = " | "
    )
    
  } else {
    ""
  }
  
  
  # ============================================================
  # Flat aRborist-compatible row
  # ============================================================
  
  out <- data.frame(
    Accession =
      biosample_accession,
    
    EntrezUID =
      biosample_id,
    
    ncbi_database =
      "biosample",
    
    accession_title =
      title,
    
    organism =
      organism,
    
    TaxID =
      taxid,
    
    BioSample_access =
      biosample_access,
    
    BioSample_submission_date =
      submission_date,
    
    BioSample_publication_date =
      publication_date,
    
    BioSample_last_update =
      last_update,
    
    BioSample_owner =
      owner,
    
    BioSample_package =
      package,
    
    BioSample_models =
      biosample_models,
    
    BioSample_description =
      biosample_description,
    
    BioSample_other_ids =
      biosample_other_ids,
    
    sample_name =
      get_biosample_attribute(
        "sample_name"
      ),
    
    sample_type =
      get_biosample_attribute(
        "sample_type"
      ),
    
    strain =
      get_biosample_attribute(
        "strain"
      ),
    
    isolate =
      get_biosample_attribute(
        "isolate"
      ),
    
    host =
      get_biosample_attribute(
        "host"
      ),
    
    host_taxid =
      get_biosample_attribute(
        "host_taxid"
      ),
    
    host_disease =
      get_biosample_attribute(
        "host_disease"
      ),
    
    isolation_source =
      get_biosample_attribute(
        "isolation_source"
      ),
    
    collection_date =
      get_biosample_attribute(
        "collection_date"
      ),
    
    geo_loc_name =
      get_biosample_attribute(
        c(
          "geo_loc_name"
        )
      ),
    
    country =
      get_biosample_attribute(
        "country"
      ),
    
    lat_lon =
      get_biosample_attribute(
        "lat_lon"
      ),
    
    tissue =
      get_biosample_attribute(
        "tissue"
      ),
    
    culture_collection =
      get_biosample_attribute(
        "culture_collection"
      ),
    
    specimen_voucher =
      get_biosample_attribute(
        c(
          "specimen_voucher",
          "specimen voucher"
        )
      ),
    
    type_material =
      get_biosample_attribute(
        "type_material"
      ),
    
    identified_by =
      get_biosample_attribute(
        "identified_by"
      ),
    
    BioSample_all_attributes =
      all_attributes,
    
    stringsAsFactors = FALSE
  )
  
  
  out
}

retry_failed_metadata_accessions <- function(
    checkpoint_dir = "./metadata_files/metadata_checkpoints",
    metadata_batch_size = 250,
    max_retry_passes = 2,
    ncbi_database = get0(
      "ncbi_database",
      envir = .GlobalEnv,
      ifnotfound = default_ncbi_database
    )
) {
  
  ncbi_database <- normalize_ncbi_database(
    ncbi_database
  )
  
  if (
    ncbi_database == "biosample" &&
    identical(
      checkpoint_dir,
      "./metadata_files/metadata_checkpoints"
    )
  ) {
    
    checkpoint_dir <-
      "./metadata_files/metadata_checkpoints_biosample"
  }
  
  failed_path <- file.path(
    checkpoint_dir,
    "metadata_failed_accessions.csv"
  )
  
  recovered_path <- file.path(
    checkpoint_dir,
    "metadata_failed_accessions_recovered.csv"
  )
  
  normalize_accession <- function(x) {
    x <- trimws(as.character(x))
    sub("\\.[0-9]+$", "", x)
  }
  
  
  # ------------------------------------------------------------
  # Nothing to retry
  # ------------------------------------------------------------
  
  if (!file.exists(failed_path)) {
    message("\nNo failed-accession file found. No metadata retries needed.")
    return(invisible(NULL))
  }
  
  failed_df <- read.csv(
    failed_path,
    stringsAsFactors = FALSE,
    colClasses = "character"
  )
  
  if (
    nrow(failed_df) == 0 ||
    !"Accession" %in% names(failed_df)
  ) {
    message("\nNo failed accessions remain. No metadata retries needed.")
    return(invisible(NULL))
  }
  
  
  # ------------------------------------------------------------
  # Collapse duplicate failure entries
  # ------------------------------------------------------------
  
  failed_df$Accession_normalized <- normalize_accession(
    failed_df$Accession
  )
  
  failed_df <- failed_df[
    !is.na(failed_df$Accession_normalized) &
      nzchar(failed_df$Accession_normalized),
    ,
    drop = FALSE
  ]
  
  failed_df <- failed_df[
    !duplicated(failed_df$Accession_normalized),
    ,
    drop = FALSE
  ]
  
  
  # ------------------------------------------------------------
  # Determine which failed accessions may already exist in
  # completed metadata checkpoints
  # ------------------------------------------------------------
  
  checkpoint_files <- list.files(
    checkpoint_dir,
    pattern = "^metadata_.*\\.csv$",
    full.names = TRUE
  )
  
  checkpoint_files <- checkpoint_files[
    !grepl(
      "metadata_failed_accessions\\.csv$",
      checkpoint_files
    )
  ]
  
  already_recovered <- character(0)
  
  for (checkpoint_file in checkpoint_files) {
    
    this_checkpoint <- tryCatch(
      {
        read.csv(
          checkpoint_file,
          stringsAsFactors = FALSE,
          colClasses = "character"
        )
      },
      error = function(e) {
        NULL
      }
    )
    
    if (
      is.null(this_checkpoint) ||
      nrow(this_checkpoint) == 0
    ) {
      next
    }
    
    accession_col <- if (
      "GBSeq_accession.version" %in% names(this_checkpoint)
    ) {
      "GBSeq_accession.version"
    } else if (
      "Accession" %in% names(this_checkpoint)
    ) {
      "Accession"
    } else {
      NULL
    }
    
    if (!is.null(accession_col)) {
      
      already_recovered <- c(
        already_recovered,
        normalize_accession(
          this_checkpoint[[accession_col]]
        )
      )
    }
  }
  
  already_recovered <- unique(
    already_recovered[
      !is.na(already_recovered) &
        nzchar(already_recovered)
    ]
  )
  
  if (length(already_recovered) > 0) {
    
    before_n <- nrow(failed_df)
    
    failed_df <- failed_df[
      !failed_df$Accession_normalized %in% already_recovered,
      ,
      drop = FALSE
    ]
    
    removed_n <- before_n - nrow(failed_df)
    
    if (removed_n > 0) {
      message(
        "\nRemoved ",
        removed_n,
        " accession(s) from the failure list because metadata already exists."
      )
    }
  }
  
  
  # ------------------------------------------------------------
  # If everything was already recovered
  # ------------------------------------------------------------
  
  if (nrow(failed_df) == 0) {
    
    empty_failed <- data.frame(
      Taxon = character(),
      Accession = character(),
      Error = character(),
      Time = character(),
      stringsAsFactors = FALSE
    )
    
    write.csv(
      empty_failed,
      failed_path,
      row.names = FALSE
    )
    
    message("\nAll previously failed accessions already have metadata.")
    
    return(invisible(NULL))
  }
  
  
  # ------------------------------------------------------------
  # Recursive retrieval helper
  #
  # Try the whole batch first.
  # If it errors, divide the batch in half.
  # Continue until successful or down to one accession.
  # ------------------------------------------------------------
  
  fetch_with_fallback <- function(accessions) {
    
    accessions <- unique(
      trimws(as.character(accessions))
    )
    
    accessions <- accessions[
      !is.na(accessions) &
        nzchar(accessions)
    ]
    
    if (length(accessions) == 0) {
      
      return(
        list(
          metadata = NULL,
          failed = data.frame(
            Accession = character(),
            Error = character(),
            stringsAsFactors = FALSE
          )
        )
      )
    }
    
    
    # ----------------------------------------------------------
    # Individual accession
    # ----------------------------------------------------------
    
    if (length(accessions) == 1) {
      
      acc <- accessions[1]
      
      entry <- tryCatch(
        {
          fetch_metadata_for_accession(
            acc,
            ncbi_database = ncbi_database
          )
        },
        error = function(e) {
          
          return(
            structure(
              NULL,
              retrieval_error = conditionMessage(e)
            )
          )
        }
      )
      
      if (is.null(entry)) {
        
        err <- attr(
          entry,
          "retrieval_error"
        )
        
        if (
          is.null(err) ||
          !nzchar(err)
        ) {
          err <- "Metadata retrieval failed."
        }
        
        return(
          list(
            metadata = NULL,
            failed = data.frame(
              Accession = acc,
              Error = err,
              stringsAsFactors = FALSE
            )
          )
        )
      }
      
      return(
        list(
          metadata = entry,
          failed = data.frame(
            Accession = character(),
            Error = character(),
            stringsAsFactors = FALSE
          )
        )
      )
    }
    
    
    # ----------------------------------------------------------
    # Attempt batch
    # ----------------------------------------------------------
    
    batch_result <- tryCatch(
      {
        fetch_metadata_for_accession_batch(
          accessions,
          ncbi_database = ncbi_database
        )
      },
      error = function(e) {
        NULL
      }
    )
    
    
    # ----------------------------------------------------------
    # Batch succeeded
    # ----------------------------------------------------------
    
    if (
      !is.null(batch_result) &&
      nrow(batch_result) > 0
    ) {
      
      returned_col <- if (
        "GBSeq_accession.version" %in% names(batch_result)
      ) {
        "GBSeq_accession.version"
      } else {
        "Accession"
      }
      
      returned <- normalize_accession(
        batch_result[[returned_col]]
      )
      
      requested <- normalize_accession(
        accessions
      )
      
      missing <- accessions[
        !requested %in% returned
      ]
      
      
      # Everything returned
      if (length(missing) == 0) {
        
        return(
          list(
            metadata = batch_result,
            failed = data.frame(
              Accession = character(),
              Error = character(),
              stringsAsFactors = FALSE
            )
          )
        )
      }
      
      
      # Some records were omitted by NCBI.
      # Retry only the omitted accessions.
      missing_result <- fetch_with_fallback(
        missing
      )
      
      combined_metadata <- batch_result
      
      if (!is.null(missing_result$metadata)) {
        
        combined_metadata <- plyr::rbind.fill(
          list(
            combined_metadata,
            missing_result$metadata
          )
        )
      }
      
      return(
        list(
          metadata = combined_metadata,
          failed = missing_result$failed
        )
      )
    }
    
    
    # ----------------------------------------------------------
    # Batch failed completely:
    # divide it and retry smaller batches
    # ----------------------------------------------------------
    
    midpoint <- floor(
      length(accessions) / 2
    )
    
    left_accessions <- accessions[
      seq_len(midpoint)
    ]
    
    right_accessions <- accessions[
      (midpoint + 1):length(accessions)
    ]
    
    left_result <- fetch_with_fallback(
      left_accessions
    )
    
    Sys.sleep(
      get_sleep_duration()
    )
    
    right_result <- fetch_with_fallback(
      right_accessions
    )
    
    
    combined_metadata <- NULL
    
    metadata_parts <- list()
    
    if (!is.null(left_result$metadata)) {
      metadata_parts[[length(metadata_parts) + 1L]] <-
        left_result$metadata
    }
    
    if (!is.null(right_result$metadata)) {
      metadata_parts[[length(metadata_parts) + 1L]] <-
        right_result$metadata
    }
    
    if (length(metadata_parts) > 0) {
      combined_metadata <- plyr::rbind.fill(
        metadata_parts
      )
    }
    
    combined_failed <- rbind(
      left_result$failed,
      right_result$failed
    )
    
    list(
      metadata = combined_metadata,
      failed = combined_failed
    )
  }
  
  
  # ============================================================
  # RETRY PASSES
  # ============================================================
  
  for (retry_pass in seq_len(max_retry_passes)) {
    
    if (nrow(failed_df) == 0) {
      break
    }
    
    message(
      "\n============================================================"
    )
    
    message(
      "FAILED METADATA RETRY PASS ",
      retry_pass,
      " / ",
      max_retry_passes
    )
    
    message(
      "Retrying ",
      nrow(failed_df),
      " unique accession(s)."
    )
    
    message(
      "============================================================\n"
    )
    
    retry_accessions <- failed_df$Accession
    
    retry_batches <- split(
      retry_accessions,
      ceiling(
        seq_along(retry_accessions) /
          metadata_batch_size
      )
    )
    
    recovered_this_pass <- list()
    failures_this_pass <- list()
    
    num_batches <- length(
      retry_batches
    )
    
    
    for (batch_index in seq_along(retry_batches)) {
      
      accession_batch <- unname(
        retry_batches[[batch_index]]
      )
      
      message(
        "Retry pass ",
        retry_pass,
        " | batch ",
        batch_index,
        " / ",
        num_batches,
        " | ",
        length(accession_batch),
        " accession(s)"
      )
      
      retry_result <- fetch_with_fallback(
        accession_batch
      )
      
      
      if (!is.null(retry_result$metadata)) {
        
        recovered_this_pass[[
          length(recovered_this_pass) + 1L
        ]] <- retry_result$metadata
      }
      
      
      if (
        !is.null(retry_result$failed) &&
        nrow(retry_result$failed) > 0
      ) {
        
        failures_this_pass[[
          length(failures_this_pass) + 1L
        ]] <- retry_result$failed
      }
      
      Sys.sleep(
        get_sleep_duration()
      )
    }
    
    
    # ----------------------------------------------------------
    # Save recovered records
    # ----------------------------------------------------------
    
    recovered_df <- NULL
    
    if (length(recovered_this_pass) > 0) {
      
      recovered_df <- plyr::rbind.fill(
        recovered_this_pass
      )
      
      accession_col <- if (
        "GBSeq_accession.version" %in% names(recovered_df)
      ) {
        "GBSeq_accession.version"
      } else {
        "Accession"
      }
      
      recovered_df$.normalized_accession <-
        normalize_accession(
          recovered_df[[accession_col]]
        )
      
      recovered_df <- recovered_df[
        !duplicated(
          recovered_df$.normalized_accession
        ),
        ,
        drop = FALSE
      ]
      
      recovered_df$.normalized_accession <- NULL
      
      
      if (file.exists(recovered_path)) {
        
        existing_recovered <- read.csv(
          recovered_path,
          stringsAsFactors = FALSE,
          colClasses = "character"
        )
        
        recovered_df <- plyr::rbind.fill(
          list(
            existing_recovered,
            recovered_df
          )
        )
        
        accession_col <- if (
          "GBSeq_accession.version" %in% names(recovered_df)
        ) {
          "GBSeq_accession.version"
        } else {
          "Accession"
        }
        
        recovered_df$.normalized_accession <-
          normalize_accession(
            recovered_df[[accession_col]]
          )
        
        recovered_df <- recovered_df[
          !duplicated(
            recovered_df$.normalized_accession
          ),
          ,
          drop = FALSE
        ]
        
        recovered_df$.normalized_accession <- NULL
      }
      
      
      write.csv(
        recovered_df,
        recovered_path,
        row.names = FALSE
      )
    }
    
    
    # ----------------------------------------------------------
    # Determine what still failed
    # ----------------------------------------------------------
    
    if (length(failures_this_pass) > 0) {
      
      remaining_failures <- do.call(
        rbind,
        failures_this_pass
      )
      
      remaining_failures$Accession_normalized <-
        normalize_accession(
          remaining_failures$Accession
        )
      
      remaining_failures <- remaining_failures[
        !duplicated(
          remaining_failures$Accession_normalized
        ),
        ,
        drop = FALSE
      ]
      
      
      original_taxon_map <- failed_df[
        !duplicated(
          failed_df$Accession_normalized
        ),
        c(
          "Accession_normalized",
          "Taxon"
        ),
        drop = FALSE
      ]
      
      remaining_failures <- merge(
        remaining_failures,
        original_taxon_map,
        by = "Accession_normalized",
        all.x = TRUE,
        sort = FALSE
      )
      
      failed_df <- data.frame(
        Taxon = remaining_failures$Taxon,
        Accession = remaining_failures$Accession,
        Error = remaining_failures$Error,
        Time = format(
          Sys.time(),
          "%Y-%m-%d %H:%M:%S"
        ),
        Accession_normalized =
          remaining_failures$Accession_normalized,
        stringsAsFactors = FALSE
      )
      
    } else {
      
      failed_df <- data.frame(
        Taxon = character(),
        Accession = character(),
        Error = character(),
        Time = character(),
        Accession_normalized = character(),
        stringsAsFactors = FALSE
      )
    }
    
    
    # ----------------------------------------------------------
    # Rewrite the failure log after every pass
    # ----------------------------------------------------------
    
    failed_to_write <- failed_df
    
    failed_to_write$Accession_normalized <- NULL
    
    write.csv(
      failed_to_write,
      failed_path,
      row.names = FALSE
    )
    
    
    recovered_count <- length(
      retry_accessions
    ) - nrow(
      failed_df
    )
    
    message(
      "\nRetry pass ",
      retry_pass,
      " complete."
    )
    
    message(
      "  Recovered: ",
      recovered_count
    )
    
    message(
      "  Still failed: ",
      nrow(failed_df)
    )
    
    
    if (nrow(failed_df) == 0) {
      
      message(
        "\nAll failed metadata accessions were recovered."
      )
      
      break
    }
  }
  
  
  # ============================================================
  # FINAL STATUS
  # ============================================================
  
  if (nrow(failed_df) > 0) {
    
    message(
      "\nMetadata retry limit reached."
    )
    
    message(
      nrow(failed_df),
      " accession(s) still failed after ",
      max_retry_passes,
      " retry pass(es)."
    )
    
    message(
      "No additional automatic retries will be attempted."
    )
    
  } else {
    
    message(
      "\nNo unresolved failed metadata accessions remain."
    )
  }
  
  
  invisible(
    failed_df
  )
}

# Metadata retrieval, using list of pulled accessions
fetch_metadata_for_accession <- function(
    accession,
    ncbi_database = get0(
      "ncbi_database",
      envir = .GlobalEnv,
      ifnotfound = default_ncbi_database
    )
) {
  
  result <- fetch_metadata_for_accession_batch(
    accessions = accession,
    ncbi_database = ncbi_database
  )
  
  
  if (
    is.null(result) ||
    nrow(result) == 0
  ) {
    return(NULL)
  }
  
  
  result[
    1,
    ,
    drop = FALSE
  ]
}

# metadata retrieval
retrieve_ncbi_metadata <- function(
    project_name,
    resume = TRUE,
    overwrite_checkpoints = FALSE,
    checkpoint_dir = "./metadata_files/metadata_checkpoints",
    ncbi_database = get0(
      "ncbi_database",
      envir = .GlobalEnv,
      ifnotfound = default_ncbi_database
    ),
    checkpoint_every = 500,
    progress_every = 50,
    metadata_batch_size = 250,
    batch_max_retries = 2,
    batch_retry_wait = 5,
    min_batch_size = 1,
    failed_retry_passes = 2
) {
  
  # ============================================================
  # Resolve NCBI database
  # ============================================================
  
  ncbi_database <- normalize_ncbi_database(
    ncbi_database
  )
  
  
  # Keep the existing metadata checkpoint directory for
  # Nucleotide searches, but use a separate directory for
  # BioSample searches so the two databases cannot accidentally
  # reuse each other's checkpoints.
  if (
    ncbi_database == "biosample" &&
    identical(
      checkpoint_dir,
      "./metadata_files/metadata_checkpoints"
    )
  ) {
    
    checkpoint_dir <-
      "./metadata_files/metadata_checkpoints_biosample"
  }
  
  
  message(
    "\nNCBI metadata source: ",
    ncbi_database
  )
  
  accession_path <- "./intermediate_files/all_pulled_accessions.csv"
  
  if (!file.exists(accession_path)) {
    stop("Accession list not found at: ", accession_path)
  }
  
  if (!dir.exists("./metadata_files")) {
    dir.create("./metadata_files", recursive = TRUE)
  }
  
  if (!dir.exists(checkpoint_dir)) {
    dir.create(checkpoint_dir, recursive = TRUE)
  }
  
  accession_list <- read.csv(
    accession_path,
    header = TRUE,
    stringsAsFactors = FALSE
  )
  
  if (!"Accession" %in% names(accession_list)) {
    stop("Accession list must contain an 'Accession' column.")
  }
  
  # ============================================================
  # Normalize accession versions for matching
  # ============================================================
  
  normalize_accession <- function(x) {
    x <- trimws(as.character(x))
    sub("\\.[0-9]+$", "", x)
  }
  
  # ============================================================
  # Split metadata retrieval into search groups
  # ============================================================
  
  if ("search_group" %in% names(accession_list)) {
    
    taxa_groups <- split(
      accession_list,
      accession_list$search_group
    )
    
  } else if ("genus" %in% names(accession_list)) {
    
    message(
      "Using legacy 'genus' column as the metadata checkpoint group."
    )
    
    taxa_groups <- split(
      accession_list,
      accession_list$genus
    )
    
  } else {
    
    taxa_groups <- list(
      ALL = accession_list
    )
  }
  
  # ============================================================
  # Output / logging paths
  # ============================================================
  
  failed_path <- file.path(
    checkpoint_dir,
    "metadata_failed_accessions.csv"
  )
  
  recovered_failed_path <- file.path(
    checkpoint_dir,
    "metadata_failed_accessions_recovered.csv"
  )
  
  timing_path <- "./intermediate_files/fetch_times_metadata_by_taxon.csv"
  
  # ------------------------------------------------------------
  # Preserve an existing failure log when resuming
  #
  # This is important when all normal search-group checkpoints
  # already exist and we only want to retry prior failures.
  # ------------------------------------------------------------
  
  if (
    resume &&
    !overwrite_checkpoints &&
    file.exists(failed_path)
  ) {
    
    failed_log <- read.csv(
      failed_path,
      stringsAsFactors = FALSE,
      colClasses = "character"
    )
    
    required_failed_columns <- c(
      "Taxon",
      "Accession",
      "Error",
      "Time"
    )
    
    missing_failed_columns <- setdiff(
      required_failed_columns,
      names(failed_log)
    )
    
    if (length(missing_failed_columns) > 0) {
      stop(
        "Existing failed accession file is missing column(s): ",
        paste(
          missing_failed_columns,
          collapse = ", "
        )
      )
    }
    
    failed_log <- failed_log[
      !is.na(failed_log$Accession) &
        nzchar(failed_log$Accession),
      ,
      drop = FALSE
    ]
    
    if (nrow(failed_log) > 0) {
      
      failed_log$.normalized_accession <- normalize_accession(
        failed_log$Accession
      )
      
      failed_log <- failed_log[
        !duplicated(failed_log$.normalized_accession),
        ,
        drop = FALSE
      ]
      
      failed_log$.normalized_accession <- NULL
    }
    
    message(
      "\nLoaded ",
      nrow(failed_log),
      " previously failed accession(s) for possible retry."
    )
    
  } else {
    
    failed_log <- data.frame(
      Taxon = character(),
      Accession = character(),
      Error = character(),
      Time = character(),
      stringsAsFactors = FALSE
    )
  }
  
  timing_log <- data.frame(
    Taxon = character(),
    Num_accessions = integer(),
    Num_successful = integer(),
    Num_failed = integer(),
    Start_time = character(),
    End_time = character(),
    Elapsed_minutes = numeric(),
    Checkpoint_file = character(),
    stringsAsFactors = FALSE
  )
  
  overall_start <- Sys.time()
  
  message(
    "\nStarting metadata retrieval for ",
    nrow(accession_list),
    " accessions across ",
    length(taxa_groups),
    " taxon group(s)."
  )
  
  # ============================================================
  # Process each search group
  # ============================================================
  
  for (tx in names(taxa_groups)) {
    
    safe_tx <- gsub(
      "[^A-Za-z0-9_.-]+",
      "_",
      tx
    )
    
    checkpoint_path <- file.path(
      checkpoint_dir,
      paste0("metadata_", safe_tx, ".csv")
    )
    
    partial_checkpoint_path <- file.path(
      checkpoint_dir,
      paste0("metadata_", safe_tx, ".partial.csv")
    )
    
    block <- taxa_groups[[tx]]
    
    # ----------------------------------------------------------
    # Optional checkpoint overwrite
    # ----------------------------------------------------------
    
    if (overwrite_checkpoints) {
      
      if (file.exists(checkpoint_path)) {
        message(
          "\n--- Removing completed checkpoint for ",
          tx,
          " ---"
        )
        
        file.remove(checkpoint_path)
      }
      
      if (file.exists(partial_checkpoint_path)) {
        message(
          "\n--- Removing partial checkpoint for ",
          tx,
          " ---"
        )
        
        file.remove(partial_checkpoint_path)
      }
    }
    
    # ----------------------------------------------------------
    # Skip completed groups
    # ----------------------------------------------------------
    
    if (
      file.exists(checkpoint_path) &&
      resume &&
      !overwrite_checkpoints
    ) {
      
      message(
        "\n--- Skipping ",
        tx,
        ": completed checkpoint already exists ---"
      )
      
      next
    }
    
    # ----------------------------------------------------------
    # Resume partial metadata checkpoint
    # ----------------------------------------------------------
    
    completed_accessions <- character(0)
    
    if (
      file.exists(partial_checkpoint_path) &&
      resume &&
      !overwrite_checkpoints
    ) {
      
      partial_existing <- read.csv(
        partial_checkpoint_path,
        stringsAsFactors = FALSE,
        colClasses = "character"
      )
      
      if (!"Accession" %in% names(partial_existing)) {
        stop(
          "Partial metadata checkpoint is missing the Accession column: ",
          partial_checkpoint_path
        )
      }
      
      checkpoint_accession_col <- if (
        "GBSeq_accession.version" %in% names(partial_existing)
      ) {
        "GBSeq_accession.version"
      } else {
        "Accession"
      }
      
      completed_accessions <- unique(
        partial_existing[[checkpoint_accession_col]][
          !is.na(partial_existing[[checkpoint_accession_col]]) &
            partial_existing[[checkpoint_accession_col]] != ""
        ]
      )
      
      message(
        "\n--- Resuming ",
        tx,
        " from partial checkpoint with ",
        length(completed_accessions),
        " completed accession(s) ---"
      )
    }
    
    # ----------------------------------------------------------
    # Match partial checkpoint accessions
    # ----------------------------------------------------------
    
    if (length(completed_accessions) > 0) {
      
      completed_accessions_normalized <- normalize_accession(
        completed_accessions
      )
      
      block_accessions_normalized <- normalize_accession(
        block$Accession
      )
      
      already_completed <- block_accessions_normalized %in%
        completed_accessions_normalized
      
      message(
        "Matched ",
        sum(already_completed),
        " of ",
        length(completed_accessions),
        " completed accession(s) to the current download manifest."
      )
      
      block <- block[
        !already_completed,
        ,
        drop = FALSE
      ]
    }
    
    # ----------------------------------------------------------
    # Start group
    # ----------------------------------------------------------
    
    taxon_start <- Sys.time()
    
    message(
      "\n--- ",
      tx,
      ": ",
      nrow(block),
      " accession(s) remaining ---"
    )
    
    metadata_buffer <- list()
    buffer_success_count <- 0L
    
    success_count <- length(completed_accessions)
    fail_count <- 0L
    
    # ----------------------------------------------------------
    # Split remaining accessions into main metadata batches
    # ----------------------------------------------------------
    
    accession_batches <- split(
      block$Accession,
      ceiling(
        seq_along(block$Accession) /
          metadata_batch_size
      )
    )
    
    num_batches <- length(accession_batches)
    processed_count <- 0L
    
    # ==========================================================
    # Metadata batch loop
    # ==========================================================
    
    for (batch_index in seq_along(accession_batches)) {
      
      accession_batch <- unname(
        accession_batches[[batch_index]]
      )
      
      batch_start_position <- processed_count + 1L
      
      batch_end_position <- processed_count +
        length(accession_batch)
      
      # --------------------------------------------------------
      # Progress reporting
      # --------------------------------------------------------
      
      if (
        batch_index == 1 ||
        batch_index %% progress_every == 0 ||
        batch_index == num_batches
      ) {
        
        message(
          "[",
          format_progress(
            current = batch_end_position,
            total = nrow(block),
            start_time = taxon_start
          ),
          "] Fetching batch ",
          batch_index,
          " / ",
          num_batches,
          " (",
          length(accession_batch),
          " accession(s))"
        )
      }
      
      # --------------------------------------------------------
      # Resilient metadata retrieval
      # --------------------------------------------------------
      
      batch_result <- fetch_metadata_batch_resilient(
        accessions = accession_batch,
        min_batch_size = min_batch_size,
        max_retries = batch_max_retries,
        retry_wait = batch_retry_wait,
        ncbi_database = ncbi_database
      )
      
      # --------------------------------------------------------
      # Add successful metadata rows to buffer
      # --------------------------------------------------------
      
      if (
        !is.null(batch_result) &&
        nrow(batch_result) > 0
      ) {
        
        for (row_index in seq_len(nrow(batch_result))) {
          
          entry <- batch_result[
            row_index,
            ,
            drop = FALSE
          ]
          
          accession_display <- if (
            "GBSeq_accession.version" %in% names(entry) &&
            !is.na(entry$GBSeq_accession.version[1]) &&
            nzchar(entry$GBSeq_accession.version[1])
          ) {
            entry$GBSeq_accession.version[1]
          } else {
            entry$Accession[1]
          }
          
          success_count <- success_count + 1L
          buffer_success_count <- buffer_success_count + 1L
          
          metadata_buffer[[length(metadata_buffer) + 1L]] <- entry
          
          message(
            "  OK | ",
            accession_display,
            " | Species: ",
            entry$organism[1],
            " | Strain: ",
            dplyr::coalesce(
              entry$strain[1],
              entry$specimen_voucher[1],
              entry$isolate[1],
              ""
            ),
            " | Host: ",
            entry$host[1]
          )
        }
      }
      
      # --------------------------------------------------------
      # Determine accessions that ultimately were not returned
      # --------------------------------------------------------
      
      returned_accessions <- character(0)
      
      if (
        !is.null(batch_result) &&
        nrow(batch_result) > 0
      ) {
        
        returned_column <- if (
          "GBSeq_accession.version" %in%
          names(batch_result)
        ) {
          "GBSeq_accession.version"
        } else {
          "Accession"
        }
        
        returned_accessions <- normalize_accession(
          batch_result[[returned_column]]
        )
      }
      
      missing_accessions <- accession_batch[
        !normalize_accession(accession_batch) %in%
          returned_accessions
      ]
      
      # --------------------------------------------------------
      # Log only accessions that failed even after progressive
      # batch subdivision
      # --------------------------------------------------------
      
      if (length(missing_accessions) > 0) {
        
        fail_count <- fail_count +
          length(missing_accessions)
        
        missing_rows <- data.frame(
          Taxon = rep(
            tx,
            length(missing_accessions)
          ),
          Accession = missing_accessions,
          Error = rep(
            paste0(
              "Accession was not returned after ",
              "progressive metadata batch retries."
            ),
            length(missing_accessions)
          ),
          Time = rep(
            format(
              Sys.time(),
              "%Y-%m-%d %H:%M:%S"
            ),
            length(missing_accessions)
          ),
          stringsAsFactors = FALSE
        )
        
        failed_log <- rbind(
          failed_log,
          missing_rows
        )
        
        # Deduplicate the failure log as we go
        failed_log$.normalized_accession <- normalize_accession(
          failed_log$Accession
        )
        
        failed_log <- failed_log[
          !duplicated(
            failed_log$.normalized_accession,
            fromLast = TRUE
          ),
          ,
          drop = FALSE
        ]
        
        failed_log$.normalized_accession <- NULL
        
        message(
          "  WARNING: ",
          length(missing_accessions),
          " accession(s) were not returned after ",
          "progressive batch retries."
        )
        
        if (length(missing_accessions) <= 10) {
          message(
            "    ",
            paste(
              missing_accessions,
              collapse = ", "
            )
          )
        }
      }
      
      processed_count <- batch_end_position
      
      # --------------------------------------------------------
      # Write partial checkpoint when buffer reaches threshold
      # --------------------------------------------------------
      
      if (buffer_success_count >= checkpoint_every) {
        
        metadata_chunk <- plyr::rbind.fill(
          metadata_buffer
        )
        
        data.table::fwrite(
          metadata_chunk,
          partial_checkpoint_path,
          append = file.exists(
            partial_checkpoint_path
          ),
          col.names = !file.exists(
            partial_checkpoint_path
          )
        )
        
        message(
          "Partial metadata checkpoint updated: ",
          partial_checkpoint_path,
          " | ",
          success_count,
          " successful accession(s) total"
        )
        
        metadata_buffer <- list()
        buffer_success_count <- 0L
        
        gc()
      }
      
      Sys.sleep(
        get_sleep_duration()
      )
    }
    
    # ==========================================================
    # Flush remaining records
    # ==========================================================
    
    if (length(metadata_buffer) > 0) {
      
      metadata_chunk <- plyr::rbind.fill(
        metadata_buffer
      )
      
      data.table::fwrite(
        metadata_chunk,
        partial_checkpoint_path,
        append = file.exists(
          partial_checkpoint_path
        ),
        col.names = !file.exists(
          partial_checkpoint_path
        )
      )
      
      message(
        "Final partial metadata chunk written: ",
        partial_checkpoint_path
      )
      
      metadata_buffer <- list()
      buffer_success_count <- 0L
      
      gc()
    }
    
    # ==========================================================
    # Finalize search-group checkpoint
    # ==========================================================
    
    if (file.exists(partial_checkpoint_path)) {
      
      completed_metadata <- read.csv(
        partial_checkpoint_path,
        stringsAsFactors = FALSE,
        colClasses = "character"
      )
      
      completed_metadata <- dplyr::distinct(
        completed_metadata,
        Accession,
        .keep_all = TRUE
      )
      
      write.csv(
        completed_metadata,
        checkpoint_path,
        row.names = FALSE
      )
      
      file.remove(
        partial_checkpoint_path
      )
      
      message(
        "Completed metadata checkpoint written: ",
        checkpoint_path
      )
      
    } else if (file.exists(checkpoint_path)) {
      
      message(
        "Completed checkpoint already exists: ",
        checkpoint_path
      )
      
    } else {
      
      warning(
        "No metadata successfully retrieved for search group: ",
        tx
      )
    }
    
    # ==========================================================
    # Failed accession log
    # ==========================================================
    
    if (nrow(failed_log) > 0) {
      
      write.csv(
        failed_log,
        failed_path,
        row.names = FALSE
      )
      
      message(
        "Failed accession log written: ",
        failed_path
      )
    }
    
    # ==========================================================
    # Timing log for this search group
    # ==========================================================
    
    taxon_end <- Sys.time()
    
    taxon_elapsed <- as.numeric(
      difftime(
        taxon_end,
        taxon_start,
        units = "mins"
      )
    )
    
    timing_log <- rbind(
      timing_log,
      data.frame(
        Taxon = tx,
        Num_accessions = nrow(block),
        Num_successful = success_count,
        Num_failed = fail_count,
        Start_time = format(
          taxon_start,
          "%Y-%m-%d %H:%M:%S"
        ),
        End_time = format(
          taxon_end,
          "%Y-%m-%d %H:%M:%S"
        ),
        Elapsed_minutes = round(
          taxon_elapsed,
          2
        ),
        Checkpoint_file = checkpoint_path,
        stringsAsFactors = FALSE
      )
    )
    
    write.csv(
      timing_log,
      timing_path,
      row.names = FALSE
    )
    
    message(
      tx,
      " completed in ",
      round(taxon_elapsed, 2),
      " minutes. Successful: ",
      success_count,
      ". Failed: ",
      fail_count,
      "."
    )
  }
  
  # ============================================================
  # Retry accessions that failed during normal retrieval
  #
  # Pass 1:
  #   retry every accession currently in the failure log.
  #
  # Pass 2:
  #   retry only the accessions that still failed after pass 1.
  #
  # After failed_retry_passes, stop automatically.
  # ============================================================
  
  if (
    failed_retry_passes > 0 &&
    file.exists(failed_path)
  ) {
    
    retry_failed_log <- read.csv(
      failed_path,
      stringsAsFactors = FALSE,
      colClasses = "character"
    )
    
    if (
      "Accession" %in% names(retry_failed_log) &&
      nrow(retry_failed_log) > 0
    ) {
      
      retry_failed_log <- retry_failed_log[
        !is.na(retry_failed_log$Accession) &
          nzchar(retry_failed_log$Accession),
        ,
        drop = FALSE
      ]
      
      retry_failed_log$.normalized_accession <- normalize_accession(
        retry_failed_log$Accession
      )
      
      retry_failed_log <- retry_failed_log[
        !duplicated(retry_failed_log$.normalized_accession),
        ,
        drop = FALSE
      ]
      
      retry_failed_log$.normalized_accession <- NULL
      
      message(
        "\n============================================================"
      )
      message(
        "Beginning final failed-accession retry stage."
      )
      message(
        nrow(retry_failed_log),
        " unique failed accession(s) currently remain."
      )
      message(
        "Maximum additional retry passes: ",
        failed_retry_passes
      )
      message(
        "============================================================"
      )
      
      # --------------------------------------------------------
      # Load metadata that may already have been recovered by
      # a previous retry run
      # --------------------------------------------------------
      
      recovered_accessions <- character(0)
      
      if (file.exists(recovered_failed_path)) {
        
        existing_recovered <- read.csv(
          recovered_failed_path,
          stringsAsFactors = FALSE,
          colClasses = "character"
        )
        
        if (nrow(existing_recovered) > 0) {
          
          recovered_col <- if (
            "GBSeq_accession.version" %in%
            names(existing_recovered)
          ) {
            "GBSeq_accession.version"
          } else {
            "Accession"
          }
          
          recovered_accessions <- normalize_accession(
            existing_recovered[[recovered_col]]
          )
          
          retry_failed_log <- retry_failed_log[
            !normalize_accession(
              retry_failed_log$Accession
            ) %in% recovered_accessions,
            ,
            drop = FALSE
          ]
          
          if (nrow(retry_failed_log) == 0) {
            message(
              "\nAll failed accessions were already recovered ",
              "by an earlier retry run."
            )
          }
        }
      }
      
      # ========================================================
      # Retry passes
      # ========================================================
      
      for (
        retry_pass in seq_len(failed_retry_passes)
      ) {
        
        if (nrow(retry_failed_log) == 0) {
          break
        }
        
        retry_start <- Sys.time()
        
        retry_accessions <- retry_failed_log$Accession
        
        retry_taxon_lookup <- retry_failed_log[, c(
          "Accession",
          "Taxon"
        ), drop = FALSE]
        
        retry_taxon_lookup$.normalized_accession <-
          normalize_accession(
            retry_taxon_lookup$Accession
          )
        
        retry_taxon_lookup <- retry_taxon_lookup[
          !duplicated(
            retry_taxon_lookup$.normalized_accession
          ),
          ,
          drop = FALSE
        ]
        
        retry_batches <- split(
          retry_accessions,
          ceiling(
            seq_along(retry_accessions) /
              metadata_batch_size
          )
        )
        
        num_retry_batches <- length(retry_batches)
        
        message(
          "\n------------------------------------------------------------"
        )
        message(
          "FAILED ACCESSION RETRY PASS ",
          retry_pass,
          " / ",
          failed_retry_passes
        )
        message(
          "Retrying ",
          length(retry_accessions),
          " accession(s) in ",
          num_retry_batches,
          " batch(es)."
        )
        message(
          "------------------------------------------------------------"
        )
        
        remaining_failed_accessions <- character(0)
        
        recovered_buffer <- list()
        recovered_buffer_count <- 0L
        
        retry_success_count <- 0L
        retry_fail_count <- 0L
        
        # ======================================================
        # Retry batch loop
        # ======================================================
        
        for (
          retry_batch_index in seq_along(
            retry_batches
          )
        ) {
          
          accession_batch <- unname(
            retry_batches[[retry_batch_index]]
          )
          
          if (
            retry_batch_index == 1 ||
            retry_batch_index %% progress_every == 0 ||
            retry_batch_index == num_retry_batches
          ) {
            
            message(
              "Retry pass ",
              retry_pass,
              " | batch ",
              retry_batch_index,
              " / ",
              num_retry_batches,
              " | ",
              length(accession_batch),
              " accession(s)"
            )
          }
          
          batch_result <- fetch_metadata_batch_resilient(
            accessions = accession_batch,
            min_batch_size = min_batch_size,
            max_retries = batch_max_retries,
            retry_wait = batch_retry_wait,
            ncbi_database = ncbi_database
          )
          
          # ----------------------------------------------------
          # Identify returned accessions
          # ----------------------------------------------------
          
          returned_accessions <- character(0)
          
          if (
            !is.null(batch_result) &&
            nrow(batch_result) > 0
          ) {
            
            returned_column <- if (
              "GBSeq_accession.version" %in%
              names(batch_result)
            ) {
              "GBSeq_accession.version"
            } else {
              "Accession"
            }
            
            returned_accessions <- normalize_accession(
              batch_result[[returned_column]]
            )
            
            retry_success_count <- retry_success_count +
              nrow(batch_result)
            
            # --------------------------------------------------
            # Add recovered metadata to retry buffer
            # --------------------------------------------------
            
            for (
              row_index in seq_len(
                nrow(batch_result)
              )
            ) {
              
              recovered_buffer[[
                length(recovered_buffer) + 1L
              ]] <- batch_result[
                row_index,
                ,
                drop = FALSE
              ]
              
              recovered_buffer_count <-
                recovered_buffer_count + 1L
            }
          }
          
          # ----------------------------------------------------
          # Determine what is still missing
          # ----------------------------------------------------
          
          missing_accessions <- accession_batch[
            !normalize_accession(accession_batch) %in%
              returned_accessions
          ]
          
          if (length(missing_accessions) > 0) {
            
            retry_fail_count <- retry_fail_count +
              length(missing_accessions)
            
            remaining_failed_accessions <- c(
              remaining_failed_accessions,
              missing_accessions
            )
          }
          
          # ----------------------------------------------------
          # Incrementally save recovered retry metadata
          # ----------------------------------------------------
          
          if (
            recovered_buffer_count >= checkpoint_every
          ) {
            
            recovered_chunk <- plyr::rbind.fill(
              recovered_buffer
            )
            
            data.table::fwrite(
              recovered_chunk,
              recovered_failed_path,
              append = file.exists(
                recovered_failed_path
              ),
              col.names = !file.exists(
                recovered_failed_path
              )
            )
            
            message(
              "Recovered-accession checkpoint updated: ",
              recovered_failed_path
            )
            
            recovered_buffer <- list()
            recovered_buffer_count <- 0L
            
            gc()
          }
          
          Sys.sleep(
            get_sleep_duration()
          )
        }
        
        # ======================================================
        # Flush remaining recovered metadata
        # ======================================================
        
        if (length(recovered_buffer) > 0) {
          
          recovered_chunk <- plyr::rbind.fill(
            recovered_buffer
          )
          
          data.table::fwrite(
            recovered_chunk,
            recovered_failed_path,
            append = file.exists(
              recovered_failed_path
            ),
            col.names = !file.exists(
              recovered_failed_path
            )
          )
          
          message(
            "Final recovered-accession chunk written: ",
            recovered_failed_path
          )
          
          recovered_buffer <- list()
          recovered_buffer_count <- 0L
          
          gc()
        }
        
        # ======================================================
        # Deduplicate recovered metadata checkpoint
        # ======================================================
        
        if (file.exists(recovered_failed_path)) {
          
          recovered_metadata <- read.csv(
            recovered_failed_path,
            stringsAsFactors = FALSE,
            colClasses = "character"
          )
          
          if (nrow(recovered_metadata) > 0) {
            
            recovered_col <- if (
              "GBSeq_accession.version" %in%
              names(recovered_metadata)
            ) {
              "GBSeq_accession.version"
            } else {
              "Accession"
            }
            
            recovered_metadata$.normalized_accession <-
              normalize_accession(
                recovered_metadata[[recovered_col]]
              )
            
            recovered_metadata <- recovered_metadata[
              !duplicated(
                recovered_metadata$.normalized_accession
              ),
              ,
              drop = FALSE
            ]
            
            recovered_metadata$.normalized_accession <- NULL
            
            write.csv(
              recovered_metadata,
              recovered_failed_path,
              row.names = FALSE
            )
          }
        }
        
        # ======================================================
        # Build failure list for next retry pass
        # ======================================================
        
        remaining_failed_accessions <- unique(
          remaining_failed_accessions
        )
        
        if (
          length(remaining_failed_accessions) > 0
        ) {
          
          remaining_normalized <- normalize_accession(
            remaining_failed_accessions
          )
          
          taxon_match <- match(
            remaining_normalized,
            retry_taxon_lookup$.normalized_accession
          )
          
          remaining_taxa <- retry_taxon_lookup$Taxon[
            taxon_match
          ]
          
          remaining_taxa[
            is.na(remaining_taxa)
          ] <- "FAILED_RETRY"
          
          retry_failed_log <- data.frame(
            Taxon = remaining_taxa,
            Accession = remaining_failed_accessions,
            Error = rep(
              paste0(
                "Accession still not returned after failed-accession ",
                "retry pass ",
                retry_pass,
                "."
              ),
              length(remaining_failed_accessions)
            ),
            Time = rep(
              format(
                Sys.time(),
                "%Y-%m-%d %H:%M:%S"
              ),
              length(remaining_failed_accessions)
            ),
            stringsAsFactors = FALSE
          )
          
        } else {
          
          retry_failed_log <- data.frame(
            Taxon = character(),
            Accession = character(),
            Error = character(),
            Time = character(),
            stringsAsFactors = FALSE
          )
        }
        
        # ------------------------------------------------------
        # Rewrite failure file after every retry pass
        #
        # It now contains ONLY accessions that remain unresolved.
        # ------------------------------------------------------
        
        write.csv(
          retry_failed_log,
          failed_path,
          row.names = FALSE
        )
        
        retry_end <- Sys.time()
        
        retry_elapsed <- as.numeric(
          difftime(
            retry_end,
            retry_start,
            units = "mins"
          )
        )
        
        message(
          "\nRetry pass ",
          retry_pass,
          " complete in ",
          round(retry_elapsed, 2),
          " minutes."
        )
        
        message(
          "Recovered this pass: ",
          length(retry_accessions) -
            nrow(retry_failed_log)
        )
        
        message(
          "Still failed: ",
          nrow(retry_failed_log)
        )
      }
      
      # ========================================================
      # Final failed-retry status
      # ========================================================
      
      if (nrow(retry_failed_log) == 0) {
        
        message(
          "\nAll failed metadata accessions were recovered."
        )
        
      } else {
        
        message(
          "\n============================================================"
        )
        
        message(
          "FAILED ACCESSION RETRY LIMIT REACHED"
        )
        
        message(
          nrow(retry_failed_log),
          " accession(s) remain unresolved after ",
          failed_retry_passes,
          " additional retry pass(es)."
        )
        
        message(
          "These accessions will not be retried again automatically ",
          "during this call."
        )
        
        message(
          "Remaining failures are recorded in: ",
          failed_path
        )
        
        message(
          "============================================================"
        )
      }
    }
  }
  
  # ============================================================
  # Combine completed metadata checkpoints
  # ============================================================
  
  checkpoint_files <- list.files(
    checkpoint_dir,
    pattern = "^metadata_.*\\.csv$",
    full.names = TRUE
  )
  
  checkpoint_files <- checkpoint_files[
    !grepl(
      "metadata_failed_accessions\\.csv$",
      checkpoint_files
    )
  ]
  
  if (length(checkpoint_files) == 0) {
    stop(
      "No checkpoint metadata files found in: ",
      checkpoint_dir
    )
  }
  
  message(
    "\nCombining ",
    length(checkpoint_files),
    " checkpoint file(s)."
  )
  
  metadata_database <- plyr::rbind.fill(
    lapply(
      checkpoint_files,
      function(x) {
        read.csv(
          x,
          stringsAsFactors = FALSE
        )
      }
    )
  )
  
  metadata_database <- dplyr::distinct(
    metadata_database,
    Accession,
    .keep_all = TRUE
  )
  
  # ============================================================
  # Write final metadata
  # ============================================================
  
  final_path <- paste0(
    "./metadata_files/all_accessions_pulled_metadata_",
    project_name,
    ".csv"
  )
  
  write.csv(
    metadata_database,
    final_path,
    row.names = FALSE
  )
  
  # ============================================================
  # Final timing summary
  # ============================================================
  
  overall_end <- Sys.time()
  
  total_elapsed <- as.numeric(
    difftime(
      overall_end,
      overall_start,
      units = "mins"
    )
  )
  
  final_failed_count <- 0L
  
  if (file.exists(failed_path)) {
    
    final_failed_file <- read.csv(
      failed_path,
      stringsAsFactors = FALSE
    )
    
    final_failed_count <- nrow(
      final_failed_file
    )
  }
  
  timing_log <- rbind(
    timing_log,
    data.frame(
      Taxon = "TOTAL",
      Num_accessions = nrow(accession_list),
      Num_successful = nrow(metadata_database),
      Num_failed = final_failed_count,
      Start_time = format(
        overall_start,
        "%Y-%m-%d %H:%M:%S"
      ),
      End_time = format(
        overall_end,
        "%Y-%m-%d %H:%M:%S"
      ),
      Elapsed_minutes = round(
        total_elapsed,
        2
      ),
      Checkpoint_file = final_path,
      stringsAsFactors = FALSE
    )
  )
  
  write.csv(
    timing_log,
    timing_path,
    row.names = FALSE
  )
  
  message("\nMetadata retrieval complete.")
  
  message(
    "Final metadata written to: ",
    final_path
  )
  
  message(
    "Timing log written to: ",
    timing_path
  )
  
  if (file.exists(failed_path)) {
    
    message(
      "Remaining failed accession log written to: ",
      failed_path
    )
    
    message(
      "Final unresolved accession count: ",
      final_failed_count
    )
  }
  
  if (file.exists(recovered_failed_path)) {
    
    message(
      "Recovered failed-accession metadata written to: ",
      recovered_failed_path
    )
  }
  
  invisible(metadata_database)
}



# Custom sequences merge
merge_metadata_with_custom_file <- function(project_name,
                                            my_lab_sequences = get0("my_lab_sequences", envir = .GlobalEnv, ifnotfound = ""),
                                            metadata_dir = "./metadata_files") {
  
  metadata_file_path <- file.path(
    metadata_dir,
    paste0("all_accessions_pulled_metadata_", project_name, ".csv")
  )
  
  if (!file.exists(metadata_file_path)) {
    stop("Metadata file not found: ", metadata_file_path)
  }
  
  if (is.null(my_lab_sequences) || length(my_lab_sequences) == 0) {
    message("No custom sequences file provided; skipping custom sequence merge.")
    return(invisible(NULL))
  }
  
  if (!is.character(my_lab_sequences) || length(my_lab_sequences) != 1) {
    stop(
      "my_lab_sequences must be a single file path, not a loaded data frame/table.\n",
      "Use:\n",
      '  my_lab_sequences <- "/Users/scott/Desktop/dactylonectria_extracted_sequences.tsv"'
    )
  }
  
  if (!nzchar(my_lab_sequences)) {
    message("No custom sequences file provided; skipping custom sequence merge.")
    return(invisible(NULL))
  }
  
  if (!file.exists(my_lab_sequences)) {
    stop("Custom sequences file not found: ", my_lab_sequences)
  }
  
  metadata_database <- read.csv(
    metadata_file_path,
    stringsAsFactors = FALSE
  )
  
  custom_sequences <- readr::read_delim(
    my_lab_sequences,
    delim = ifelse(
      grepl("\\.tsv$|\\.txt$", my_lab_sequences, ignore.case = TRUE),
      "\t",
      ","
    ),
    show_col_types = FALSE,
    col_types = readr::cols(.default = readr::col_character())
  ) |>
    as.data.frame()
  
  required_cols <- c(
    "Accession",
    "strain",
    "sequence",
    "organism",
    "gene"
  )
  
  missing_required <- setdiff(required_cols, names(custom_sequences))
  
  if (length(missing_required) > 0) {
    stop(
      "Custom file is missing required column(s): ",
      paste(missing_required, collapse = ", "),
      "\nRequired columns are: ",
      paste(required_cols, collapse = ", "),
      "\nDetected columns are: ",
      paste(names(custom_sequences), collapse = ", ")
    )
  }
  
  recommended_cols <- c(
    "product",
    "accession_title",
    "host",
    "isolation_source",
    "GBSeq_taxonomy",
    "Strain.taxonomy"
  )
  
  for (col in recommended_cols) {
    if (!col %in% names(custom_sequences)) {
      custom_sequences[[col]] <- NA_character_
    }
  }
  
  custom_sequences <- custom_sequences %>%
    dplyr::mutate(
      
      # Create globally unique accession names for lab sequences
      Accession = paste(strain, gene, Accession, sep = "_"),
      
      # Populate these fields if absent so region curation can recognize them
      product = ifelse(
        is.na(product) | product == "",
        gene,
        product
      ),
      
      accession_title = ifelse(
        is.na(accession_title) | accession_title == "",
        gene,
        accession_title
      ),
      
      custom_sequence = TRUE
    )
  
  if (!"custom_sequence" %in% names(metadata_database)) {
    metadata_database$custom_sequence <- FALSE
  }
  
  merged_data <- plyr::rbind.fill(
    metadata_database,
    custom_sequences
  )
  
  # Remove only true duplicate records
  merged_data <- merged_data %>%
    dplyr::distinct(
      strain,
      organism,
      gene,
      Accession,
      .keep_all = TRUE
    )
  
  write.csv(
    merged_data,
    metadata_file_path,
    row.names = FALSE
  )
  
  message("Merged custom sequences into: ", metadata_file_path)
  message("Custom rows added from: ", my_lab_sequences)
  message("Custom accessions renamed as strain_gene_originalAccession")
  
  invisible(merged_data)
}


# helper fucntion for strain taxonomy pull (part of basic curation)
add_strain_taxonomy_columns_to_df <- function(meta, overwrite = TRUE) {
  
  trim_ws <- function(x) {
    x <- as.character(x)
    x <- gsub("^\\s+|\\s+$", "", x)
    x <- clean_taxon_name(x)
    x[x == ""] <- NA_character_
    x
  }
  
  split_tax <- function(x) {
    if (is.na(x) || x == "") return(character(0))
    parts <- unlist(strsplit(x, ";"))
    parts <- trim_ws(parts)
    parts[!is.na(parts)]
  }
  
  extract_rank <- function(parts, pattern) {
    hit <- grep(pattern, parts, ignore.case = TRUE, value = TRUE)
    if (length(hit) == 0) return(NA_character_)
    hit[1]
  }
  
  parse_genus_species_from_organism <- function(org) {
    if (is.na(org) || org == "") {
      return(list(genus = NA_character_, species = NA_character_))
    }
    
    org <- trim_ws(org)
    parts <- unlist(strsplit(org, "\\s+"))
    
    genus <- if (length(parts) >= 1) parts[1] else NA_character_
    species <- NA_character_
    
    bad_species_terms <- c(
      "sp.", "sp", "cf.", "cf", "aff.", "aff",
      "nr.", "nr", "complex", "group"
    )
    
    if (length(parts) >= 2 && !(tolower(parts[2]) %in% bad_species_terms)) {
      species <- parts[2]
    }
    
    list(genus = genus, species = species)
  }
  
  get_tax_string <- function(i) {
    candidates <- c("Strain.taxonomy", "GBSeq_taxonomy", "taxonomy")
    
    for (col in candidates) {
      if (col %in% names(meta)) {
        val <- trim_ws(meta[[col]][i])
        if (!is.na(val) && val != "") return(val)
      }
    }
    
    NA_character_
  }
  
  parsed <- lapply(seq_len(nrow(meta)), function(i) {
    tax_string <- get_tax_string(i)
    parts <- split_tax(tax_string)
    
    org <- if ("organism" %in% names(meta)) meta$organism[i] else NA_character_
    org_parsed <- parse_genus_species_from_organism(org)
    
    phylum <- extract_rank(parts, "mycota$|mycotina$")
    class  <- extract_rank(parts, "mycetes$")
    order  <- extract_rank(parts, "ales$")
    family <- extract_rank(parts, "aceae$")
    
    genus <- org_parsed$genus
    if (is.na(genus) && length(parts) > 0) {
      genus <- tail(parts, 1)
    }
    
    species <- org_parsed$species
    
    data.frame(
      Strain.taxonomy = tax_string,
      Strain.phylum   = phylum,
      Strain.class    = class,
      Strain.order    = order,
      Strain.family   = family,
      Strain.genus    = genus,
      Strain.species  = species,
      stringsAsFactors = FALSE
    )
  })
  
  parsed_df <- do.call(rbind, parsed)
  
  strain_tax_cols <- c(
    "Strain.taxonomy",
    "Strain.phylum",
    "Strain.class",
    "Strain.order",
    "Strain.family",
    "Strain.genus",
    "Strain.species"
  )
  
  for (col in strain_tax_cols) {
    if (col %in% names(parsed_df)) {
      parsed_df[[col]] <- clean_taxon_name(parsed_df[[col]])
    }
  }
  
  new_cols <- names(parsed_df)
  
  if (overwrite) {
    meta <- meta[, setdiff(names(meta), new_cols), drop = FALSE]
  }
  
  cbind(meta, parsed_df)
}


# Metadata curation and region selection
curate_metadata_basic <- function(project_name,
                                  taxa_of_interest = NULL,
                                  add_strain_taxonomy = TRUE,
                                  overwrite_strain_taxonomy = TRUE) {
  
  metadata_file_path <- paste0(
    "./metadata_files/all_accessions_pulled_metadata_",
    project_name,
    ".csv"
  )
  
  if (!file.exists(metadata_file_path)) {
    stop("Metadata file not found at: ", metadata_file_path)
  }
  
  accession_list <- read.csv(
    metadata_file_path,
    header = TRUE,
    stringsAsFactors = FALSE
  )
  
  # Ensure expected columns exist
  if (!"specimen_voucher" %in% names(accession_list)) accession_list$specimen_voucher <- NA_character_
  if (!"strain" %in% names(accession_list))           accession_list$strain           <- NA_character_
  if (!"isolate" %in% names(accession_list))          accession_list$isolate          <- NA_character_
  if (!"type_material" %in% names(accession_list))    accession_list$type_material    <- NA_character_
  if (!"geo_loc_name" %in% names(accession_list))     accession_list$geo_loc_name     <- NA_character_
  if (!"organism" %in% names(accession_list))         accession_list$organism         <- NA_character_
  if (!"Accession" %in% names(accession_list)) {
    stop("Metadata file must contain an 'Accession' column.")
  }
  
  # Optional organism filtering
  if (!is.null(taxa_of_interest) &&
      length(taxa_of_interest) > 0 &&
      any(!is.na(taxa_of_interest) & nzchar(taxa_of_interest))) {
    
    taxa_of_interest <- taxa_of_interest[!is.na(taxa_of_interest) & nzchar(taxa_of_interest)]
    
    pat <- paste0("^(", paste(taxa_of_interest, collapse = "|"), ")\\b")
    accession_list <- accession_list[grepl(pat, accession_list$organism), ]
  }
  
  # Turn empty strings into NA for strain-name source columns
  name_cols <- c("specimen_voucher", "strain", "isolate")
  accession_list[name_cols] <- lapply(accession_list[name_cols], function(x) {
    x <- as.character(x)
    x[x == ""] <- NA_character_
    x
  })
  
  # Choose the best available strain-like identifier
  accession_list <- accession_list %>%
    dplyr::mutate(
      strain.standard = dplyr::coalesce(
        specimen_voucher,
        strain,
        isolate,
        Accession
      )
    )
  
  # Clean strain names for safe FASTA headers / plotting labels
  remove_char_pattern <- "[><\\s:;_\\-\\.()&|#/\\\\,'\"!?\\[\\]{}+=%\\*\\^~@$]"
  
  accession_list$strain.standard <- stringr::str_remove_all(
    accession_list$strain.standard,
    remove_char_pattern
  )
  
  # Add TYPE suffix when type material is present
  accession_list$strain.standard.type <- ifelse(
    !is.na(accession_list$type_material) & accession_list$type_material != "",
    paste0(accession_list$strain.standard, ".TYPE"),
    accession_list$strain.standard
  )
  
  # Dot-separated organism name for FASTA headers
  accession_list$org_name <- gsub("\\s+", "\\.", accession_list$organism)
  
  # Add Strain.taxonomy and parsed Strain.* taxonomy columns
  if (add_strain_taxonomy) {
    accession_list <- add_strain_taxonomy_columns_to_df(
      accession_list,
      overwrite = overwrite_strain_taxonomy
    )
  }
  
  out_path <- paste0(
    "./metadata_files/all_accessions_pulled_metadata_",
    project_name,
    "_curated.csv"
  )
  
  write.csv(
    accession_list,
    out_path,
    row.names = FALSE
  )
  
  cat("Wrote basic curated metadata to:", out_path, "\n")
  
  if (add_strain_taxonomy) {
    cat("Added/updated Strain.taxonomy and parsed Strain.* taxonomy columns.\n")
  }
  
  invisible(accession_list)
}

# for legacy projects where I didn't have strain taxonomy lookup implemented yet
add_strain_taxonomy_columns <- function(metadata_file, overwrite = TRUE) {
  if (!file.exists(metadata_file)) {
    stop("Metadata file not found: ", metadata_file)
  }

  meta <- read.csv(metadata_file, stringsAsFactors = FALSE)

  meta <- add_strain_taxonomy_columns_to_df(
    meta,
    overwrite = overwrite
  )

  write.csv(meta, metadata_file, row.names = FALSE)

  message("Added/updated fungal taxonomy columns in: ", metadata_file)

  invisible(meta)
}


# ============================================================
# Host assessment functions
# ============================================================


initialize_host_standardized <- function(
    project_name,
    metadata_dir = "./metadata_files",
    use_isolation_source = FALSE,
    overwrite = FALSE
) {
  metadata_path <- file.path(
    metadata_dir,
    paste0("all_accessions_pulled_metadata_", project_name, "_curated.csv")
  )
  
  if (!file.exists(metadata_path)) {
    stop("Curated metadata file not found at: ", metadata_path)
  }
  
  df <- read.csv(metadata_path, stringsAsFactors = FALSE)
  
  # If host.standardized already exists and we're not overwriting, just return
  if ("host.standardized" %in% names(df) && !overwrite) {
    message("host.standardized already exists and overwrite = FALSE; leaving as-is.")
    return(invisible(metadata_path))
  }
  
  # Ensure host + isolation_source columns exist
  if (!"host" %in% names(df)) {
    df$host <- NA_character_
  }
  if (!"isolation_source" %in% names(df)) {
    df$isolation_source <- NA_character_
  }
  
  # Start from host
  host_std <- df$host
  
  # Optionally fill NAs/empties with isolation_source
  if (use_isolation_source) {
    host_is_na_or_empty <- is.na(host_std) | host_std == ""
    iso_vals <- df$isolation_source
    iso_vals[iso_vals == ""] <- NA
    host_std[host_is_na_or_empty] <- iso_vals[host_is_na_or_empty]
  }
  
  # Assign into dataframe
  df$host.standardized <- host_std
  
  # Write back
  write.csv(df, metadata_path, row.names = FALSE)
  message("Initialized host.standardized in: ", metadata_path)
  
  invisible(metadata_path)
}


prepare_host_terms <- function(
    project_name,
    metadata_dir = "./metadata_files",
    host_dir     = "./host_assessment",
    use_isolation_source = FALSE
) {
  # ensure host_assessment folder exists
  if (!dir.exists(host_dir)) dir.create(host_dir, recursive = TRUE)
  
  metadata_path <- file.path(
    metadata_dir,
    paste0("all_accessions_pulled_metadata_", project_name, "_curated.csv")
  )
  
  if (!file.exists(metadata_path)) {
    stop("Curated metadata file not found at: ", metadata_path)
  }
  
  # Make sure host.standardized exists; do NOT overwrite if user already curated it
  initialize_host_standardized(
    project_name         = project_name,
    metadata_dir         = metadata_dir,
    use_isolation_source = use_isolation_source,
    overwrite            = FALSE
  )
  
  # Re-read after possible initialization
  df <- read.csv(metadata_path, stringsAsFactors = FALSE)
  
  if (!"host.standardized" %in% names(df)) {
    stop("host.standardized column is missing even after initialization; something went wrong.")
  }
  
  host_terms <- df$host.standardized
  
  # Drop NA, empty strings, and literal "NA" (any capitalization, with/without spaces)
  bad <- is.na(host_terms) |
    host_terms == "" |
    toupper(trimws(host_terms)) == "NA"
  
  host_terms <- unique(host_terms[!bad])
  
  out_path <- file.path(
    host_dir,
    paste0("host_terms_for_taxonomy_", project_name, ".csv")
  )
  
  # Keep column name 'host' for compatibility with current run_host_taxonomy_lookup()
  write.csv(
    data.frame(host = host_terms, stringsAsFactors = FALSE),
    out_path,
    row.names = FALSE
  )
  
  message("Wrote standardized host terms list to: ", out_path)
  return(out_path)
}



# to count the number of accessions per failed term
.add_failed_host_term_counts <- function(failed_df,
                                         metadata_file) {
  
  if (nrow(failed_df) == 0) {
    failed_df$accession_count <- integer(0)
    return(failed_df)
  }
  
  meta <- read.csv(
    metadata_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  if (!"host.standardized" %in% names(meta)) {
    stop(
      "host.standardized column not found in metadata: ",
      metadata_file
    )
  }
  
  accession_col <- intersect(
    c("Accession", "accession"),
    names(meta)
  )
  
  if (length(accession_col) == 0) {
    stop("Could not find accession column in metadata.")
  }
  
  accession_col <- accession_col[1]
  
  # Normalize exactly the same way as the taxonomy search
  standardized_terms <- .host_clean_term(
    meta$host.standardized
  )
  
  failed_terms <- .host_clean_term(
    failed_df$original_term
  )
  
  failed_df$accession_count <- vapply(
    failed_terms,
    function(term) {
      
      if (is.na(term) || term == "") {
        return(0L)
      }
      
      idx <- which(
        !is.na(standardized_terms) &
          standardized_terms == term
      )
      
      if (length(idx) == 0) {
        return(0L)
      }
      
      length(
        unique(
          meta[[accession_col]][idx]
        )
      )
    },
    integer(1)
  )
  
  failed_df
}

run_host_taxonomy_lookup <- function(
    project_name,
    host_dir = "./host_assessment",
    db = "ncbi",
    sleep_sec = 0.1
) {
  if (!dir.exists(host_dir)) dir.create(host_dir, recursive = TRUE)
  
  terms_path <- file.path(
    host_dir,
    paste0("host_terms_for_taxonomy_", project_name, ".csv")
  )
  if (!file.exists(terms_path)) {
    stop("Host terms file not found at: ", terms_path,
         "\nRun prepare_host_terms() first.")
  }
  
  # host_terms_for_taxonomy has a column named 'host'
  host_terms <- read.csv(terms_path, stringsAsFactors = FALSE)$host
  
  # Sanity clean: drop NA / "" / literal "NA"
  bad <- is.na(host_terms) |
    host_terms == "" |
    toupper(trimws(host_terms)) == "NA"
  host_terms <- unique(host_terms[!bad])
  
  # Paths for taxonomy + failed mapping
  taxonomy_path <- file.path(
    host_dir,
    paste0("host_taxonomy_", project_name, ".csv")
  )
  failed_path <- file.path(
    host_dir,
    paste0("host_failed_terms_", project_name, ".csv")
  )
  
  # ---- helper to validate taxize result ----
  is_valid_tax_table <- function(x) {
    is.data.frame(x) &&
      nrow(x) > 0 &&
      all(c("name", "rank") %in% names(x))
  }
  
  # Load existing taxonomy (if any), with legacy compatibility
  if (file.exists(taxonomy_path)) {
    host_taxonomy <- read.csv(taxonomy_path, stringsAsFactors = FALSE)
    
    # Legacy: older pipeline used "Host.standard"
    if (!"Host.standardized" %in% names(host_taxonomy)) {
      if ("Host.standard" %in% names(host_taxonomy)) {
        message("Detected legacy taxonomy file. Renaming 'Host.standard' to 'Host.standardized'.")
        names(host_taxonomy)[names(host_taxonomy) == "Host.standard"] <- "Host.standardized"
      } else {
        stop(
          "Existing host_taxonomy file is missing 'Host.standardized' and 'Host.standard' columns: ",
          taxonomy_path,
          "\nIf this is an old file and you don't care about it, you can delete it and rerun."
        )
      }
    }
    
    already_done <- unique(host_taxonomy$Host.standardized)
  } else {
    host_taxonomy <- data.frame(
      Host.standardized = character(0),
      Host.kingdom      = character(0),
      Host.phylum       = character(0),
      Host.class        = character(0),
      Host.order        = character(0),
      Host.family       = character(0),
      Host.genus        = character(0),
      Host.species      = character(0),
      stringsAsFactors  = FALSE
    )
    already_done <- character(0)
  }
  
  # Terms that still need to be queried
  to_query <- setdiff(host_terms, already_done)
  
  if (length(to_query) == 0) {
    message("All host terms in ", terms_path, " already have taxonomy.")
    # still show unresolved failed terms if any
    if (file.exists(failed_path)) {
      failed_df <- read.csv(failed_path, stringsAsFactors = FALSE)
      
      failed_df <- add_failed_host_term_counts(
        failed_df,
        project_name = project_name
      )
      
      write.csv(failed_df, failed_path, row.names = FALSE)
      unresolved <- subset(failed_df,
                           is.na(replacement_term) | replacement_term == "")
      message("Unresolved failed terms in mapping file: ", nrow(unresolved))
      if (nrow(unresolved) > 0) {
        n_show <- min(10, nrow(unresolved))
        message("Example unresolved terms (showing ", n_show, "):")
        message("  - ", paste(unresolved$original_term[1:n_show],
                              collapse = "\n  - "))
      }
    }
    return(invisible(list(
      taxonomy = host_taxonomy,
      newly_failed = character(0)
    )))
  }
  
  message("Host taxonomy lookup starting for ", length(to_query),
          " new standardized host terms.")
  
  # Containers for this run
  successful_rows <- list()
  failed_terms    <- character(0)
  
  # Ranks we care about
  ranks_of_interest <- c("kingdom", "phylum", "class",
                         "order", "family", "genus", "species")
  
  for (term in to_query) {
    message("  Querying: ", term)
    
    # extra guard against "NA" / empty
    if (is.na(term) || term == "" || toupper(trimws(term)) == "NA") {
      message("    Skipping invalid term: ", term)
      next
    }
    
    res_list <- tryCatch({
      taxize::classification(term, db = db)
    }, error = function(e) {
      NULL
    })
    
    # If classification() itself failed or returned something weird
    if (is.null(res_list) ||
        length(res_list) == 0 ||
        is.atomic(res_list)) {
      message("    FAILED to retrieve taxonomy for: ", term, " (no usable result)")
      failed_terms <- c(failed_terms, term)
      Sys.sleep(sleep_sec)
      next
    }
    
    res <- res_list[[1]]
    
    # If the first element is not a proper tax table, treat as failure
    if (!is_valid_tax_table(res)) {
      message("    FAILED to retrieve taxonomy for: ", term, " (invalid tax table)")
      failed_terms <- c(failed_terms, term)
      Sys.sleep(sleep_sec)
      next
    }
    
    # wide format: one row, columns Host.kingdom, Host.phylum, ...
    this_row <- setNames(
      as.list(rep(NA_character_, length(ranks_of_interest))),
      paste0("Host.", ranks_of_interest)
    )
    
    for (rk in ranks_of_interest) {
      hit <- res$name[res$rank == rk]
      if (length(hit) > 0) {
        this_row[[paste0("Host.", rk)]] <- hit[1]
      }
    }
    
    df_row <- data.frame(
      Host.standardized = term,
      as.data.frame(this_row, stringsAsFactors = FALSE),
      stringsAsFactors = FALSE
    )
    
    successful_rows[[length(successful_rows) + 1L]] <- df_row
    message("    OK")
    
    Sys.sleep(sleep_sec)
  }
  
  # Bind new successes and append to existing taxonomy
  if (length(successful_rows) > 0) {
    new_tax_rows <- do.call(rbind, successful_rows)
    
    host_taxonomy <- dplyr::bind_rows(host_taxonomy, new_tax_rows) %>%
      dplyr::distinct(Host.standardized, .keep_all = TRUE)
  }
  
  # Write updated taxonomy file
  write.csv(host_taxonomy, taxonomy_path, row.names = FALSE)
  message("\nHost taxonomy file written to: ", taxonomy_path)
  
  # Update failed mapping file
  if (file.exists(failed_path)) {
    failed_df <- read.csv(failed_path, stringsAsFactors = FALSE)
  } else {
    failed_df <- data.frame(
      original_term    = character(0),
      replacement_term = character(0),
      stringsAsFactors = FALSE
    )
  }
  
  # Add new failed terms if they are not already present
  for (ft in unique(failed_terms)) {
    if (!ft %in% failed_df$original_term) {
      failed_df <- rbind(
        failed_df,
        data.frame(
          original_term    = ft,
          replacement_term = NA_character_,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  # Add accession counts and sort by importance
  failed_df <- add_failed_host_term_counts(
    failed_df,
    project_name = project_name
  )
  
  # Write failed terms mapping
  write.csv(failed_df, failed_path, row.names = FALSE)
  message("Failed terms mapping file written to: ", failed_path)
  
  # Console summary
  newly_success <- if (length(successful_rows) > 0)
    nrow(do.call(rbind, successful_rows)) else 0
  
  message("\nHost taxonomy lookup completed.")
  message("  Newly successful terms this run: ", newly_success)
  message("  Total successful terms (cumulative): ", nrow(host_taxonomy))
  message("  Newly failed terms this run: ", length(unique(failed_terms)))
  message("  Total failed terms in mapping file: ", nrow(failed_df))
  
  unresolved <- subset(failed_df,
                       is.na(replacement_term) | replacement_term == "")
  if (nrow(unresolved) > 0) {
    n_show <- min(10, nrow(unresolved))
    message("  Unresolved failed terms needing manual curation (showing ",
            n_show, "):")
    message("    - ",
            paste(unresolved$original_term[1:n_show], collapse = "\n    - "))
  } else {
    message("  No unresolved failed terms; all failures have a replacement_term assigned.")
  }
  
  invisible(list(
    taxonomy     = host_taxonomy,
    newly_failed = unique(failed_terms),
    failed_table = failed_df
  ))
}


apply_host_standardization_mapping <- function(
    project_name,
    metadata_dir = "./metadata_files",
    host_dir     = "./host_assessment"
) {
  metadata_path <- file.path(
    metadata_dir,
    paste0("all_accessions_pulled_metadata_", project_name, "_curated.csv")
  )
  if (!file.exists(metadata_path)) {
    stop("Curated metadata file not found at: ", metadata_path)
  }
  
  failed_path <- file.path(
    host_dir,
    paste0("host_failed_terms_", project_name, ".csv")
  )
  if (!file.exists(failed_path)) {
    stop("Failed terms mapping file not found at: ", failed_path,
         "\nRun run_host_taxonomy_lookup() at least once first.")
  }
  
  df_meta   <- read.csv(metadata_path, stringsAsFactors = FALSE)
  failed_df <- read.csv(failed_path,  stringsAsFactors = FALSE)
  
  if (!"host.standardized" %in% names(df_meta)) {
    stop("Metadata is missing 'host.standardized'. ",
         "Run initialize_host_standardized() / prepare_host_terms() first.")
  }
  
  if (!all(c("original_term", "replacement_term") %in% names(failed_df))) {
    stop("host_failed_terms file must have 'original_term' and 'replacement_term' columns.")
  }
  
  n_changed_to_na  <- 0L
  n_changed_to_new <- 0L
  
  for (i in seq_len(nrow(failed_df))) {
    orig <- failed_df$original_term[i]
    repl <- failed_df$replacement_term[i]
    
    idx <- which(df_meta$host.standardized == orig)
    
    if (length(idx) == 0) next
    
    # Interpret replacement_term:
    # - NA / "" / "NA" (any capitalization) => drop (set to NA)
    if (is.na(repl) || repl == "" || toupper(trimws(repl)) == "NA") {
      df_meta$host.standardized[idx] <- NA_character_
      n_changed_to_na <- n_changed_to_na + length(idx)
    } else {
      df_meta$host.standardized[idx] <- repl
      n_changed_to_new <- n_changed_to_new + length(idx)
    }
  }
  
  write.csv(df_meta, metadata_path, row.names = FALSE)
  
  message("Applied host standardization mapping to metadata:")
  message("  Rows set to NA (ignored hosts): ", n_changed_to_na)
  message("  Rows set to new standardized terms: ", n_changed_to_new)
  message("Updated metadata written to: ", metadata_path)
  
  invisible(list(
    changed_to_na  = n_changed_to_na,
    changed_to_new = n_changed_to_new
  ))
}


merge_host_taxonomy_into_metadata <- function(
    project_name,
    host_dir = "./host_assessment"
) {
  # 1) Load curated metadata
  meta_path <- paste0(
    "./metadata_files/all_accessions_pulled_metadata_",
    project_name,
    "_curated.csv"
  )
  if (!file.exists(meta_path)) {
    stop("Curated metadata file not found at: ", meta_path,
         "\nRun curate_metadata_basic() first.")
  }
  
  meta <- read.csv(meta_path, stringsAsFactors = FALSE)
  
  # 2) Ensure host.standardized exists (default = original 'host')
  if (!"host.standardized" %in% names(meta)) {
    if ("host" %in% names(meta)) {
      meta$host.standardized <- meta$host
    } else {
      stop("Metadata does not contain a 'host' column to initialize 'host.standardized'.")
    }
  }
  
  # 3) Load host taxonomy table
  host_tax_path <- file.path(host_dir, paste0("host_taxonomy_", project_name, ".csv"))
  if (!file.exists(host_tax_path)) {
    stop("Host taxonomy file not found at: ", host_tax_path,
         "\nRun run_host_taxonomy_lookup(project_name) first.")
  }
  
  host_tax <- read.csv(host_tax_path, stringsAsFactors = FALSE)
  
  if (!"Host.standardized" %in% names(host_tax)) {
    stop(
      "Host taxonomy file does not contain 'Host.standardized' column: ",
      host_tax_path
    )
  }
  
  # 3b) Deduplicate host_taxonomy by Host.standardized (keep first)
  host_tax <- host_tax %>%
    dplyr::filter(!is.na(Host.standardized) & Host.standardized != "") %>%
    dplyr::distinct(Host.standardized, .keep_all = TRUE)
  
  # 4) DROP any existing Host.* columns from metadata to avoid .x/.y/.x.x zoo
  existing_host_cols <- grep("^Host\\.", names(meta), value = TRUE)
  if (length(existing_host_cols) > 0) {
    message("Removing existing Host.* columns from metadata before merge: ",
            paste(existing_host_cols, collapse = ", "))
    meta <- meta[, setdiff(names(meta), existing_host_cols), drop = FALSE]
  }
  
  # 5) Left-join taxonomy on host.standardized
  meta_merged <- meta %>%
    dplyr::left_join(
      host_tax,
      by = c("host.standardized" = "Host.standardized")
    )
  
  # 6) Ensure Strain.taxonomy exists, populated from GBSeq_taxonomy if available
  if (!"Strain.taxonomy" %in% names(meta_merged)) {
    if ("GBSeq_taxonomy" %in% names(meta_merged)) {
      meta_merged$Strain.taxonomy <- meta_merged$GBSeq_taxonomy
      message("Created 'Strain.taxonomy' column from 'GBSeq_taxonomy'.")
    } else {
      meta_merged$Strain.taxonomy <- NA_character_
      warning("Neither 'Strain.taxonomy' nor 'GBSeq_taxonomy' found; ",
              "created empty 'Strain.taxonomy' column.")
    }
  }
  
  # 7) Write updated metadata (overwriting previous curated file)
  write.csv(
    meta_merged,
    meta_path,
    row.names = FALSE
  )
  
  message("Merged host taxonomy into metadata and wrote updated file:\n  ", meta_path)
  
  # 8) Summary
  n_total <- nrow(meta_merged)
  
  n_host_raw <- if ("host" %in% names(meta_merged)) {
    sum(!is.na(meta_merged$host) & meta_merged$host != "")
  } else 0L
  
  n_host_std <- sum(!is.na(meta_merged$host.standardized) & meta_merged$host.standardized != "")
  
  host_phylum_col <- "Host.phylum"
  if (host_phylum_col %in% names(meta_merged)) {
    n_host_phylum <- sum(!is.na(meta_merged[[host_phylum_col]]) &
                           meta_merged[[host_phylum_col]] != "")
    n_both_std_phylum <- sum(
      !is.na(meta_merged$host.standardized) & meta_merged$host.standardized != "" &
        !is.na(meta_merged[[host_phylum_col]]) & meta_merged[[host_phylum_col]] != ""
    )
  } else {
    n_host_phylum <- NA_integer_
    n_both_std_phylum <- NA_integer_
  }
  
  message("Host taxonomy summary:")
  message("  Total accessions: ", n_total)
  message("  Accessions with non-empty raw 'host': ", n_host_raw)
  message("  Accessions with non-empty 'host.standardized': ", n_host_std)
  if (!is.na(n_host_phylum)) {
    message("  Accessions with annotated Host.phylum: ", n_host_phylum)
    message("  Accessions with BOTH non-empty host.standardized and Host.phylum: ",
            n_both_std_phylum)
  }
  
  invisible(meta_merged)
}


summarize_host_usage <- function(
    project_name,
    fungal_rank = "species",
    host_rank   = "phylum",
    keep_NAs    = FALSE,
    host_dir    = "./host_assessment",
    metadata_file = NULL
) {
  if (!dir.exists(host_dir)) dir.create(host_dir, recursive = TRUE)
  
  if (is.null(metadata_file)) {
    metadata_file <- paste0(
      "./metadata_files/all_accessions_pulled_metadata_",
      project_name,
      "_curated.csv"
    )
  }
  
  if (!file.exists(metadata_file)) {
    stop("Metadata file not found at: ", metadata_file)
  }
  
  meta <- read.csv(metadata_file, stringsAsFactors = FALSE)
  
  if (fungal_rank == "species") {
    fungal_col <- "org_name"
    fungal_output_col <- "Strain.species"
    
    if (!fungal_col %in% names(meta)) {
      stop(
        "fungal_rank = 'species' requires column 'org_name' in metadata.\n",
        "Available columns: ",
        paste(names(meta), collapse = ", ")
      )
    }
  } else {
    fungal_col <- paste0("Strain.", fungal_rank)
    fungal_output_col <- fungal_col
    
    if (!fungal_col %in% names(meta)) {
      stop(
        "Column '", fungal_col, "' not found in metadata.\n",
        "Available Strain.* columns: ",
        paste(grep("^Strain\\.", names(meta), value = TRUE), collapse = ", ")
      )
    }
  }
  
  host_col <- paste0("Host.", host_rank)
  
  if (!host_col %in% names(meta)) {
    stop(
      "Column '", host_col, "' not found in metadata.\n",
      "Available Host.* columns: ",
      paste(grep("^Host\\.", names(meta), value = TRUE), collapse = ", ")
    )
  }
  
  df <- meta %>%
    dplyr::select(
      fungal_value = dplyr::all_of(fungal_col),
      host_value   = dplyr::all_of(host_col)
    ) %>%
    dplyr::mutate(
      fungal_value = as.character(fungal_value),
      host_value   = as.character(host_value)
    )
  
  if (keep_NAs) {
    df <- df %>%
      dplyr::mutate(
        host_value = dplyr::if_else(
          is.na(host_value) | host_value == "",
          "NoData",
          host_value
        ),
        fungal_value = dplyr::if_else(
          is.na(fungal_value) | fungal_value == "",
          "NoData",
          fungal_value
        )
      )
  } else {
    df <- df %>%
      dplyr::filter(
        !is.na(host_value) & host_value != "",
        !is.na(fungal_value) & fungal_value != ""
      )
  }
  
  if (nrow(df) == 0) {
    warning("No rows remain after filtering.")
    return(invisible(NULL))
  }
  
  counts_long <- df %>%
    dplyr::count(fungal_value, host_value, name = "accession_count")
  
  counts_wide <- counts_long %>%
    tidyr::pivot_wider(
      names_from = host_value,
      values_from = accession_count,
      values_fill = list(accession_count = 0)
    )
  
  names(counts_wide)[1] <- fungal_output_col
  host_cols <- setdiff(names(counts_wide), fungal_output_col)
  
  counts_wide <- counts_wide %>%
    dplyr::mutate(
      accession.count = rowSums(dplyr::across(dplyr::all_of(host_cols)))
    )
  
  result <- counts_wide %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(host_cols),
        ~ ifelse(accession.count > 0, .x / accession.count * 100, 0),
        .names = "{.col}_percentage"
      )
    )
  
  perc_cols <- paste0(host_cols, "_percentage")
  host_cols_local <- host_cols
  perc_cols_local <- perc_cols
  
  result <- result %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      top_host_categories = {
        vals <- c_across(dplyr::all_of(host_cols_local))
        if (all(is.na(vals)) || all(vals == 0, na.rm = TRUE)) {
          NA_character_
        } else {
          max_val <- max(vals, na.rm = TRUE)
          paste0(host_cols_local[!is.na(vals) & vals == max_val], collapse = ";")
        }
      },
      host_profile = {
        percs <- c_across(dplyr::all_of(perc_cols_local))
        names(percs) <- host_cols_local
        nz <- percs[!is.na(percs) & percs > 0]
        
        if (length(nz) == 0) {
          NA_character_
        } else {
          ord <- order(nz, decreasing = TRUE)
          paste(
            paste0(names(nz)[ord], "(", round(nz[ord], 1), "%)"),
            collapse = ";"
          )
        }
      }
    ) %>%
    dplyr::ungroup()
  
  out_path <- file.path(
    host_dir,
    paste0(
      "host_usage_fungal_",
      fungal_rank,
      "_by_host_",
      host_rank,
      "_",
      project_name,
      ".csv"
    )
  )
  
  write.csv(result, out_path, row.names = FALSE)
  
  message("Host usage summary written to: ", out_path)
  message("Metadata source: ", metadata_file)
  message("Rows in metadata: ", nrow(meta))
  message("Rows used after filtering: ", nrow(df))
  
  invisible(result)
}



##############################
# Host assessment wrappers
##############################

# first pass only
# -----------------------------
# Host assessment helpers
# -----------------------------

.host_metadata_path <- function(project_name) {
  file.path(
    "metadata_files",
    paste0("all_accessions_pulled_metadata_", project_name, "_curated.csv")
  )
}

.host_taxonomy_path <- function(project_name, host_dir = "./host_assessment") {
  file.path(host_dir, paste0("host_taxonomy_", project_name, ".csv"))
}

.host_failed_path <- function(project_name, host_dir = "./host_assessment") {
  file.path(host_dir, paste0("host_failed_terms_", project_name, ".csv"))
}

.host_is_blank <- function(x) {
  is.na(x) | trimws(as.character(x)) == ""
}

.host_is_invalid_tax_term <- function(x) {
  .host_is_blank(x) | toupper(trimws(as.character(x))) == "NA"
}

.host_clean_term <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)] <- ""
  x
}

.host_display_term <- function(x) {
  x <- .host_clean_term(x)
  ifelse(x == "" | toupper(x) == "NA", "NA", x)
}

.host_taxonomy_cols <- function() {
  c(
    "Host.standardized",
    "Host.taxid",
    "Host.matched_name",
    "Host.matched_name_class",
    "Host.lookup_source",
    "Host.superkingdom",
    "Host.kingdom",
    "Host.phylum",
    "Host.class",
    "Host.order",
    "Host.family",
    "Host.genus",
    "Host.species"
  )
}

.empty_host_taxonomy <- function() {
  cols <- .host_taxonomy_cols()
  
  as.data.frame(
    setNames(
      replicate(length(cols), character(0), simplify = FALSE),
      cols
    ),
    stringsAsFactors = FALSE
  )
}

.is_valid_tax_table <- function(x) {
  is.data.frame(x) &&
    nrow(x) > 0 &&
    all(c("name", "rank") %in% names(x))
}

.ensure_failed_columns <- function(failed_df) {
  required_cols <- c(
    "original_term",
    "replacement_term",
    "accession_count",
    "term_type",
    "parent_original_term",
    "lookup_status",
    "replacement_lookup_status",
    "notes"
  )
  
  for (col in required_cols) {
    if (!col %in% names(failed_df)) {
      failed_df[[col]] <- rep(NA_character_, nrow(failed_df))
    }
  }
  
  failed_df <- failed_df[, required_cols, drop = FALSE]
  
  failed_df$original_term <- .host_display_term(failed_df$original_term)
  failed_df$replacement_term <- .host_clean_term(failed_df$replacement_term)
  failed_df$parent_original_term <- .host_clean_term(failed_df$parent_original_term)
  
  failed_df
}

.read_failed_terms <- function(project_name, host_dir = "./host_assessment") {
  failed_path <- .host_failed_path(project_name, host_dir)
  
  if (file.exists(failed_path)) {
    failed_df <- read.csv(failed_path, stringsAsFactors = FALSE)
  } else {
    failed_df <- data.frame(
      original_term = character(0),
      replacement_term = character(0),
      stringsAsFactors = FALSE
    )
  }
  
  .ensure_failed_columns(failed_df)
}

.write_failed_terms <- function(failed_df, project_name, host_dir = "./host_assessment") {
  failed_path <- .host_failed_path(project_name, host_dir)
  failed_df <- .ensure_failed_columns(failed_df)
  
  write.csv(failed_df, failed_path, row.names = FALSE)
  message("Failed terms file written to: ", failed_path)
  
  invisible(failed_df)
}

.read_host_taxonomy <- function(
    project_name,
    host_dir = "./host_assessment"
) {
  taxonomy_path <- .host_taxonomy_path(
    project_name,
    host_dir
  )
  
  taxonomy_cols <- .host_taxonomy_cols()
  
  if (file.exists(taxonomy_path)) {
    
    host_taxonomy <- read.csv(
      taxonomy_path,
      stringsAsFactors = FALSE
    )
    
    if (!"Host.standardized" %in% names(host_taxonomy)) {
      if ("Host.standard" %in% names(host_taxonomy)) {
        
        message(
          "Detected legacy taxonomy file. Renaming 'Host.standard' to 'Host.standardized'."
        )
        
        names(host_taxonomy)[
          names(host_taxonomy) == "Host.standard"
        ] <- "Host.standardized"
        
      } else {
        stop(
          "Existing host taxonomy file is missing 'Host.standardized': ",
          taxonomy_path
        )
      }
    }
    
    for (col in taxonomy_cols) {
      if (!col %in% names(host_taxonomy)) {
        host_taxonomy[[col]] <- NA_character_
      }
    }
    
    host_taxonomy <- host_taxonomy[
      ,
      taxonomy_cols,
      drop = FALSE
    ]
    
    # Keep TaxID type consistent with new taxonomy lookup rows
    host_taxonomy$Host.taxid <- as.character(
      host_taxonomy$Host.taxid
    )
    
  } else {
    
    host_taxonomy <- .empty_host_taxonomy()
  }
  
  host_taxonomy$Host.standardized <-
    .host_clean_term(
      host_taxonomy$Host.standardized
    )
  
  host_taxonomy
}

.write_host_taxonomy <- function(host_taxonomy, project_name, host_dir = "./host_assessment") {
  taxonomy_path <- .host_taxonomy_path(project_name, host_dir)
  taxonomy_cols <- .host_taxonomy_cols()
  
  for (col in taxonomy_cols) {
    if (!col %in% names(host_taxonomy)) {
      host_taxonomy[[col]] <- NA_character_
    }
  }
  
  host_taxonomy <- host_taxonomy[, taxonomy_cols, drop = FALSE]
  host_taxonomy <- host_taxonomy[!duplicated(host_taxonomy$Host.standardized), ]
  
  write.csv(host_taxonomy, taxonomy_path, row.names = FALSE)
  message("Host taxonomy file written to: ", taxonomy_path)
  
  invisible(host_taxonomy)
}

.ensure_host_standardized <- function(
    meta,
    use_isolation_source = FALSE,
    overwrite_host_standardized = FALSE
) {
  
  host_cols <- c(
    "host",
    "Host",
    "host.name",
    "host_name",
    "host.scientific_name",
    "host.scientific.name",
    "host.standard",
    "Host.standard"
  )
  
  isolation_cols <- c(
    "isolation_source",
    "Isolation Source",
    "isolation.source",
    "isolation-source",
    "source_material"
  )
  
  host_cols <- host_cols[host_cols %in% names(meta)]
  isolation_cols <- isolation_cols[isolation_cols %in% names(meta)]
  
  if (!"host.standardized" %in% names(meta) || overwrite_host_standardized) {
    
    if (length(host_cols) == 0 && (!use_isolation_source || length(isolation_cols) == 0)) {
      stop(
        "Could not create 'host.standardized'. No usable host column found",
        if (use_isolation_source) " and no usable isolation source column found." else ".",
        "\nAvailable columns are:\n",
        paste(names(meta), collapse = ", ")
      )
    }
    
    meta$host.standardized <- NA_character_
    
    if (length(host_cols) > 0) {
      host_source <- host_cols[1]
      message("Creating host.standardized from host column: ", host_source)
      meta$host.standardized <- meta[[host_source]]
    }
    
    if (use_isolation_source && length(isolation_cols) > 0) {
      isolation_source <- isolation_cols[1]
      message(
        "Using isolation source column as fallback where host information is missing: ",
        isolation_source
      )
      
      missing_host <- .host_is_invalid_tax_term(meta$host.standardized)
      
      meta$host.standardized[missing_host] <- meta[[isolation_source]][missing_host]
    }
    
  } else {
    message("Using existing host.standardized column.")
  }
  
  meta$host.standardized <- .host_clean_term(meta$host.standardized)
  
  if (!"host.standardized.original" %in% names(meta) || overwrite_host_standardized) {
    meta$host.standardized.original <- meta$host.standardized
  }
  
  meta
}

.count_host_terms <- function(meta) {
  meta <- .ensure_host_standardized(meta)
  
  terms <- .host_display_term(meta$host.standardized)
  
  as.data.frame(table(terms), stringsAsFactors = FALSE) |>
    stats::setNames(c("original_term", "accession_count"))
}

.add_failed_host_term_counts <- function(failed_df, project_name, metadata_file = NULL) {
  if (is.null(metadata_file)) {
    metadata_file <- .host_metadata_path(project_name)
  }
  
  if (!file.exists(metadata_file)) {
    warning("Metadata file not found while adding failed-term counts: ", metadata_file)
    return(failed_df)
  }
  
  meta <- read.csv(metadata_file, stringsAsFactors = FALSE)
  counts <- .count_host_terms(meta)
  
  failed_df <- .ensure_failed_columns(failed_df)
  
  failed_df$accession_count <- counts$accession_count[
    match(failed_df$original_term, counts$original_term)
  ]
  
  failed_df$accession_count[is.na(failed_df$accession_count)] <- 0L
  
  failed_df
}

.lookup_one_host_taxonomy_taxize <- function(
    term,
    db = "ncbi"
) {
  ranks_of_interest <- c(
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species"
  )
  
  res_list <- tryCatch(
    {
      taxize::classification(
        term,
        db = db
      )
    },
    error = function(e) {
      NULL
    }
  )
  
  if (
    is.null(res_list) ||
    length(res_list) == 0 ||
    is.atomic(res_list)
  ) {
    return(NULL)
  }
  
  res <- res_list[[1]]
  
  if (!.is_valid_tax_table(res)) {
    return(NULL)
  }
  
  this_row <- setNames(
    as.list(
      rep(
        NA_character_,
        length(ranks_of_interest)
      )
    ),
    paste0(
      "Host.",
      ranks_of_interest
    )
  )
  
  for (rk in ranks_of_interest) {
    hit <- res$name[
      res$rank == rk
    ]
    
    if (length(hit) > 0) {
      this_row[[paste0("Host.", rk)]] <- hit[1]
    }
  }
  
  df_row <- data.frame(
    Host.standardized = term,
    Host.taxid = NA_character_,
    Host.matched_name = term,
    Host.matched_name_class = NA_character_,
    Host.lookup_source = "taxize_fallback",
    Host.superkingdom = NA_character_,
    as.data.frame(
      this_row,
      stringsAsFactors = FALSE
    ),
    stringsAsFactors = FALSE
  )
  
  df_row[
    ,
    .host_taxonomy_cols(),
    drop = FALSE
  ]
}

.append_failed_terms <- function(
    failed_df,
    terms,
    term_type = "initial_lookup",
    parent_original_term = NA_character_,
    lookup_status = "failed"
) {
  failed_df <- .ensure_failed_columns(failed_df)
  
  terms <- unique(.host_display_term(terms))
  terms <- terms[!is.na(terms) & terms != ""]
  
  for (term in terms) {
    already_present <- term %in% failed_df$original_term
    
    if (!already_present) {
      failed_df <- rbind(
        failed_df,
        data.frame(
          original_term = term,
          replacement_term = NA_character_,
          accession_count = NA_character_,
          term_type = term_type,
          parent_original_term = parent_original_term,
          lookup_status = lookup_status,
          replacement_lookup_status = NA_character_,
          notes = NA_character_,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  .ensure_failed_columns(failed_df)
}


.search_host_terms <- function(
    terms,
    project_name,
    host_dir = "./host_assessment",
    db = "ncbi",
    overwrite = FALSE,
    term_type = "initial_lookup",
    parent_map = NULL,
    use_taxize_fallback = TRUE,
    name_lookup_file = NULL,
    ranked_taxonomy_file = NULL
) {
  if (!dir.exists(host_dir)) {
    dir.create(
      host_dir,
      recursive = TRUE
    )
  }
  
  # ----------------------------------------------------------
  # Load project outputs
  # ----------------------------------------------------------
  
  host_taxonomy <- .read_host_taxonomy(
    project_name,
    host_dir
  )
  
  failed_df <- .read_failed_terms(
    project_name,
    host_dir
  )
  
  ambiguous_df <- .read_ambiguous_terms(
    project_name,
    host_dir
  )
  
  # ----------------------------------------------------------
  # Load local NCBI databases ONCE
  # ----------------------------------------------------------
  
  message("Loading local NCBI taxonomy databases...")
  
  name_lookup <- .read_ncbi_name_lookup(
    name_lookup_file
  )
  
  ranked_taxonomy <- .read_ncbi_ranked_taxonomy(
    ranked_taxonomy_file
  )
  
  message(
    "Local NCBI name records: ",
    format(
      nrow(name_lookup),
      big.mark = ","
    )
  )
  
  message(
    "Local NCBI taxonomy records: ",
    format(
      nrow(ranked_taxonomy),
      big.mark = ","
    )
  )
  
  # ----------------------------------------------------------
  # Clean terms
  # ----------------------------------------------------------
  
  terms <- unique(
    .host_clean_term(terms)
  )
  
  terms <- terms[
    !is.na(terms) &
      nzchar(trimws(terms))
  ]
  
  invalid_terms <- terms[
    .host_is_invalid_tax_term(terms)
  ]
  
  valid_terms <- terms[
    !.host_is_invalid_tax_term(terms)
  ]
  
  if (length(invalid_terms) > 0) {
    failed_df <- .append_failed_terms(
      failed_df,
      invalid_terms,
      term_type = term_type,
      lookup_status = "not_searched_invalid"
    )
  }
  
  already_done <- unique(
    .host_clean_term(
      host_taxonomy$Host.standardized
    )
  )
  
  if (!overwrite) {
    valid_terms <- setdiff(
      valid_terms,
      already_done
    )
  }
  
  # ----------------------------------------------------------
  # Nothing to do
  # ----------------------------------------------------------
  
  if (length(valid_terms) == 0) {
    message(
      "No new valid host terms to query."
    )
    
    failed_df <- .add_failed_host_term_counts(
      failed_df,
      project_name
    )
    
    .write_failed_terms(
      failed_df,
      project_name,
      host_dir
    )
    
    .write_ambiguous_terms(
      ambiguous_df,
      project_name,
      host_dir
    )
    
    return(
      invisible(
        list(
          taxonomy = host_taxonomy,
          failed_table = failed_df,
          ambiguous_table = ambiguous_df,
          newly_failed = character(0),
          newly_ambiguous = character(0),
          newly_successful = character(0)
        )
      )
    )
  }
  
  message(
    "Host taxonomy lookup starting for ",
    length(valid_terms),
    " term(s)."
  )
  
  # ----------------------------------------------------------
  # Containers
  # ----------------------------------------------------------
  
  successful_rows <- list()
  failed_terms <- character(0)
  ambiguous_terms <- character(0)
  
  # ----------------------------------------------------------
  # Lookup
  # ----------------------------------------------------------
  
  for (term in valid_terms) {
    
    message("  Looking up: ", term)
    
    local_match <- .lookup_local_taxid(
      term = term,
      name_lookup = name_lookup
    )
    
    # --------------------------------------------------------
    # UNIQUE LOCAL MATCH
    # --------------------------------------------------------
    
    if (local_match$status == "unique") {
      
      df_row <- .lookup_local_ranked_taxonomy(
        term = term,
        taxid = local_match$taxid,
        matched_name = local_match$matched_name,
        matched_name_class = local_match$matched_name_class,
        ranked_taxonomy = ranked_taxonomy
      )
      
      if (!is.null(df_row)) {
        
        message(
          "    LOCAL OK | TaxID ",
          local_match$taxid
        )
        
        successful_rows[[length(successful_rows) + 1L]] <- df_row
        
        next
      }
      
      # Extremely unusual case:
      # name exists but TaxID absent from ranked taxonomy.
      message(
        "    Local TaxID found but ranked taxonomy was missing."
      )
    }
    
    # --------------------------------------------------------
    # AMBIGUOUS LOCAL MATCH
    # --------------------------------------------------------
    
    if (local_match$status == "ambiguous") {
      
      message(
        "    AMBIGUOUS | ",
        length(
          unique(local_match$hits$TaxID)
        ),
        " TaxIDs"
      )
      
      ambiguous_terms <- c(
        ambiguous_terms,
        term
      )
      
      ambiguous_df <- .append_ambiguous_lookup(
        ambiguity_df = ambiguous_df,
        term = term,
        name_hits = local_match$hits,
        ranked_taxonomy = ranked_taxonomy
      )
      
      # IMPORTANT:
      # never use taxize to guess an ambiguous local name.
      next
    }
    
    # --------------------------------------------------------
    # NO LOCAL MATCH
    # --------------------------------------------------------
    
    if (local_match$status == "not_found") {
      
      message(
        "    No local NCBI name match."
      )
      
    }
    
    # --------------------------------------------------------
    # OPTIONAL TAXIZE FALLBACK
    # --------------------------------------------------------
    
    if (isTRUE(use_taxize_fallback)) {
      
      message(
        "    Trying taxize fallback..."
      )
      
      fallback_row <- .lookup_one_host_taxonomy_taxize(
        term,
        db = db
      )
      
      if (!is.null(fallback_row)) {
        
        message(
          "    TAXIZE OK"
        )
        
        successful_rows[[length(successful_rows) + 1L]] <- fallback_row
        next
      }
    }
    
    # --------------------------------------------------------
    # COMPLETE FAILURE
    # --------------------------------------------------------
    
    message(
      "    FAILED"
    )
    
    failed_terms <- c(
      failed_terms,
      term
    )
  }
  
  # ----------------------------------------------------------
  # Add successes
  # ----------------------------------------------------------
  
  if (length(successful_rows) > 0) {
    
    new_tax_rows <- dplyr::bind_rows(
      successful_rows
    )
    
    if (overwrite) {
      
      host_taxonomy <- host_taxonomy[
        !host_taxonomy$Host.standardized %in%
          new_tax_rows$Host.standardized,
        ,
        drop = FALSE
      ]
    }
    
    host_taxonomy <- dplyr::bind_rows(
      host_taxonomy,
      new_tax_rows
    )
    
    host_taxonomy <- host_taxonomy[
      !duplicated(
        host_taxonomy$Host.standardized
      ),
      ,
      drop = FALSE
    ]
  }
  
  # ----------------------------------------------------------
  # Add failures
  # ----------------------------------------------------------
  
  if (length(failed_terms) > 0) {
    
    failed_df <- .append_failed_terms(
      failed_df,
      failed_terms,
      term_type = term_type,
      lookup_status = "failed"
    )
  }
  
  # ----------------------------------------------------------
  # Counts
  # ----------------------------------------------------------
  
  failed_df <- .add_failed_host_term_counts(
    failed_df,
    project_name
  )
  
  # ----------------------------------------------------------
  # Write everything
  # ----------------------------------------------------------
  
  .write_host_taxonomy(
    host_taxonomy,
    project_name,
    host_dir
  )
  
  .write_failed_terms(
    failed_df,
    project_name,
    host_dir
  )
  
  .write_ambiguous_terms(
    ambiguous_df,
    project_name,
    host_dir
  )
  
  # ----------------------------------------------------------
  # Summary
  # ----------------------------------------------------------
  
  newly_successful <- if (
    length(successful_rows) > 0
  ) {
    unique(
      dplyr::bind_rows(
        successful_rows
      )$Host.standardized
    )
  } else {
    character(0)
  }
  
  message("")
  message("Host taxonomy lookup complete.")
  message(
    "  Successful: ",
    length(newly_successful)
  )
  message(
    "  Ambiguous: ",
    length(unique(ambiguous_terms))
  )
  message(
    "  Failed: ",
    length(unique(failed_terms))
  )
  
  invisible(
    list(
      taxonomy = host_taxonomy,
      failed_table = failed_df,
      ambiguous_table = ambiguous_df,
      newly_failed = unique(failed_terms),
      newly_ambiguous = unique(ambiguous_terms),
      newly_successful = newly_successful
    )
  )
}



.host_example_data_dir <- function() {
  arborist_root <- get0(
    "arborist_repo",
    envir = .GlobalEnv,
    ifnotfound = normalizePath(
      "~/github/aRborist",
      mustWork = FALSE
    )
  )
  
  file.path(
    arborist_root,
    "example_data"
  )
}


.host_name_lookup_path <- function(
    name_lookup_file = NULL
) {
  if (!is.null(name_lookup_file) &&
      length(name_lookup_file) == 1 &&
      nzchar(name_lookup_file)) {
    return(name_lookup_file)
  }
  
  file.path(
    .host_example_data_dir(),
    "ncbi_name_lookup.rds"
  )
}


.host_ranked_taxonomy_path <- function(
    ranked_taxonomy_file = NULL
) {
  if (!is.null(ranked_taxonomy_file) &&
      length(ranked_taxonomy_file) == 1 &&
      nzchar(ranked_taxonomy_file)) {
    return(ranked_taxonomy_file)
  }
  
  file.path(
    .host_example_data_dir(),
    "ncbi_ranked_taxonomy.rds"
  )
}


.host_shared_replacement_path <- function(
    shared_replacement_file = NULL
) {
  if (!is.null(shared_replacement_file) &&
      length(shared_replacement_file) == 1 &&
      nzchar(shared_replacement_file)) {
    return(shared_replacement_file)
  }
  
  file.path(
    .host_example_data_dir(),
    "host_term_replacements.csv"
  )
}


.host_ambiguous_path <- function(
    project_name,
    host_dir = "./host_assessment"
) {
  file.path(
    host_dir,
    paste0(
      "host_ambiguous_terms_",
      project_name,
      ".csv"
    )
  )
}

.host_normalize_name <- function(x) {
  x <- as.character(x)
  
  x <- trimws(x)
  
  # Collapse repeated internal whitespace
  x <- gsub(
    "[[:space:]]+",
    " ",
    x
  )
  
  tolower(x)
}


.read_ncbi_name_lookup <- function(
    name_lookup_file = NULL
) {
  path <- .host_name_lookup_path(
    name_lookup_file
  )
  
  if (!file.exists(path)) {
    stop(
      "NCBI name lookup database not found:\n  ",
      path
    )
  }
  
  x <- readRDS(path)
  
  required <- c(
    "TaxID",
    "name",
    "name_class"
  )
  
  missing <- setdiff(
    required,
    names(x)
  )
  
  if (length(missing) > 0) {
    stop(
      "ncbi_name_lookup.rds is missing required column(s): ",
      paste(missing, collapse = ", ")
    )
  }
  
  if (!"normalized_name" %in% names(x)) {
    x$normalized_name <- .host_normalize_name(
      x$name
    )
  }
  
  x$TaxID <- as.character(x$TaxID)
  x$name <- as.character(x$name)
  x$name_class <- as.character(x$name_class)
  x$normalized_name <- as.character(x$normalized_name)
  
  x
}


.read_ncbi_ranked_taxonomy <- function(
    ranked_taxonomy_file = NULL
) {
  path <- .host_ranked_taxonomy_path(
    ranked_taxonomy_file
  )
  
  if (!file.exists(path)) {
    stop(
      "NCBI ranked taxonomy database not found:\n  ",
      path
    )
  }
  
  x <- readRDS(path)
  
  
  # ----------------------------------------------------------
  # Standardize column names from the current local taxonomy
  # database to the names used internally by aRborist
  # ----------------------------------------------------------
  
  rename_map <- c(
    tax_name = "Scientific.name",
    superkingdom = "Superkingdom",
    kingdom = "Kingdom",
    phylum = "Phylum",
    class = "Class",
    order = "Order",
    family = "Family",
    genus = "Genus",
    species = "Species"
  )
  
  for (old_name in names(rename_map)) {
    
    new_name <- rename_map[[old_name]]
    
    if (
      old_name %in% names(x) &&
      !new_name %in% names(x)
    ) {
      names(x)[names(x) == old_name] <- new_name
    }
  }
  
  
  # ----------------------------------------------------------
  # Confirm required columns
  # ----------------------------------------------------------
  
  required <- c(
    "TaxID",
    "Scientific.name"
  )
  
  missing <- setdiff(
    required,
    names(x)
  )
  
  if (length(missing) > 0) {
    stop(
      "ncbi_ranked_taxonomy.rds is missing required column(s): ",
      paste(missing, collapse = ", ")
    )
  }
  
  
  # ----------------------------------------------------------
  # Standardize TaxID
  # ----------------------------------------------------------
  
  x$TaxID <- as.character(x$TaxID)
  
  x
}

.read_shared_host_replacements <- function(
    shared_replacement_file = NULL
) {
  path <- .host_shared_replacement_path(
    shared_replacement_file
  )
  
  if (!file.exists(path)) {
    stop(
      "Shared host replacement file not found:\n  ",
      path
    )
  }
  
  x <- read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = character(0)
  )
  
  required <- c(
    "original_term",
    "replacement_term"
  )
  
  missing <- setdiff(
    required,
    names(x)
  )
  
  if (length(missing) > 0) {
    stop(
      "Shared replacement file is missing required column(s): ",
      paste(missing, collapse = ", ")
    )
  }
  
  x$original_term <- trimws(
    as.character(x$original_term)
  )
  
  x$replacement_term <- trimws(
    as.character(x$replacement_term)
  )
  
  x$normalized_original_term <- .host_normalize_name(
    x$original_term
  )
  
  # The same original term cannot point to different replacements.
  conflict_check <- x |>
    dplyr::filter(
      !is.na(normalized_original_term),
      normalized_original_term != ""
    ) |>
    dplyr::group_by(
      normalized_original_term
    ) |>
    dplyr::summarise(
      n_replacements = dplyr::n_distinct(
        replacement_term,
        na.rm = FALSE
      ),
      .groups = "drop"
    ) |>
    dplyr::filter(
      n_replacements > 1
    )
  
  if (nrow(conflict_check) > 0) {
    stop(
      "Conflicting duplicate terms were found in the shared host ",
      "replacement table.\nFirst examples:\n  ",
      paste(
        head(conflict_check$normalized_original_term, 10),
        collapse = "\n  "
      )
    )
  }
  
  x <- x[
    !duplicated(x$normalized_original_term),
    ,
    drop = FALSE
  ]
  
  x
}

.apply_shared_host_replacements <- function(
    meta,
    use_shared_replacements = TRUE,
    shared_replacement_file = NULL
) {
  if (!isTRUE(use_shared_replacements)) {
    message("Shared host replacements disabled.")
    return(meta)
  }
  
  replacements <- .read_shared_host_replacements(
    shared_replacement_file
  )
  
  current_key <- .host_normalize_name(
    meta$host.standardized
  )
  
  idx <- match(
    current_key,
    replacements$normalized_original_term
  )
  
  matched <- !is.na(idx)
  
  if (!any(matched)) {
    message("No terms matched the shared host replacement table.")
    return(meta)
  }
  
  replacement_values <- replacements$replacement_term[
    idx[matched]
  ]
  
  replacement_values[
    is.na(replacement_values) |
      trimws(replacement_values) == "" |
      toupper(trimws(replacement_values)) == "NA"
  ] <- NA_character_
  
  meta$host.standardized[matched] <- replacement_values
  
  message(
    "Shared replacement terms applied to ",
    sum(matched),
    " metadata row(s)."
  )
  
  meta
}

.lookup_local_taxid <- function(
    term,
    name_lookup
) {
  normalized_term <- .host_normalize_name(
    term
  )
  
  hits <- name_lookup[
    name_lookup$normalized_name == normalized_term,
    ,
    drop = FALSE
  ]
  
  if (nrow(hits) == 0) {
    return(list(
      status = "not_found",
      term = term,
      taxid = NA_character_,
      matched_name = NA_character_,
      matched_name_class = NA_character_,
      hits = hits
    ))
  }
  
  unique_taxids <- unique(
    hits$TaxID[
      !is.na(hits$TaxID) &
        hits$TaxID != ""
    ]
  )
  
  if (length(unique_taxids) == 0) {
    return(list(
      status = "not_found",
      term = term,
      taxid = NA_character_,
      matched_name = NA_character_,
      matched_name_class = NA_character_,
      hits = hits
    ))
  }
  
  if (length(unique_taxids) > 1) {
    return(list(
      status = "ambiguous",
      term = term,
      taxid = NA_character_,
      matched_name = NA_character_,
      matched_name_class = NA_character_,
      hits = hits
    ))
  }
  
  taxid <- unique_taxids[1]
  
  taxid_hits <- hits[
    hits$TaxID == taxid,
    ,
    drop = FALSE
  ]
  
  matched_names <- unique(
    taxid_hits$name[
      !is.na(taxid_hits$name) &
        taxid_hits$name != ""
    ]
  )
  
  matched_classes <- unique(
    taxid_hits$name_class[
      !is.na(taxid_hits$name_class) &
        taxid_hits$name_class != ""
    ]
  )
  
  list(
    status = "unique",
    term = term,
    taxid = taxid,
    matched_name = paste(
      matched_names,
      collapse = "; "
    ),
    matched_name_class = paste(
      matched_classes,
      collapse = "; "
    ),
    hits = taxid_hits
  )
}


# ============================================================
# Safe automatic host-name cleanup
# ============================================================

.propose_automatic_host_cleanup <- function(term) {
  
  term <- .host_clean_term(term)
  
  
  if (.host_is_invalid_tax_term(term)) {
    return(
      list(
        candidate = NA_character_,
        rule = NA_character_
      )
    )
  }
  
  
  original <- term
  rules <- character(0)
  
  
  # ----------------------------------------------------------
  # Rule 1:
  # Remove trailing botanical authority "L."
  #
  # Examples:
  #   Quercus robur L.        -> Quercus robur
  #   Heliconia aurantiaca L. -> Heliconia aurantiaca
  # ----------------------------------------------------------
  
  if (
    grepl(
      "\\s+L\\.$",
      term
    )
  ) {
    
    term <- sub(
      "\\s+L\\.$",
      "",
      term
    )
    
    term <- trimws(term)
    
    rules <- c(
      rules,
      "remove_trailing_L."
    )
  }
  
  
  # ----------------------------------------------------------
  # Rule 2:
  # Remove trailing "sp." or "sp"
  #
  # Examples:
  #   Sirex sp. -> Sirex
  #   Canna sp. -> Canna
  #
  # NOTE:
  # This merely generates a candidate.
  # The candidate is NOT accepted unless it uniquely matches
  # the local NCBI name database.
  # ----------------------------------------------------------
  
  if (
    grepl(
      "\\s+sp\\.?$",
      term,
      ignore.case = TRUE
    )
  ) {
    
    term <- sub(
      "\\s+sp\\.?$",
      "",
      term,
      ignore.case = TRUE
    )
    
    term <- trimws(term)
    
    rules <- c(
      rules,
      "remove_trailing_sp."
    )
  }
  
  
  # ----------------------------------------------------------
  # No transformation occurred
  # ----------------------------------------------------------
  
  if (
    identical(
      term,
      original
    ) ||
    !nzchar(term)
  ) {
    
    return(
      list(
        candidate = NA_character_,
        rule = NA_character_
      )
    )
  }
  
  
  list(
    candidate = term,
    rule = paste(
      rules,
      collapse = ";"
    )
  )
}


# ============================================================
# Build and validate automatic cleanup suggestions
#
# Validation requirement:
# transformed candidate must match exactly ONE TaxID in the
# local NCBI name database.
# ============================================================

.build_automatic_host_cleanup_map <- function(
    terms,
    name_lookup_file = NULL
) {
  
  terms <- unique(
    .host_clean_term(
      terms
    )
  )
  
  
  terms <- terms[
    !.host_is_invalid_tax_term(
      terms
    )
  ]
  
  
  if (length(terms) == 0) {
    
    return(
      data.frame(
        original_term = character(0),
        suggested_replacement = character(0),
        cleanup_rule = character(0),
        validation_status = character(0),
        TaxID = character(0),
        matched_name = character(0),
        matched_name_class = character(0),
        stringsAsFactors = FALSE
      )
    )
  }
  
  
  message(
    "Loading local NCBI name lookup for automatic host cleanup..."
  )
  
  
  name_lookup <- .read_ncbi_name_lookup(
    name_lookup_file
  )
  
  
  rows <- list()
  
  
  for (term in terms) {
    
    proposal <- .propose_automatic_host_cleanup(
      term
    )
    
    
    if (
      is.na(proposal$candidate) ||
      !nzchar(proposal$candidate)
    ) {
      next
    }
    
    
    candidate <- proposal$candidate
    
    
    # --------------------------------------------------------
    # Critical safety check:
    # candidate must uniquely resolve in local NCBI taxonomy
    # --------------------------------------------------------
    
    local_match <- .lookup_local_taxid(
      term = candidate,
      name_lookup = name_lookup
    )
    
    
    if (local_match$status == "unique") {
      
      validation_status <- "unique_ncbi_match"
      
      taxid <- local_match$taxid
      matched_name <- local_match$matched_name
      matched_name_class <- local_match$matched_name_class
      
      
    } else if (local_match$status == "ambiguous") {
      
      validation_status <- "ambiguous_ncbi_match"
      
      taxid <- NA_character_
      matched_name <- NA_character_
      matched_name_class <- NA_character_
      
      
    } else {
      
      validation_status <- "no_ncbi_match"
      
      taxid <- NA_character_
      matched_name <- NA_character_
      matched_name_class <- NA_character_
    }
    
    
    rows[[length(rows) + 1L]] <- data.frame(
      original_term = term,
      suggested_replacement = candidate,
      cleanup_rule = proposal$rule,
      validation_status = validation_status,
      TaxID = taxid,
      matched_name = matched_name,
      matched_name_class = matched_name_class,
      stringsAsFactors = FALSE
    )
  }
  
  
  if (length(rows) == 0) {
    
    return(
      data.frame(
        original_term = character(0),
        suggested_replacement = character(0),
        cleanup_rule = character(0),
        validation_status = character(0),
        TaxID = character(0),
        matched_name = character(0),
        matched_name_class = character(0),
        stringsAsFactors = FALSE
      )
    )
  }
  
  
  dplyr::bind_rows(
    rows
  )
}


# ============================================================
# Apply uniquely validated automatic cleanup to metadata
#
# This modifies ONLY host.standardized.
# It does not create another host-name column.
#
# host.standardized.original remains untouched.
# ============================================================

.apply_automatic_host_cleanup <- function(
    meta,
    project_name,
    host_dir = "./host_assessment",
    name_lookup_file = NULL,
    write_audit = TRUE
) {
  
  if (!"host.standardized" %in% names(meta)) {
    stop(
      "host.standardized must exist before automatic host cleanup."
    )
  }
  
  
  # ----------------------------------------------------------
  # Work only on UNIQUE terms.
  #
  # This is important for very large metadata tables.
  # We do NOT perform NCBI matching 1.5 million or 3.5 million
  # times.
  # ----------------------------------------------------------
  
  current_hosts <- .host_clean_term(
    meta$host.standardized
  )
  
  
  unique_terms <- unique(
    current_hosts[
      !.host_is_invalid_tax_term(
        current_hosts
      )
    ]
  )
  
  
  message(
    "Checking ",
    format(
      length(unique_terms),
      big.mark = ","
    ),
    " unique host terms for safe automatic cleanup."
  )
  
  
  cleanup_map <- .build_automatic_host_cleanup_map(
    terms = unique_terms,
    name_lookup_file = name_lookup_file
  )
  
  
  if (nrow(cleanup_map) == 0) {
    
    message(
      "No host terms matched the automatic cleanup patterns."
    )
    
    return(meta)
  }
  
  
  # ----------------------------------------------------------
  # Add accession counts for audit purposes
  # ----------------------------------------------------------
  
  term_counts <- table(
    current_hosts
  )
  
  
  cleanup_map$accession_count <- as.integer(
    term_counts[
      cleanup_map$original_term
    ]
  )
  
  
  cleanup_map$accession_count[
    is.na(cleanup_map$accession_count)
  ] <- 0L
  
  
  cleanup_map <- cleanup_map[
    order(
      -cleanup_map$accession_count,
      cleanup_map$original_term
    ),
    ,
    drop = FALSE
  ]
  
  
  # ----------------------------------------------------------
  # Write audit table
  #
  # This is outside the metadata table so it does not create
  # another host column.
  # ----------------------------------------------------------
  
  if (isTRUE(write_audit)) {
    
    if (!dir.exists(host_dir)) {
      dir.create(
        host_dir,
        recursive = TRUE
      )
    }
    
    
    audit_path <- file.path(
      host_dir,
      paste0(
        "host_automatic_cleanup_",
        project_name,
        ".csv"
      )
    )
    
    
    write.csv(
      cleanup_map,
      audit_path,
      row.names = FALSE
    )
    
    
    message(
      "Automatic host-cleanup audit written to: ",
      audit_path
    )
  }
  
  
  # ----------------------------------------------------------
  # Accept ONLY unique NCBI matches
  # ----------------------------------------------------------
  
  accepted <- cleanup_map[
    cleanup_map$validation_status ==
      "unique_ncbi_match",
    ,
    drop = FALSE
  ]
  
  
  if (nrow(accepted) == 0) {
    
    message(
      "No automatic cleanup candidates had unique NCBI matches."
    )
    
    return(meta)
  }
  
  
  # ----------------------------------------------------------
  # Vectorized mapping back onto metadata
  # ----------------------------------------------------------
  
  map_index <- match(
    current_hosts,
    accepted$original_term
  )
  
  
  rows_to_change <- !is.na(
    map_index
  )
  
  
  meta$host.standardized[
    rows_to_change
  ] <- accepted$suggested_replacement[
    map_index[
      rows_to_change
    ]
  ]
  
  
  changed_accessions <- sum(
    rows_to_change
  )
  
  
  message("")
  message(
    "Safe automatic host cleanup complete."
  )
  
  message(
    "  Unique terms changed: ",
    nrow(accepted)
  )
  
  message(
    "  Metadata rows changed: ",
    format(
      changed_accessions,
      big.mark = ","
    )
  )
  
  message(
    "  Ambiguous candidates left unchanged: ",
    sum(
      cleanup_map$validation_status ==
        "ambiguous_ncbi_match"
    )
  )
  
  message(
    "  Candidates without NCBI matches left unchanged: ",
    sum(
      cleanup_map$validation_status ==
        "no_ncbi_match"
    )
  )
  
  
  meta
}


.lookup_local_ranked_taxonomy <- function(
    term,
    taxid,
    matched_name,
    matched_name_class,
    ranked_taxonomy
) {
  hit <- ranked_taxonomy[
    ranked_taxonomy$TaxID == as.character(taxid),
    ,
    drop = FALSE
  ]
  
  if (nrow(hit) == 0) {
    return(NULL)
  }
  
  hit <- hit[1, , drop = FALSE]
  
  get_value <- function(col) {
    if (!col %in% names(hit)) {
      return(NA_character_)
    }
    
    x <- as.character(hit[[col]][1])
    
    if (is.na(x) || trimws(x) == "") {
      NA_character_
    } else {
      x
    }
  }
  
  data.frame(
    Host.standardized = term,
    Host.taxid = as.character(taxid),
    Host.matched_name = matched_name,
    Host.matched_name_class = matched_name_class,
    Host.lookup_source = "local_ncbi",
    Host.superkingdom = get_value("Superkingdom"),
    Host.kingdom = get_value("Kingdom"),
    Host.phylum = get_value("Phylum"),
    Host.class = get_value("Class"),
    Host.order = get_value("Order"),
    Host.family = get_value("Family"),
    Host.genus = get_value("Genus"),
    Host.species = get_value("Species"),
    stringsAsFactors = FALSE
  )
}


.empty_host_ambiguity <- function() {
  data.frame(
    original_term = character(0),
    normalized_name = character(0),
    TaxID = character(0),
    matched_name = character(0),
    name_class = character(0),
    scientific_name = character(0),
    stringsAsFactors = FALSE
  )
}


.read_ambiguous_terms <- function(
    project_name,
    host_dir = "./host_assessment"
) {
  path <- .host_ambiguous_path(
    project_name,
    host_dir
  )
  
  if (!file.exists(path)) {
    return(
      .empty_host_ambiguity()
    )
  }
  
  x <- read.csv(
    path,
    stringsAsFactors = FALSE
  )
  
  # Keep TaxID type consistent with the rest of the host pipeline
  if ("TaxID" %in% names(x)) {
    x$TaxID <- as.character(x$TaxID)
  }
  
  x
}


.write_ambiguous_terms <- function(
    ambiguity_df,
    project_name,
    host_dir = "./host_assessment"
) {
  path <- .host_ambiguous_path(
    project_name,
    host_dir
  )
  
  if (nrow(ambiguity_df) > 0) {
    ambiguity_df <- unique(
      ambiguity_df
    )
    
    ambiguity_df <- ambiguity_df[
      order(
        ambiguity_df$original_term,
        ambiguity_df$TaxID
      ),
      ,
      drop = FALSE
    ]
  }
  
  write.csv(
    ambiguity_df,
    path,
    row.names = FALSE
  )
  
  message(
    "Ambiguous host terms file written to: ",
    path
  )
  
  invisible(
    ambiguity_df
  )
}


.append_ambiguous_lookup <- function(
    ambiguity_df,
    term,
    name_hits,
    ranked_taxonomy
) {
  if (nrow(name_hits) == 0) {
    return(ambiguity_df)
  }
  
  name_hits <- name_hits |>
    dplyr::distinct(
      TaxID,
      name,
      name_class,
      .keep_all = TRUE
    )
  
  scientific_lookup <- ranked_taxonomy |>
    dplyr::select(
      TaxID,
      Scientific.name
    )
  
  scientific_lookup$TaxID <- as.character(
    scientific_lookup$TaxID
  )
  
  detail <- name_hits |>
    dplyr::transmute(
      original_term = term,
      normalized_name = .host_normalize_name(term),
      TaxID = as.character(TaxID),
      matched_name = as.character(name),
      name_class = as.character(name_class)
    ) |>
    dplyr::left_join(
      scientific_lookup,
      by = "TaxID"
    ) |>
    dplyr::rename(
      scientific_name = Scientific.name
    )
  
  dplyr::bind_rows(
    ambiguity_df,
    detail
  ) |>
    dplyr::distinct(
      original_term,
      TaxID,
      matched_name,
      name_class,
      .keep_all = TRUE
    )
}



# -----------------------------
# Merge taxonomy into metadata
# -----------------------------

merge_host_taxonomy_into_metadata <- function(
    project_name,
    host_dir = "./host_assessment",
    metadata_file = NULL
) {
  if (is.null(metadata_file)) {
    metadata_file <- .host_metadata_path(project_name)
  }
  
  if (!file.exists(metadata_file)) {
    stop("Metadata file not found: ", metadata_file)
  }
  
  taxonomy_path <- .host_taxonomy_path(project_name, host_dir)
  
  if (!file.exists(taxonomy_path)) {
    stop("Host taxonomy file not found: ", taxonomy_path)
  }
  
  meta <- read.csv(metadata_file, stringsAsFactors = FALSE)
  meta <- .ensure_host_standardized(meta)
  
  host_taxonomy <- .read_host_taxonomy(project_name, host_dir)
  
  taxonomy_cols <- .host_taxonomy_cols()
  taxonomy_value_cols <- setdiff(taxonomy_cols, "Host.standardized")
  
  for (col in taxonomy_value_cols) {
    if (col %in% names(meta)) {
      meta[[col]] <- NULL
    }
  }
  
  match_idx <- match(meta$host.standardized, host_taxonomy$Host.standardized)
  
  for (col in taxonomy_value_cols) {
    meta[[col]] <- host_taxonomy[[col]][match_idx]
  }
  
  write.csv(meta, metadata_file, row.names = FALSE)
  
  message("Host taxonomy merged into metadata: ", metadata_file)
  
  invisible(meta)
}


# -----------------------------
# Initial pass
# -----------------------------
run_host_assessment_initial_pass <- function(
    project_name,
    host_dir = "./host_assessment",
    db = "ncbi",
    overwrite = FALSE,
    metadata_file = NULL,
    use_isolation_source = FALSE,
    overwrite_host_standardized = FALSE,
    use_shared_replacements = TRUE,
    shared_replacement_file = NULL,
    use_automatic_cleanup = TRUE,
    use_taxize_fallback = TRUE,
    name_lookup_file = NULL,
    ranked_taxonomy_file = NULL
) {
  
  if (!dir.exists(host_dir)) {
    dir.create(
      host_dir,
      recursive = TRUE
    )
  }
  
  
  if (is.null(metadata_file)) {
    metadata_file <- .host_metadata_path(
      project_name
    )
  }
  
  
  if (!file.exists(metadata_file)) {
    stop(
      "Metadata file not found: ",
      metadata_file
    )
  }
  
  
  message(
    "Running host assessment initial pass..."
  )
  
  
  # ============================================================
  # 1. Read metadata
  # ============================================================
  
  meta <- read.csv(
    metadata_file,
    stringsAsFactors = FALSE
  )
  
  
  # ============================================================
  # 2. Initialize host.standardized
  #
  # host.standardized.original is preserved by
  # .ensure_host_standardized().
  # ============================================================
  
  meta <- .ensure_host_standardized(
    meta,
    use_isolation_source = use_isolation_source,
    overwrite_host_standardized =
      overwrite_host_standardized
  )
  
  
  # ============================================================
  # 3. Apply shared reusable replacements
  #
  # These remain the highest-priority curated replacements.
  # ============================================================
  
  meta <- .apply_shared_host_replacements(
    meta = meta,
    use_shared_replacements =
      use_shared_replacements,
    shared_replacement_file =
      shared_replacement_file
  )
  
  
  # ============================================================
  # 4. Safe automatic cleanup
  #
  # Current automatic rules:
  #
  #   Quercus robur L. -> Quercus robur
  #   Sirex sp.        -> Sirex
  #
  # A transformed candidate is accepted ONLY when it maps to
  # exactly one TaxID in the local NCBI name database.
  #
  # Failed or ambiguous transformed strings are left completely
  # unchanged.
  # ============================================================
  
  if (isTRUE(use_automatic_cleanup)) {
    
    meta <- .apply_automatic_host_cleanup(
      meta = meta,
      project_name = project_name,
      host_dir = host_dir,
      name_lookup_file = name_lookup_file,
      write_audit = TRUE
    )
    
  } else {
    
    message(
      "Automatic host-name cleanup disabled."
    )
  }
  
  
  # ============================================================
  # 5. Save standardized metadata
  # ============================================================
  
  write.csv(
    meta,
    metadata_file,
    row.names = FALSE
  )
  
  
  # ============================================================
  # 6. Extract unique valid standardized host terms
  # ============================================================
  
  host_terms <- unique(
    .host_clean_term(
      meta$host.standardized
    )
  )
  
  
  host_terms <- host_terms[
    !.host_is_invalid_tax_term(
      host_terms
    )
  ]
  
  
  terms_path <- file.path(
    host_dir,
    paste0(
      "host_terms_for_taxonomy_",
      project_name,
      ".csv"
    )
  )
  
  
  write.csv(
    data.frame(
      host = host_terms,
      stringsAsFactors = FALSE
    ),
    terms_path,
    row.names = FALSE
  )
  
  
  message(
    "Unique standardized host terms: ",
    length(host_terms)
  )
  
  
  message(
    "Host terms file written to: ",
    terms_path
  )
  
  
  # ============================================================
  # 7. Local NCBI taxonomy lookup
  #
  # The local name database already contains multiple NCBI name
  # classes. Therefore scientific names, synonyms, common names,
  # etc. can resolve here when they uniquely identify one TaxID.
  #
  # taxize remains an optional fallback.
  # ============================================================
  
  lookup_result <- .search_host_terms(
    terms = host_terms,
    project_name = project_name,
    host_dir = host_dir,
    db = db,
    overwrite = overwrite,
    term_type = "initial_lookup",
    use_taxize_fallback =
      use_taxize_fallback,
    name_lookup_file =
      name_lookup_file,
    ranked_taxonomy_file =
      ranked_taxonomy_file
  )
  
  
  # ============================================================
  # 8. Merge successful taxonomy into metadata
  # ============================================================
  
  merge_host_taxonomy_into_metadata(
    project_name = project_name,
    host_dir = host_dir,
    metadata_file = metadata_file
  )
  
  
  message(
    "Initial host assessment pass completed."
  )
  
  
  invisible(
    lookup_result
  )
}




# -----------------------------
# Refinement pass
# -----------------------------

.has_usable_host_taxonomy <- function(host_taxonomy) {
  rank_cols <- c(
    "Host.kingdom",
    "Host.phylum",
    "Host.class",
    "Host.order",
    "Host.family",
    "Host.genus",
    "Host.species"
  )
  
  rank_cols <- rank_cols[rank_cols %in% names(host_taxonomy)]
  
  if (length(rank_cols) == 0) {
    return(rep(FALSE, nrow(host_taxonomy)))
  }
  
  apply(
    host_taxonomy[, rank_cols, drop = FALSE],
    1,
    function(x) any(!is.na(x) & trimws(as.character(x)) != "")
  )
}

run_host_assessment_refinement_pass <- function(
    project_name,
    host_dir = "./host_assessment",
    db = "ncbi",
    overwrite = FALSE,
    metadata_file = NULL,
    use_taxize_fallback = TRUE,
    name_lookup_file = NULL,
    ranked_taxonomy_file = NULL
) {
  
  if (!dir.exists(host_dir)) {
    dir.create(
      host_dir,
      recursive = TRUE
    )
  }
  
  
  if (is.null(metadata_file)) {
    metadata_file <- .host_metadata_path(
      project_name
    )
  }
  
  
  if (!file.exists(metadata_file)) {
    stop(
      "Metadata file not found: ",
      metadata_file
    )
  }
  
  
  failed_path <- .host_failed_path(
    project_name,
    host_dir
  )
  
  
  if (!file.exists(failed_path)) {
    stop(
      "Failed terms file not found: ",
      failed_path,
      "\nRun run_host_assessment_initial_pass() first."
    )
  }
  
  
  message(
    "Running host assessment refinement pass..."
  )
  
  
  # ============================================================
  # 1. Read failed terms and existing taxonomy
  # ============================================================
  
  failed_df <- .read_failed_terms(
    project_name,
    host_dir
  )
  
  
  host_taxonomy <- .read_host_taxonomy(
    project_name,
    host_dir
  )
  
  
  replacement_ok <- !.host_is_invalid_tax_term(
    failed_df$replacement_term
  )
  
  
  replacement_rows <- failed_df[
    replacement_ok,
    ,
    drop = FALSE
  ]
  
  
  # ============================================================
  # Nothing supplied
  # ============================================================
  
  if (nrow(replacement_rows) == 0) {
    
    message(
      "No replacement terms supplied yet. Nothing to refine."
    )
    
    
    failed_df <- .add_failed_host_term_counts(
      failed_df,
      project_name,
      metadata_file
    )
    
    
    .write_failed_terms(
      failed_df,
      project_name,
      host_dir
    )
    
    
    merge_host_taxonomy_into_metadata(
      project_name = project_name,
      host_dir = host_dir,
      metadata_file = metadata_file
    )
    
    
    return(
      invisible(
        list(
          taxonomy = host_taxonomy,
          failed_table = failed_df
        )
      )
    )
  }
  
  
  # ============================================================
  # 2. Determine replacement terms requiring taxonomy lookup
  # ============================================================
  
  replacement_terms <- unique(
    .host_clean_term(
      replacement_rows$replacement_term
    )
  )
  
  
  replacement_terms <- replacement_terms[
    !.host_is_invalid_tax_term(
      replacement_terms
    )
  ]
  
  
  usable_taxonomy_rows <- .has_usable_host_taxonomy(
    host_taxonomy
  )
  
  
  already_have_taxonomy <- unique(
    .host_clean_term(
      host_taxonomy$Host.standardized[
        usable_taxonomy_rows
      ]
    )
  )
  
  
  if (overwrite) {
    
    terms_to_query <- replacement_terms
    
  } else {
    
    terms_to_query <- setdiff(
      replacement_terms,
      already_have_taxonomy
    )
  }
  
  
  message(
    "Replacement terms supplied: ",
    length(replacement_terms)
  )
  
  
  message(
    "Replacement terms already have usable taxonomy: ",
    length(
      intersect(
        replacement_terms,
        already_have_taxonomy
      )
    )
  )
  
  
  message(
    "Replacement terms to query this pass: ",
    length(terms_to_query)
  )
  
  
  # ============================================================
  # 3. Look up replacement taxonomy
  # ============================================================
  
  if (length(terms_to_query) > 0) {
    
    lookup_result <- .search_host_terms(
      terms = terms_to_query,
      project_name = project_name,
      host_dir = host_dir,
      db = db,
      overwrite = overwrite,
      term_type = "replacement_lookup",
      use_taxize_fallback = use_taxize_fallback,
      name_lookup_file = name_lookup_file,
      ranked_taxonomy_file = ranked_taxonomy_file
    )
    
  } else {
    
    lookup_result <- list(
      taxonomy = host_taxonomy,
      failed_table = failed_df,
      newly_failed = character(0),
      newly_successful = character(0)
    )
  }
  
  
  # ============================================================
  # 4. Reload taxonomy after replacement lookups
  # ============================================================
  
  host_taxonomy <- .read_host_taxonomy(
    project_name,
    host_dir
  )
  
  
  usable_taxonomy_rows <- .has_usable_host_taxonomy(
    host_taxonomy
  )
  
  
  successful_taxonomy_terms <- unique(
    .host_clean_term(
      host_taxonomy$Host.standardized[
        usable_taxonomy_rows
      ]
    )
  )
  
  
  # ============================================================
  # 5. Read large metadata table
  #
  # fread is substantially faster than read.csv for this dataset.
  # ============================================================
  
  message(
    "Reading metadata for replacement application..."
  )
  
  
  meta <- data.table::fread(
    metadata_file,
    data.table = FALSE
  )
  
  
  meta <- .ensure_host_standardized(
    meta
  )
  
  
  # ============================================================
  # 6. Build validated replacement map
  #
  # ONLY replacements with successful taxonomy are allowed.
  # ============================================================
  
  replacement_map <- replacement_rows
  
  
  replacement_map$original_key <- .host_display_term(
    replacement_map$original_term
  )
  
  
  replacement_map$replacement_clean <- .host_clean_term(
    replacement_map$replacement_term
  )
  
  
  replacement_map <- replacement_map[
    !.host_is_invalid_tax_term(
      replacement_map$replacement_clean
    ) &
      replacement_map$replacement_clean %in%
      successful_taxonomy_terms,
    ,
    drop = FALSE
  ]
  
  
  # ------------------------------------------------------------
  # Make sure one original term cannot point to two replacements
  # ------------------------------------------------------------
  
  if (nrow(replacement_map) > 0) {
    
    conflict_check <- replacement_map |>
      dplyr::group_by(
        original_key
      ) |>
      dplyr::summarise(
        n_replacements = dplyr::n_distinct(
          replacement_clean
        ),
        .groups = "drop"
      ) |>
      dplyr::filter(
        n_replacements > 1
      )
    
    
    if (nrow(conflict_check) > 0) {
      
      stop(
        "Conflicting replacement terms were found for the same original host term:\n  ",
        paste(
          head(
            conflict_check$original_key,
            20
          ),
          collapse = "\n  "
        )
      )
    }
    
    
    replacement_map <- replacement_map[
      !duplicated(
        replacement_map$original_key
      ),
      ,
      drop = FALSE
    ]
  }
  
  
  message(
    "Validated replacement mappings available: ",
    nrow(replacement_map)
  )
  
  
  # ============================================================
  # 7. Apply replacements VECTORIZED
  #
  # This replaces the old per-replacement/per-1.5-million-row
  # loop.
  # ============================================================
  
  message(
    "Applying successful replacement terms to metadata..."
  )
  
  
  if (nrow(replacement_map) > 0) {
    
    current_key <- .host_display_term(
      meta$host.standardized
    )
    
    
    original_key <- if (
      "host.standardized.original" %in% names(meta)
    ) {
      
      .host_display_term(
        meta$host.standardized.original
      )
      
    } else {
      
      current_key
    }
    
    
    # First prefer the CURRENT standardized value.
    current_match <- match(
      current_key,
      replacement_map$original_key
    )
    
    
    # If that does not match, try the preserved original value.
    original_match <- match(
      original_key,
      replacement_map$original_key
    )
    
    
    replacement_index <- current_match
    
    
    use_original_match <- is.na(
      replacement_index
    ) &
      !is.na(
        original_match
      )
    
    
    replacement_index[
      use_original_match
    ] <- original_match[
      use_original_match
    ]
    
    
    rows_to_change <- !is.na(
      replacement_index
    )
    
    
    n_rows_changed <- sum(
      rows_to_change
    )
    
    
    n_terms_used <- length(
      unique(
        replacement_index[
          rows_to_change
        ]
      )
    )
    
    
    if (n_rows_changed > 0) {
      
      meta$host.standardized[
        rows_to_change
      ] <- replacement_map$replacement_clean[
        replacement_index[
          rows_to_change
        ]
      ]
    }
    
    
    message(
      "Replacement application complete."
    )
    
    
    message(
      "  Unique replacement mappings used: ",
      n_terms_used
    )
    
    
    message(
      "  Metadata rows changed: ",
      format(
        n_rows_changed,
        big.mark = ","
      )
    )
    
  } else {
    
    message(
      "No validated replacement mappings were available to apply."
    )
  }
  
  
  # ============================================================
  # 8. Update failed-term statuses
  # ============================================================
  
  failed_df <- .read_failed_terms(
    project_name,
    host_dir
  )
  
  
  host_taxonomy <- .read_host_taxonomy(
    project_name,
    host_dir
  )
  
  
  usable_taxonomy_rows <- .has_usable_host_taxonomy(
    host_taxonomy
  )
  
  
  successful_taxonomy_terms <- unique(
    .host_clean_term(
      host_taxonomy$Host.standardized[
        usable_taxonomy_rows
      ]
    )
  )
  
  
  replacement_clean <- .host_clean_term(
    failed_df$replacement_term
  )
  
  
  has_replacement <- !.host_is_invalid_tax_term(
    replacement_clean
  )
  
  
  replacement_success <- has_replacement &
    replacement_clean %in%
    successful_taxonomy_terms
  
  
  replacement_failed <- has_replacement &
    !replacement_success
  
  
  failed_df$replacement_lookup_status[
    replacement_success
  ] <- "replacement_successful"
  
  
  failed_df$lookup_status[
    replacement_success
  ] <- "resolved_by_replacement"
  
  
  failed_df$replacement_lookup_status[
    replacement_failed
  ] <- "replacement_failed"
  
  
  # ============================================================
  # 9. Add failed replacement terms if needed
  # ============================================================
  
  failed_replacements <- unique(
    replacement_clean[
      replacement_failed
    ]
  )
  
  
  failed_replacements <- failed_replacements[
    !.host_is_invalid_tax_term(
      failed_replacements
    )
  ]
  
  
  if (length(failed_replacements) > 0) {
    
    message(
      "Replacement terms that still failed taxonomy lookup:"
    )
    
    
    message(
      "  - ",
      paste(
        failed_replacements,
        collapse = "\n  - "
      )
    )
    
    
    failed_df <- .append_failed_terms(
      failed_df,
      terms = failed_replacements,
      term_type = "replacement_lookup",
      lookup_status = "failed"
    )
  }
  
  
  # ============================================================
  # 10. Update failed-term accession counts IN MEMORY
  #
  # Avoid rereading the giant metadata CSV.
  # ============================================================
  
  message(
    "Updating failed-term accession counts..."
  )
  
  
  host_display <- .host_display_term(
    meta$host.standardized
  )
  
  
  count_table <- as.data.frame(
    table(
      host_display
    ),
    stringsAsFactors = FALSE
  )
  
  
  names(count_table) <- c(
    "original_term",
    "accession_count"
  )
  
  
  failed_df <- .ensure_failed_columns(
    failed_df
  )
  
  
  failed_df$accession_count <- count_table$accession_count[
    match(
      failed_df$original_term,
      count_table$original_term
    )
  ]
  
  
  failed_df$accession_count[
    is.na(
      failed_df$accession_count
    )
  ] <- 0L
  
  
  .write_failed_terms(
    failed_df,
    project_name,
    host_dir
  )
  
  
  # ============================================================
  # 11. Merge taxonomy into metadata IN MEMORY
  #
  # Avoid writing metadata, rereading it, then writing it again.
  # ============================================================
  
  message(
    "Merging host taxonomy into metadata..."
  )
  
  
  host_taxonomy <- .read_host_taxonomy(
    project_name,
    host_dir
  )
  
  
  taxonomy_cols <- .host_taxonomy_cols()
  
  
  taxonomy_value_cols <- setdiff(
    taxonomy_cols,
    "Host.standardized"
  )
  
  
  host_key <- .host_clean_term(
    meta$host.standardized
  )
  
  
  taxonomy_key <- .host_clean_term(
    host_taxonomy$Host.standardized
  )
  
  
  taxonomy_match <- match(
    host_key,
    taxonomy_key
  )
  
  
  for (col in taxonomy_value_cols) {
    
    meta[[col]] <- host_taxonomy[[col]][
      taxonomy_match
    ]
  }
  
  
  # ============================================================
  # 12. Write metadata ONCE
  # ============================================================
  
  message(
    "Writing updated metadata to:"
  )
  
  message(
    "  ",
    metadata_file
  )
  
  
  data.table::fwrite(
    meta,
    metadata_file,
    na = "NA"
  )
  
  
  message(
    "Updated metadata written successfully."
  )
  
  
  message(
    "Refinement pass completed."
  )
  
  
  invisible(
    list(
      taxonomy = host_taxonomy,
      failed_table = failed_df,
      lookup_result = lookup_result
    )
  )
}


# -----------------------------
# Summary wrapper
# -----------------------------

run_host_assessment_summary <- function(
    project_name,
    fungal_rank = "genus",
    host_rank = "phylum",
    keep_NAs = FALSE,
    host_dir = "./host_assessment",
    metadata_file = NULL
) {
  if (is.null(metadata_file)) {
    metadata_file <- .host_metadata_path(project_name)
  }
  
  host_col <- paste0("Host.", host_rank)
  
  meta <- read.csv(metadata_file, stringsAsFactors = FALSE)
  
  if (!host_col %in% names(meta)) {
    message("Requested column ", host_col, " not found. Attempting to merge host taxonomy into metadata first.")
    
    merge_host_taxonomy_into_metadata(
      project_name = project_name,
      host_dir = host_dir,
      metadata_file = metadata_file
    )
  }
  
  summarize_host_usage(
    project_name = project_name,
    fungal_rank = fungal_rank,
    host_rank = host_rank,
    keep_NAs = keep_NAs,
    host_dir = host_dir,
    metadata_file = metadata_file
  )
}


# ============================================================
# Advanced host-name extraction / cleanup
#
# Purpose:
#   Recover valid NCBI taxa from messy host metadata such as:
#
#   cucumber fruit (Cucumis sativus L.)
#       -> Cucumis sativus
#
#   on dead branches of Camellia sinensis
#       -> Camellia sinensis
#
#   pearl millet stem
#       -> Cenchrus americanus, if "pearl millet" uniquely
#          resolves to that TaxID in the local NCBI name table
#
# Safety:
#   - Candidates are generated liberally.
#   - Automatic replacement occurs ONLY when all valid
#     candidates resolve to exactly one TaxID.
#   - Multiple taxa are never automatically collapsed.
#   - Uncertain strings containing "or" or "(?)" are review-only.
#
# This does NOT create any additional host columns.
# ============================================================


.host_candidate_clean_text <- function(x) {
  
  x <- as.character(x)
  
  # HTML line breaks
  x <- gsub(
    "(?i)<br\\s*/?>",
    " ",
    x,
    perl = TRUE
  )
  
  # Non-breaking spaces
  x <- gsub(
    "\u00A0",
    " ",
    x,
    fixed = TRUE
  )
  
  # Treat some obvious accidental separators as spaces
  # Example:
  #   Panax?notoginseng
  # becomes:
  #   Panax notoginseng
  x <- gsub(
    "[?+|/\\\\]",
    " ",
    x
  )
  
  # Remove grouping punctuation while retaining the text inside
  x <- gsub(
    "[()\\[\\]\\{\\}\"']",
    " ",
    x,
    perl = TRUE
  )
  
  # Other separators
  x <- gsub(
    "[:;,]",
    " ",
    x
  )
  
  # Collapse whitespace
  x <- gsub(
    "[[:space:]]+",
    " ",
    x
  )
  
  trimws(x)
}


.host_candidate_segments <- function(x) {
  
  x <- as.character(x)
  
  # Split obvious multi-host / multi-part strings.
  #
  # Examples:
  #   Cedrus deodara and Pinus wallichiana
  #   Quercus alba or Pinus strobus
  #   Sphagnum, Betula & Picea spp.
  
  parts <- unlist(
    strsplit(
      x,
      "(?i)\\s*(?:;|,|&|\\band\\b|\\bor\\b|<br\\s*/?>)\\s*",
      perl = TRUE
    )
  )
  
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  
  unique(parts)
}


.host_descriptor_candidates <- function(x) {
  
  x <- .host_candidate_clean_text(x)
  
  if (!nzchar(x)) {
    return(character(0))
  }
  
  candidates <- x
  
  
  # ----------------------------------------------------------
  # Remove trailing cultivar notation
  # ----------------------------------------------------------
  
  y <- gsub(
    "(?i)\\s+(?:cultivar|cv\\.?)\\s+.*$",
    "",
    x,
    perl = TRUE
  )
  
  y <- trimws(y)
  
  if (nzchar(y) && y != x) {
    candidates <- c(
      candidates,
      y
    )
  }
  
  
  # ----------------------------------------------------------
  # Progressively strip sample / tissue descriptors
  #
  # Keeping intermediate versions is important.
  #
  # Curry leaf plant
  #
  # gives:
  #   Curry leaf plant
  #   Curry leaf
  #   Curry
  #
  # If "Curry leaf" is an NCBI common name, that version can
  # succeed without incorrectly reducing it all the way to Curry.
  # ----------------------------------------------------------
  
  suffix_pattern <- paste0(
    "(?i)\\s+(?:",
    paste(
      c(
        "roots?",
        "leaf",
        "leaves",
        "stems?",
        "fruits?",
        "surface",
        "bark",
        "branches?",
        "bulbs?",
        "spore[- ]capsules?",
        "ears?",
        "plants?"
      ),
      collapse = "|"
    ),
    ")$"
  )
  
  current <- x
  
  for (i in seq_len(5)) {
    
    new_value <- sub(
      suffix_pattern,
      "",
      current,
      perl = TRUE
    )
    
    new_value <- trimws(
      new_value
    )
    
    if (
      !nzchar(new_value) ||
      new_value == current
    ) {
      break
    }
    
    candidates <- c(
      candidates,
      new_value
    )
    
    current <- new_value
  }
  
  
  # ----------------------------------------------------------
  # Simple leading contextual phrases
  # ----------------------------------------------------------
  
  leading_patterns <- c(
    "(?i)^isolated from\\s+",
    "(?i)^collected from\\s+",
    "(?i)^obtained from\\s+",
    "(?i)^hosted on\\s+",
    "(?i)^from\\s+",
    "(?i)^on\\s+"
  )
  
  for (pat in leading_patterns) {
    
    new_value <- sub(
      pat,
      "",
      x,
      perl = TRUE
    )
    
    new_value <- trimws(
      new_value
    )
    
    if (
      nzchar(new_value) &&
      new_value != x
    ) {
      candidates <- c(
        candidates,
        new_value
      )
    }
  }
  
  
  unique(
    candidates[
      nzchar(candidates)
    ]
  )
}


.generate_host_candidates_one <- function(
    term,
    max_ngram = 5
) {
  
  term <- trimws(
    as.character(term)
  )
  
  rows <- list()
  
  
  add_candidate <- function(
      candidate,
      method,
      priority
  ) {
    
    candidate <- trimws(
      as.character(candidate)
    )
    
    if (
      is.na(candidate) ||
      !nzchar(candidate) ||
      toupper(candidate) == "NA"
    ) {
      return(invisible(NULL))
    }
    
    rows[[length(rows) + 1L]] <<- data.frame(
      original_term = term,
      candidate = candidate,
      candidate_method = method,
      method_priority = priority,
      stringsAsFactors = FALSE
    )
    
    invisible(NULL)
  }
  
  
  cleaned <- .host_candidate_clean_text(
    term
  )
  
  
  # ==========================================================
  # A. Embedded binomial / infraspecific names
  # ==========================================================
  
  # Examples:
  #
  #   cucumber fruit Cucumis sativus L.
  #                  ^^^^^^^^^^^^^^^^
  #
  #   Taxus chinensis Pilg. Rehder
  #   ^^^^^^^^^^^^^^^
  
  binomial_pattern <- paste0(
    "\\b",
    "[A-Z][A-Za-z-]{2,}",
    "\\s+",
    "[a-z][A-Za-z-]{1,}",
    "(?:",
      "\\s+",
      "(?:subsp\\.?|ssp\\.?|var\\.?|f\\.?)",
      "\\s+",
      "[A-Za-z-]{1,}",
    ")?"
  )
  
  binomial_hits <- stringr::str_extract_all(
    cleaned,
    binomial_pattern
  )[[1]]
  
  binomial_hits <- unique(
    binomial_hits[
      !is.na(binomial_hits) &
        nzchar(binomial_hits)
    ]
  )
  
  if (length(binomial_hits) > 0) {
    
    for (candidate in binomial_hits) {
      
      add_candidate(
        candidate,
        "embedded_scientific_name",
        1
      )
    }
  }
  
  
  # ==========================================================
  # B. Genus sp. / genus spp.
  # ==========================================================
  
  genus_sp_pattern <- paste0(
    "\\b",
    "([A-Z][A-Za-z-]{2,})",
    "\\s+",
    "spp?\\.?",
    "\\b"
  )
  
  genus_sp_hits <- stringr::str_match_all(
    cleaned,
    genus_sp_pattern
  )[[1]]
  
  if (
    !is.null(genus_sp_hits) &&
    nrow(genus_sp_hits) > 0
  ) {
    
    genus_values <- unique(
      genus_sp_hits[, 2]
    )
    
    genus_values <- genus_values[
      !is.na(genus_values) &
        nzchar(genus_values)
    ]
    
    for (candidate in genus_values) {
      
      add_candidate(
        candidate,
        "embedded_genus_sp",
        2
      )
    }
  }
  
  
  # ==========================================================
  # C. Genus followed by "hybrid"
  #
  # Vitis hybrid cultivar; Prairie Star
  # -> Vitis
  # ==========================================================
  
  hybrid_hits <- stringr::str_match_all(
    cleaned,
    "\\b([A-Z][A-Za-z-]{2,})\\s+hybrid\\b"
  )[[1]]
  
  if (
    !is.null(hybrid_hits) &&
    nrow(hybrid_hits) > 0
  ) {
    
    hybrid_values <- unique(
      hybrid_hits[, 2]
    )
    
    hybrid_values <- hybrid_values[
      !is.na(hybrid_values) &
        nzchar(hybrid_values)
    ]
    
    for (candidate in hybrid_values) {
      
      add_candidate(
        candidate,
        "hybrid_genus",
        2
      )
    }
  }
  
  
  # ==========================================================
  # D. Single scientific-looking taxon after context words
  #
  # on dead Bambusa
  # -> Bambusa
  #
  # ... of Salix
  # -> Salix
  # ==========================================================
  
  context_pattern <- paste0(
    "(?:(?i:on|of|under|with|from))",
    "\\s+",
    "(?:dead\\s+)?",
    "(?:branches?\\s+of\\s+)?",
    "([A-Z][A-Za-z-]{2,})",
    "\\b"
  )
  
  context_hits <- stringr::str_match_all(
    cleaned,
    context_pattern
  )[[1]]
  
  if (
    !is.null(context_hits) &&
    nrow(context_hits) > 0
  ) {
    
    context_values <- unique(
      context_hits[, 2]
    )
    
    context_values <- context_values[
      !is.na(context_values) &
        nzchar(context_values)
    ]
    
    for (candidate in context_values) {
      
      add_candidate(
        candidate,
        "context_single_taxon",
        3
      )
    }
  }
  
  
  # ==========================================================
  # E. Split multi-part strings
  #
  # Allows detection of:
  #
  #   Tsuga canadensis or Betula
  #
  # where Betula would otherwise be missed.
  # ==========================================================
  
  segments <- .host_candidate_segments(
    term
  )
  
  for (segment in segments) {
    
    segment_clean <- .host_candidate_clean_text(
      segment
    )
    
    segment_binomial <- stringr::str_extract_all(
      segment_clean,
      binomial_pattern
    )[[1]]
    
    segment_genus_sp <- stringr::str_match_all(
      segment_clean,
      genus_sp_pattern
    )[[1]]
    
    has_segment_taxon <-
      length(segment_binomial) > 0 ||
      (
        !is.null(segment_genus_sp) &&
        nrow(segment_genus_sp) > 0
      )
    
    
    # If an entire segment is just one capitalized word,
    # treat it as a possible taxon.
    #
    # Example:
    #
    #   Tsuga canadensis OR Betula
    #                       ^^^^^^
    
    if (
      !has_segment_taxon &&
      grepl(
        "^[A-Z][A-Za-z-]{2,}$",
        segment_clean
      )
    ) {
      
      add_candidate(
        segment_clean,
        "bare_segment_taxon",
        3
      )
    }
  }
  
  
  scientific_style_methods <- c(
    "embedded_scientific_name",
    "embedded_genus_sp",
    "hybrid_genus",
    "context_single_taxon",
    "bare_segment_taxon"
  )
  
  currently_scientific <- FALSE
  
  if (length(rows) > 0) {
    
    currently_scientific <- any(
      vapply(
        rows,
        function(x) {
          x$candidate_method[1] %in%
            scientific_style_methods
        },
        logical(1)
      )
    )
  }
  
  
  # ==========================================================
  # F. Common-name / descriptive phrase recovery
  #
  # Only use this when we have not already found a plausible
  # scientific-looking taxon.
  #
  # Examples:
  #
  #   pearl millet stem
  #   Curry leaf plant
  #   Tomato plant roots
  #   Human ear
  # ==========================================================
  
  if (!currently_scientific) {
    
    descriptor_candidates <-
      .host_descriptor_candidates(
        term
      )
    
    for (candidate in descriptor_candidates) {
      
      add_candidate(
        candidate,
        "descriptor_cleanup",
        4
      )
    }
    
    
    # --------------------------------------------------------
    # Short contiguous phrases
    #
    # These are checked directly against the local NCBI names
    # table, including common names.
    # --------------------------------------------------------
    
    stop_words <- c(
      "a",
      "an",
      "the",
      "on",
      "in",
      "of",
      "from",
      "under",
      "near",
      "with",
      "and",
      "or",
      "cultivar",
      "cv",
      "sp",
      "spp"
    )
    
    
    for (segment in segments) {
      
      segment_clean <- .host_candidate_clean_text(
        segment
      )
      
      tokens <- unlist(
        strsplit(
          segment_clean,
          "[[:space:]]+"
        )
      )
      
      tokens <- tokens[
        nzchar(tokens)
      ]
      
      if (length(tokens) < 2) {
        next
      }
      
      maximum_n <- min(
        max_ngram,
        length(tokens)
      )
      
      for (n in seq.int(2, maximum_n)) {
        
        starts <- seq_len(
          length(tokens) - n + 1L
        )
        
        for (start in starts) {
          
          end <- start + n - 1L
          
          words <- tokens[
            start:end
          ]
          
          if (
            tolower(words[1]) %in% stop_words ||
            tolower(words[length(words)]) %in% stop_words
          ) {
            next
          }
          
          candidate <- paste(
            words,
            collapse = " "
          )
          
          add_candidate(
            candidate,
            "ncbi_phrase",
            5
          )
        }
      }
    }
  }
  
  
  # ==========================================================
  # No candidates
  # ==========================================================
  
  if (length(rows) == 0) {
    
    return(
      data.frame(
        original_term = character(0),
        candidate = character(0),
        candidate_method = character(0),
        method_priority = integer(0),
        stringsAsFactors = FALSE
      )
    )
  }
  
  
  out <- dplyr::bind_rows(
    rows
  )
  
  
  out$candidate_key <- .host_normalize_name(
    out$candidate
  )
  
  
  out <- out[
    !duplicated(
      paste(
        out$candidate_key,
        out$candidate_method,
        sep = "|||"
      )
    ),
    ,
    drop = FALSE
  ]
  
  
  out
}


.build_embedded_host_replacement_map <- function(
    terms,
    name_lookup_file = NULL,
    ranked_taxonomy_file = NULL
) {
  
  terms <- unique(
    .host_clean_term(
      terms
    )
  )
  
  terms <- terms[
    !.host_is_invalid_tax_term(
      terms
    )
  ]
  
  
  if (length(terms) == 0) {
    return(data.frame())
  }
  
  
  message(
    "Generating candidate taxa from ",
    format(
      length(terms),
      big.mark = ","
    ),
    " host term(s)..."
  )
  
  
  candidate_list <- lapply(
    terms,
    .generate_host_candidates_one
  )
  
  
  candidate_df <- dplyr::bind_rows(
    candidate_list
  )
  
  
  if (nrow(candidate_df) == 0) {
    return(data.frame())
  }
  
  
  message(
    "Generated ",
    format(
      nrow(candidate_df),
      big.mark = ","
    ),
    " candidate phrase(s)."
  )
  
  
  # ==========================================================
  # Load local NCBI names
  # ==========================================================
  
  name_lookup <- .read_ncbi_name_lookup(
    name_lookup_file
  )
  
  
  ranked_taxonomy <- .read_ncbi_ranked_taxonomy(
    ranked_taxonomy_file
  )
  
  
  candidate_keys <- unique(
    candidate_df$candidate_key
  )
  
  
  # Reduce the huge NCBI table BEFORE joining.
  relevant_names <- name_lookup[
    name_lookup$normalized_name %in%
      candidate_keys,
    c(
      "TaxID",
      "name",
      "name_class",
      "normalized_name"
    ),
    drop = FALSE
  ]
  
  
  if (nrow(relevant_names) > 0) {
    
    hits <- merge(
      candidate_df,
      relevant_names,
      by.x = "candidate_key",
      by.y = "normalized_name",
      all = FALSE,
      sort = FALSE
    )
    
  } else {
    
    hits <- data.frame()
  }
  
  
  # ==========================================================
  # Single-word candidates must be actual scientific names.
  #
  # This prevents ordinary words that happen to occur in the
  # NCBI name table from being treated as taxa.
  # ==========================================================
  
  if (nrow(hits) > 0) {
    
    single_methods <- c(
      "embedded_genus_sp",
      "hybrid_genus",
      "context_single_taxon",
      "bare_segment_taxon"
    )
    
    single_rows <- hits$candidate_method %in%
      single_methods
    
    scientific_name_rows <- tolower(
      hits$name_class
    ) == "scientific name"
    
    hits <- hits[
      !single_rows |
        scientific_name_rows,
      ,
      drop = FALSE
    ]
  }
  
  
  # ==========================================================
  # Canonical scientific name for each matched TaxID
  # ==========================================================
  
  canonical_lookup <- ranked_taxonomy[
    ,
    c(
      "TaxID",
      "Scientific.name"
    ),
    drop = FALSE
  ]
  
  canonical_lookup$TaxID <- as.character(
    canonical_lookup$TaxID
  )
  
  canonical_lookup <- canonical_lookup[
    !duplicated(
      canonical_lookup$TaxID
    ),
    ,
    drop = FALSE
  ]
  
  
  if (nrow(hits) > 0) {
    
    hits$TaxID <- as.character(
      hits$TaxID
    )
    
    hits <- dplyr::left_join(
      hits,
      canonical_lookup,
      by = "TaxID"
    )
  }
  
  
  # Split hits for faster term-wise processing
  hit_split <- if (nrow(hits) > 0) {
    split(
      hits,
      hits$original_term
    )
  } else {
    list()
  }
  
  
  # ==========================================================
  # Summarize each original host term
  # ==========================================================
  
  output_rows <- lapply(
    terms,
    function(term) {
      
      this_hits <- hit_split[[term]]
      
      
      # ------------------------------------------------------
      # Flags
      # ------------------------------------------------------
      
      uncertain <- grepl(
        "\\bor\\b|\\(\\s*\\?\\s*\\)",
        term,
        ignore.case = TRUE,
        perl = TRUE
      )
      
      likely_non_taxonomic <- grepl(
        paste0(
          "(?i)\\b(",
          paste(
            c(
              "bulk soil",
              "soil",
              "sediment",
              "decaying wood",
              "submerged decaying wood",
              "scorched earth",
              "sand",
              "clay",
              "meadow",
              "forest",
              "rhizosphere"
            ),
            collapse = "|"
          ),
          ")\\b"
        ),
        term,
        perl = TRUE
      )
      
      
      # ------------------------------------------------------
      # Nothing matched NCBI
      # ------------------------------------------------------
      
      if (
        is.null(this_hits) ||
        nrow(this_hits) == 0
      ) {
        
        status <- if (
          likely_non_taxonomic
        ) {
          "likely_non_taxonomic"
        } else {
          "no_ncbi_candidate"
        }
        
        return(
          data.frame(
            original_term = term,
            suggested_replacement = NA_character_,
            resolution_status = status,
            TaxID = NA_character_,
            chosen_candidate = NA_character_,
            candidate_method = NA_character_,
            matched_name = NA_character_,
            matched_name_class = NA_character_,
            all_detected_taxids = NA_character_,
            all_detected_candidates = NA_character_,
            stringsAsFactors = FALSE
          )
        )
      }
      
      
      unique_taxids <- unique(
        this_hits$TaxID[
          !is.na(this_hits$TaxID) &
            nzchar(this_hits$TaxID)
        ]
      )
      
      
      # ------------------------------------------------------
      # Multiple different taxa detected
      # ------------------------------------------------------
      
      if (length(unique_taxids) > 1) {
        
        return(
          data.frame(
            original_term = term,
            suggested_replacement = NA_character_,
            resolution_status = "multiple_taxa_detected",
            TaxID = NA_character_,
            chosen_candidate = NA_character_,
            candidate_method = NA_character_,
            matched_name = NA_character_,
            matched_name_class = NA_character_,
            all_detected_taxids = paste(
              unique_taxids,
              collapse = "; "
            ),
            all_detected_candidates = paste(
              unique(
                this_hits$candidate
              ),
              collapse = "; "
            ),
            stringsAsFactors = FALSE
          )
        )
      }
      
      
      # ------------------------------------------------------
      # Exactly one TaxID
      # ------------------------------------------------------
      
      taxid <- unique_taxids[1]
      
      
      # Prefer stronger candidate-generation methods.
      class_priority <- ifelse(
        tolower(
          this_hits$name_class
        ) == "scientific name",
        1L,
        2L
      )
      
      ordering <- order(
        this_hits$method_priority,
        class_priority,
        -nchar(
          this_hits$candidate
        )
      )
      
      chosen <- this_hits[
        ordering[1],
        ,
        drop = FALSE
      ]
      
      
      canonical_name <- chosen$Scientific.name[1]
      
      if (
        is.na(canonical_name) ||
        !nzchar(canonical_name)
      ) {
        canonical_name <- chosen$name[1]
      }
      
      
      status <- if (
        uncertain
      ) {
        "review_uncertain"
      } else {
        "auto_safe"
      }
      
      
      data.frame(
        original_term = term,
        suggested_replacement = canonical_name,
        resolution_status = status,
        TaxID = taxid,
        chosen_candidate = chosen$candidate[1],
        candidate_method = chosen$candidate_method[1],
        matched_name = chosen$name[1],
        matched_name_class = chosen$name_class[1],
        all_detected_taxids = taxid,
        all_detected_candidates = paste(
          unique(
            this_hits$candidate
          ),
          collapse = "; "
        ),
        stringsAsFactors = FALSE
      )
    }
  )
  
  
  dplyr::bind_rows(
    output_rows
  )
}


# ============================================================
# Run advanced extraction against the FAILED TERM TABLE
#
# This is the version to use first on the current fungal
# dataset.
#
# apply_suggestions = FALSE:
#   preview only
#
# apply_suggestions = TRUE:
#   fills replacement_term ONLY for auto_safe rows
#
# The large metadata file is NOT altered by this function.
# ============================================================

suggest_embedded_host_replacements <- function(
    project_name,
    host_dir = "./host_assessment",
    name_lookup_file = NULL,
    ranked_taxonomy_file = NULL,
    apply_suggestions = FALSE,
    overwrite_existing_replacements = FALSE
) {
  
  failed_df <- .read_failed_terms(
    project_name,
    host_dir
  )
  
  
  if (nrow(failed_df) == 0) {
    message(
      "No failed host terms found."
    )
    return(invisible(NULL))
  }
  
  
  has_replacement <- !.host_is_invalid_tax_term(
    failed_df$replacement_term
  )
  
  
  eligible <- if (
    overwrite_existing_replacements
  ) {
    rep(
      TRUE,
      nrow(failed_df)
    )
  } else {
    !has_replacement
  }
  
  
  eligible <- eligible &
    !.host_is_invalid_tax_term(
      failed_df$original_term
    )
  
  
  terms <- unique(
    failed_df$original_term[
      eligible
    ]
  )
  
  
  message(
    "Unresolved failed terms being examined: ",
    format(
      length(terms),
      big.mark = ","
    )
  )
  
  
  result <- .build_embedded_host_replacement_map(
    terms = terms,
    name_lookup_file = name_lookup_file,
    ranked_taxonomy_file = ranked_taxonomy_file
  )
  
  
  if (nrow(result) == 0) {
    message(
      "No candidate replacements were generated."
    )
    return(invisible(NULL))
  }
  
  
  # ----------------------------------------------------------
  # Add accession counts
  # ----------------------------------------------------------
  
  count_lookup <- failed_df[
    ,
    c(
      "original_term",
      "accession_count"
    ),
    drop = FALSE
  ]
  
  count_lookup <- count_lookup[
    !duplicated(
      count_lookup$original_term
    ),
    ,
    drop = FALSE
  ]
  
  
  result$accession_count <- suppressWarnings(
    as.numeric(
      count_lookup$accession_count[
        match(
          result$original_term,
          count_lookup$original_term
        )
      ]
    )
  )
  
  
  result$accession_count[
    is.na(result$accession_count)
  ] <- 0
  
  
  result <- result[
    order(
      factor(
        result$resolution_status,
        levels = c(
          "auto_safe",
          "review_uncertain",
          "multiple_taxa_detected",
          "likely_non_taxonomic",
          "no_ncbi_candidate"
        )
      ),
      -result$accession_count,
      result$original_term
    ),
    ,
    drop = FALSE
  ]
  
  
  preview_path <- file.path(
    host_dir,
    paste0(
      "host_embedded_name_suggestions_",
      project_name,
      ".csv"
    )
  )
  
  
  write.csv(
    result,
    preview_path,
    row.names = FALSE
  )
  
  
  message(
    "Embedded-name suggestion table written to: ",
    preview_path
  )
  
  
  message("")
  message("Summary:")
  
  status_table <- table(
    result$resolution_status
  )
  
  for (status in names(status_table)) {
    message(
      "  ",
      status,
      ": ",
      status_table[[status]]
    )
  }
  
  
  safe <- result[
    result$resolution_status ==
      "auto_safe",
    ,
    drop = FALSE
  ]
  
  
  message(
    "  Accessions represented by auto-safe replacements: ",
    format(
      sum(
        safe$accession_count,
        na.rm = TRUE
      ),
      big.mark = ","
    )
  )
  
  
  # ----------------------------------------------------------
  # Preview only
  # ----------------------------------------------------------
  
  if (!isTRUE(apply_suggestions)) {
    
    message("")
    message(
      "Preview only. No replacement terms were changed."
    )
    
    return(
      invisible(result)
    )
  }
  
  
  # ----------------------------------------------------------
  # Apply only auto-safe replacements
  # ----------------------------------------------------------
  
  if (nrow(safe) == 0) {
    
    message(
      "No auto-safe replacements are available to apply."
    )
    
    return(
      invisible(result)
    )
  }
  
  
  for (i in seq_len(nrow(safe))) {
    
    original <- safe$original_term[i]
    replacement <- safe$suggested_replacement[i]
    
    idx <- which(
      failed_df$original_term ==
        original
    )
    
    if (!overwrite_existing_replacements) {
      
      idx <- idx[
        .host_is_invalid_tax_term(
          failed_df$replacement_term[idx]
        )
      ]
    }
    
    if (length(idx) == 0) {
      next
    }
    
    
    failed_df$replacement_term[idx] <-
      replacement
    
    
    note <- paste0(
      "automatic embedded-name recovery: ",
      safe$candidate_method[i],
      " [",
      safe$chosen_candidate[i],
      "]"
    )
    
    
    for (j in idx) {
      
      if (
        is.na(failed_df$notes[j]) ||
        !nzchar(
          trimws(
            failed_df$notes[j]
          )
        )
      ) {
        
        failed_df$notes[j] <- note
        
      } else {
        
        failed_df$notes[j] <- paste(
          failed_df$notes[j],
          note,
          sep = "; "
        )
      }
    }
  }
  
  
  .write_failed_terms(
    failed_df,
    project_name,
    host_dir
  )
  
  
  message(
    "Auto-safe replacements were added to the failed-term table."
  )
  
  message(
    "The metadata table has not been modified."
  )
  
  
  invisible(result)
}


# ============================================================
# Phylogeny creation functions
# ============================================================


# Create Projects/<project_name>/ and run normal structure inside it
if (!exists("arborist_repo", envir = .GlobalEnv)) {
  arborist_repo <- normalizePath("~/github/aRborist")
}

start_project <- function(project_name,
                          projects_dir = file.path(arborist_repo, "projects")) {
  
  if (!dir.exists(projects_dir)) {
    dir.create(projects_dir, recursive = TRUE)
    message("Created projects dir: ", projects_dir)
  }
  
  base_dir <- file.path(projects_dir, project_name)
  if (!dir.exists(base_dir)) {
    dir.create(base_dir, recursive = TRUE)
    message("Created project dir: ", base_dir)
  }
  
  # make these visible to the rest of the pipeline
  assign("project_name", project_name, envir = .GlobalEnv)
  assign("base_dir", base_dir, envir = .GlobalEnv)
  
  # optional but useful: work inside the project
  setwd(base_dir)
  
  # your existing function that makes metadata_files/, etc.
  if (exists("setup_project_structure")) {
    setup_project_structure(base_dir)
  }
  
  invisible(normalizePath(base_dir))
}

# option to flag accessions from literature
flag_literature_accessions <- function(project_name,
                                       literature_accessions = NULL,
                                       curated_metadata_file = NULL) {
  if (is.null(literature_accessions) || literature_accessions == "") {
    message("No literature accession table supplied. Skipping.")
    return(invisible(NULL))
  }
  
  if (is.null(curated_metadata_file)) {
    curated_metadata_file <- file.path(
      "metadata_files",
      paste0("all_accessions_pulled_metadata_", project_name, "_curated.csv")
    )
  }
  
  metadata <- read.csv(curated_metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
  
  lit <- read.table(
    literature_accessions,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
  
  required_cols <- c("paper_id", "region", "accession")
  missing_cols <- setdiff(required_cols, names(lit))
  
  if (length(missing_cols) > 0) {
    stop(
      "Literature accession table is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  if (!"Accession" %in% names(metadata)) {
    stop("Metadata file must contain an 'Accession' column.")
  }
  
  lit$accession <- trimws(lit$accession)
  lit$region <- trimws(lit$region)
  lit$paper_id <- trimws(lit$paper_id)
  
  metadata$Accession <- trimws(metadata$Accession)
  
  if (!"literature_accession" %in% names(metadata)) {
    metadata$literature_accession <- FALSE
  }
  
  if (!"literature_source" %in% names(metadata)) {
    metadata$literature_source <- NA_character_
  }
  
  if (!"literature_region" %in% names(metadata)) {
    metadata$literature_region <- NA_character_
  }
  
  matched <- lit[lit$accession %in% metadata$Accession, , drop = FALSE]
  missing <- lit[!lit$accession %in% metadata$Accession, , drop = FALSE]
  
  for (acc in unique(matched$accession)) {
    hit_rows <- which(metadata$Accession == acc)
    lit_rows <- matched[matched$accession == acc, , drop = FALSE]
    
    sources <- paste(unique(lit_rows$paper_id), collapse = "; ")
    regions <- paste(unique(lit_rows$region), collapse = "; ")
    
    metadata$literature_accession[hit_rows] <- TRUE
    metadata$literature_source[hit_rows] <- sources
    metadata$literature_region[hit_rows] <- regions
  }
  
  write.csv(metadata, curated_metadata_file, row.names = FALSE, na = "")
  
  missing_file <- file.path(
    "metadata_files",
    paste0("missing_literature_accessions_", project_name, ".csv")
  )
  
  write.csv(missing, missing_file, row.names = FALSE, na = "")
  
  message("Literature accession flagging complete.")
  message("Updated curated metadata: ", curated_metadata_file)
  message("Missing literature accessions written to: ", missing_file)
  
  invisible(metadata)
}



# curating region data
curate_metadata_regions <- function(project_name,
                                    mapping_file = NULL,
                                    title_priority_regions = c("ITS", "LSU", "SSU")) {
  
  if (is.null(mapping_file) || !nzchar(mapping_file)) {
    arborist_root <- get0(
      "arborist_repo",
      envir = .GlobalEnv,
      ifnotfound = normalizePath("~/github/aRborist", mustWork = FALSE)
    )
    
    mapping_file <- file.path(
      arborist_root,
      "example_data",
      "region_replacement_patterns.csv"
    )
  }
  
  infile <- paste0(
    "./metadata_files/all_accessions_pulled_metadata_",
    project_name,
    "_curated.csv"
  )
  
  acc_df <- read.csv(
    infile,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  acc_df$gene.region.components <- NA_character_
  acc_df$product.region.components <- NA_character_
  acc_df$acc_title.region.components <- NA_character_
  
  if (file.exists(mapping_file)) {
    message("Using region replacement patterns from: ", mapping_file)
    map_df <- read.csv(mapping_file, stringsAsFactors = FALSE)
  } else {
    stop("Mapping file not found: ", mapping_file)
  }
  
  append_component <- function(current, add) {
    if (is.na(current) || current == "") {
      return(add)
    }
    
    current_parts <- trimws(unlist(strsplit(current, ";")))
    
    if (add %in% current_parts) {
      return(current)
    }
    
    paste(c(current_parts, add), collapse = ";")
  }
  
  make_safe_pattern <- function(pat) {
    pat <- trimws(pat)
    
    # If the user supplied explicit regex syntax, respect it.
    if (grepl("\\\\b|\\[|\\]|\\(|\\)|\\||\\+|\\*|\\?|\\{", pat)) {
      return(pat)
    }
    
    # Very short gene symbols need word boundaries.
    if (grepl("^[A-Za-z0-9]+$", pat) && nchar(pat) <= 5) {
      return(paste0("\\b", pat, "\\b"))
    }
    
    # Plain words/phrases should not match inside larger words.
    paste0("\\b", gsub(" +", "\\\\s+", pat), "\\b")
  }
  
  safe_detect <- function(x, pat) {
    if (is.null(x)) return(rep(FALSE, length(x)))
    
    pat_safe <- make_safe_pattern(pat)
    
    out <- stringr::str_detect(
      x,
      stringr::regex(pat_safe, ignore_case = TRUE)
    )
    
    out[is.na(out)] <- FALSE
    out
  }
  
  for (i in seq_len(nrow(map_df))) {
    pat <- map_df$pattern[i]
    std <- map_df$standard[i]
    
    hit_gene <- safe_detect(acc_df$gene, pat)
    
    if (any(hit_gene)) {
      acc_df$gene.region.components[hit_gene] <- mapply(
        append_component,
        acc_df$gene.region.components[hit_gene],
        std,
        USE.NAMES = FALSE
      )
    }
    
    hit_prod <- safe_detect(acc_df$product, pat)
    
    if (any(hit_prod)) {
      acc_df$product.region.components[hit_prod] <- mapply(
        append_component,
        acc_df$product.region.components[hit_prod],
        std,
        USE.NAMES = FALSE
      )
    }
    
    hit_title <- safe_detect(acc_df$accession_title, pat)
    
    if (any(hit_title)) {
      acc_df$acc_title.region.components[hit_title] <- mapply(
        append_component,
        acc_df$acc_title.region.components[hit_title],
        std,
        USE.NAMES = FALSE
      )
    }
  }
  
  detect_components <- function(txt) {
    if (is.null(txt) || is.na(txt) || txt == "") {
      return(NA_character_)
    }
    
    patterns <- list(
      ITS = stringr::regex("internal transcribed spacer|\\bITS\\b", ignore_case = TRUE),
      SSU = stringr::regex("\\b18S\\b|small subunit ribosomal", ignore_case = TRUE),
      LSU = stringr::regex("\\b28S\\b|\\b26S\\b|large subunit ribosomal", ignore_case = TRUE)
    )
    
    found <- character(0)
    
    for (nm in names(patterns)) {
      if (stringr::str_detect(txt, patterns[[nm]])) {
        found <- c(found, nm)
      }
    }
    
    found <- unique(found)
    
    if (length(found) == 0) {
      NA_character_
    } else {
      paste(found, collapse = ";")
    }
  }
  
  for (row_i in seq_len(nrow(acc_df))) {
    
    if (is.na(acc_df$gene.region.components[row_i]) ||
        acc_df$gene.region.components[row_i] == "") {
      comp <- detect_components(acc_df$gene[row_i])
      if (!is.na(comp)) acc_df$gene.region.components[row_i] <- comp
    }
    
    if (is.na(acc_df$product.region.components[row_i]) ||
        acc_df$product.region.components[row_i] == "") {
      comp <- detect_components(acc_df$product[row_i])
      if (!is.na(comp)) acc_df$product.region.components[row_i] <- comp
    }
    
    if (is.na(acc_df$acc_title.region.components[row_i]) ||
        acc_df$acc_title.region.components[row_i] == "") {
      comp <- detect_components(acc_df$accession_title[row_i])
      if (!is.na(comp)) acc_df$acc_title.region.components[row_i] <- comp
    }
  }
  
  filter_title_components <- function(x) {
    if (is.na(x) || x == "") return(NA_character_)
    
    parts <- trimws(unlist(strsplit(x, ";")))
    parts <- parts[parts %in% title_priority_regions]
    
    if (length(parts) == 0) {
      NA_character_
    } else {
      paste(unique(parts), collapse = ";")
    }
  }
  
  acc_df$acc_title.region.components <- vapply(
    acc_df$acc_title.region.components,
    filter_title_components,
    character(1)
  )
  
  acc_df <- acc_df %>%
    dplyr::mutate(
      region.standard = dplyr::coalesce(
        gene.region.components,
        product.region.components,
        acc_title.region.components
      )
    )
  
  acc_df$fasta.header <- paste0(">", acc_df$org_name, "_", acc_df$strain.standard)
  acc_df$fasta.header.type <- paste0(">", acc_df$org_name, "_", acc_df$strain.standard.type)
  
  unmatched_idx <- which(is.na(acc_df$region.standard) | acc_df$region.standard == "")
  
  if (length(unmatched_idx) > 0) {
    desired_cols <- c(
      "accession", "Accession",
      "gene", "product", "accession_title",
      "gene.region.components",
      "product.region.components",
      "acc_title.region.components",
      "org_name", "strain.standard"
    )
    
    cols_to_log <- intersect(desired_cols, colnames(acc_df))
    unmatched_df <- acc_df[unmatched_idx, cols_to_log, drop = FALSE]
    
    unmatched_file <- paste0(
      "./metadata_files/unmatched_regions_",
      project_name,
      ".csv"
    )
    
    write.csv(unmatched_df, unmatched_file, row.names = FALSE)
    
    message(
      "Some records did not match any region pattern. These were written to: ",
      unmatched_file
    )
  }
  
  outfile <- paste0(
    "./metadata_files/all_accessions_pulled_metadata_",
    project_name,
    "_curated.csv"
  )
  
  write.csv(acc_df, outfile, row.names = FALSE)
  
  cat("Wrote region-curated metadata to:", outfile, "\n")
}


# filtering metadata to only selected regions
select_regions <- function(project_name,
                           regions_to_include,
                           acc_to_exclude = character(0),
                           min_region_requirement = length(regions_to_include),
                           allow_compound_regions_for = c("ITS"),
                           prefer_literature_accessions = FALSE) {
  
  if (!exists("base_dir", envir = .GlobalEnv)) {
    stop("`base_dir` is not defined. Run start_project() first.")
  }
  
  has_region <- function(x, rg) {
    vapply(
      strsplit(as.character(x), ";"),
      function(parts) {
        rg %in% trimws(parts)
      },
      logical(1)
    )
  }
  
  region_is_exact <- function(x, rg) {
    trimws(as.character(x)) == rg
  }
  
  regions_to_include <- sort_regions(regions_to_include)
  region_set_name <- paste(regions_to_include, collapse = ".")
  
  input_path <- file.path(
    base_dir,
    "metadata_files",
    paste0("all_accessions_pulled_metadata_", project_name, "_curated.csv")
  )
  
  if (!file.exists(input_path)) {
    stop("Curated metadata not found at: ", input_path)
  }
  
  accession_list <- read.csv(
    input_path,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  phylo_dir <- file.path(base_dir, "phylogenies", region_set_name)
  if (!dir.exists(phylo_dir)) dir.create(phylo_dir, recursive = TRUE)
  
  if (!is.null(acc_to_exclude) &&
      length(acc_to_exclude) > 0 &&
      any(acc_to_exclude != "")) {
    accession_list <- accession_list[
      !accession_list$Accession %in% acc_to_exclude,
      ,
      drop = FALSE
    ]
  }
  
  expanded_list <- list()
  
  for (rg in regions_to_include) {
    
    if (rg %in% allow_compound_regions_for) {
      rg_keep <- has_region(accession_list$region.standard, rg)
    } else {
      rg_keep <- region_is_exact(accession_list$region.standard, rg)
    }
    
    rg_rows <- accession_list[rg_keep, , drop = FALSE]
    
    if (nrow(rg_rows) == 0) {
      next
    }
    
    rg_rows$region.standard.original <- rg_rows$region.standard
    rg_rows$region.standard <- rg
    
    expanded_list[[rg]] <- rg_rows
  }
  
  if (length(expanded_list) == 0) {
    stop(
      "No accessions matched the requested regions under the current compound-region policy.\n",
      "Requested regions: ", paste(regions_to_include, collapse = ", "), "\n",
      "allow_compound_regions_for: ", paste(allow_compound_regions_for, collapse = ", ")
    )
  }
  
  multifasta_prep_expanded <- dplyr::bind_rows(expanded_list)
  
  output_long <- file.path(
    phylo_dir,
    paste0("selected_accessions_metadata_", project_name, ".", region_set_name, ".csv")
  )
  
  write.csv(multifasta_prep_expanded, output_long, row.names = FALSE)
  
  cols_to_keep <- c(
    "strain.standard.type",
    "organism",
    "Accession",
    "region.standard",
    "literature_accession",
    "literature_source",
    "literature_region"
  )
  
  multifasta_prep_complete <- multifasta_prep_expanded[
    ,
    intersect(cols_to_keep, names(multifasta_prep_expanded)),
    drop = FALSE
  ]
  
  has_lit_col <- "literature_accession" %in% names(multifasta_prep_complete)
  
  if (has_lit_col) {
    multifasta_prep_complete$literature_accession[
      is.na(multifasta_prep_complete$literature_accession)
    ] <- FALSE
  }
  
  if (prefer_literature_accessions && has_lit_col) {
    multifasta_prep_complete <- multifasta_prep_complete[
      order(
        !multifasta_prep_complete$literature_accession,
        multifasta_prep_complete$strain.standard.type,
        multifasta_prep_complete$region.standard,
        multifasta_prep_complete$Accession
      ),
      ,
      drop = FALSE
    ]
  }
  
  multifasta_prep_select <- dplyr::distinct(
    multifasta_prep_complete,
    strain.standard.type,
    region.standard,
    .keep_all = TRUE
  )
  
  select_region_attendance <- tidyr::pivot_wider(
    multifasta_prep_select,
    id_cols = c("strain.standard.type", "organism"),
    names_from = "region.standard",
    values_from = "Accession"
  )
  
  if (has_lit_col) {
    
    lit_summary <- multifasta_prep_complete %>%
      dplyr::filter(literature_accession == TRUE) %>%
      dplyr::mutate(
        literature_accession_detail = paste0(region.standard, ":", Accession)
      ) %>%
      dplyr::group_by(strain.standard.type) %>%
      dplyr::summarise(
        literature_accessions_available = paste(
          unique(literature_accession_detail),
          collapse = "; "
        ),
        literature_sources_available = if ("literature_source" %in% names(.)) {
          paste(unique(na.omit(literature_source)), collapse = "; ")
        } else {
          NA_character_
        },
        .groups = "drop"
      )
    
    select_region_attendance <- select_region_attendance %>%
      dplyr::left_join(lit_summary, by = "strain.standard.type")
  }
  
  select_region_attendance_filtered <- select_region_attendance %>%
    dplyr::mutate(
      total = rowSums(
        !is.na(dplyr::select(., tidyselect::any_of(regions_to_include))) &
          dplyr::select(., tidyselect::any_of(regions_to_include)) != ""
      )
    ) %>%
    dplyr::filter(total >= min_region_requirement) %>%
    dplyr::select(-total)
  
  output_wide <- file.path(
    phylo_dir,
    paste0("Region_attendance_sheet_", project_name, ".", region_set_name, ".csv")
  )
  
  write.csv(select_region_attendance_filtered, output_wide, row.names = FALSE)
  
  policy_path <- file.path(
    phylo_dir,
    paste0("region_selection_policy_", project_name, ".", region_set_name, ".txt")
  )
  
  writeLines(
    c(
      paste0("project_name: ", project_name),
      paste0("region_set_name: ", region_set_name),
      paste0("regions_to_include: ", paste(regions_to_include, collapse = ", ")),
      paste0("min_region_requirement: ", min_region_requirement),
      paste0("allow_compound_regions_for: ", paste(allow_compound_regions_for, collapse = ", ")),
      paste0("prefer_literature_accessions: ", prefer_literature_accessions),
      "",
      "Rule:",
      "Regions listed in allow_compound_regions_for can be selected from compound region.standard values such as ITS;LSU;SSU.",
      "All other regions require exact region.standard matches.",
      "",
      "Literature accession rule:",
      "If prefer_literature_accessions = TRUE and a literature_accession column exists, literature accessions are prioritized when choosing one accession per strain/region.",
      "The attendance sheet includes literature_accessions_available and literature_sources_available columns when literature accession flags are present."
    ),
    con = policy_path
  )
  
  message("Filtered metadata written to: ", output_long)
  message("Region attendance sheet written to: ", output_wide)
  message("Region selection policy written to: ", policy_path)
  message("Region set: ", region_set_name)
  message("Regions included: ", paste(regions_to_include, collapse = ", "))
  message("Minimum region requirement: ", min_region_requirement)
  message("Compound regions allowed for: ", paste(allow_compound_regions_for, collapse = ", "))
  message("Prefer literature accessions: ", prefer_literature_accessions)
}


# optional filtering of strains in attendance sheet
filter_strains_for_tree <- function(project_name,
                                    attendance_file = NULL,
                                    metadata_file = NULL,
                                    output_file = NULL,
                                    strain_col = "strain.standard.type",
                                    include_col = "include_in_tree") {
  
  if (is.null(attendance_file)) {
    attendance_file <- file.path(
      "metadata_files",
      paste0("strain_attendance_sheet_", project_name, ".csv")
    )
  }
  
  if (is.null(metadata_file)) {
    metadata_file <- file.path(
      "metadata_files",
      paste0("all_accessions_pulled_metadata_", project_name, "_curated.csv")
    )
  }
  
  if (is.null(output_file)) {
    output_file <- file.path(
      "metadata_files",
      paste0("all_accessions_pulled_metadata_", project_name, "_curated_treefiltered.csv")
    )
  }
  
  attendance <- read.csv(attendance_file, stringsAsFactors = FALSE, check.names = FALSE)
  metadata <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
  
  if (!include_col %in% names(attendance)) {
    stop(
      "The attendance sheet does not contain column: ", include_col, "\n",
      "Add this column and mark strains to keep with TRUE, yes, keep, or 1."
    )
  }
  
  if (!strain_col %in% names(attendance)) {
    stop("Attendance sheet does not contain strain column: ", strain_col)
  }
  
  if (!strain_col %in% names(metadata)) {
    stop("Metadata file does not contain strain column: ", strain_col)
  }
  
  keep_values <- c("TRUE", "true", "T", "t", "yes", "YES", "Yes",
                   "keep", "KEEP", "Keep", "1")
  
  strains_to_keep <- attendance[[strain_col]][
    as.character(attendance[[include_col]]) %in% keep_values
  ]
  
  filtered_metadata <- metadata[metadata[[strain_col]] %in% strains_to_keep, ]
  
  write.csv(filtered_metadata, output_file, row.names = FALSE)
  
  message("Original metadata rows: ", nrow(metadata))
  message("Filtered metadata rows: ", nrow(filtered_metadata))
  message("Strains retained: ", length(unique(filtered_metadata[[strain_col]])))
  message("Filtered metadata written to: ", output_file)
  
  return(filtered_metadata)
}


# making subfolders for unique analyses within a single region set. For easier comparisons
start_phylogeny_run <- function(project_name,
                                regions_to_include,
                                run_label = NULL) {
  if (!exists("base_dir", envir = .GlobalEnv)) {
    stop("`base_dir` is not defined. Run start_project() first.")
  }
  
  regions_to_include <- sort_regions(regions_to_include)
  region_set_name <- paste(regions_to_include, collapse = ".")
  
  region_root_dir <- file.path(base_dir, "phylogenies", region_set_name)
  runs_dir <- file.path(region_root_dir, "runs")
  
  if (!dir.exists(region_root_dir)) {
    stop("Region-set folder not found: ", region_root_dir,
         "\nDid you run select_regions() first?")
  }
  
  if (!dir.exists(runs_dir)) dir.create(runs_dir, recursive = TRUE)
  
  existing_runs <- list.dirs(runs_dir, recursive = FALSE, full.names = FALSE)
  existing_nums <- suppressWarnings(as.integer(sub("^run_([0-9]+).*", "\\1", existing_runs)))
  existing_nums <- existing_nums[!is.na(existing_nums)]
  
  next_num <- if (length(existing_nums) == 0) 1 else max(existing_nums) + 1
  run_name <- sprintf("run_%03d", next_num)
  
  if (!is.null(run_label) && nzchar(run_label)) {
    clean_label <- gsub("[^A-Za-z0-9_-]+", "_", run_label)
    run_name <- paste0(run_name, "_", clean_label)
  }
  
  run_dir <- file.path(runs_dir, run_name)
  
  dir.create(file.path(run_dir, "prep"), recursive = TRUE)
  dir.create(file.path(run_dir, "single_gene_trees"), recursive = TRUE)
  dir.create(file.path(run_dir, "multi_gene_trees"), recursive = TRUE)
  dir.create(file.path(run_dir, "logs"), recursive = TRUE)
  
  message("Created phylogeny run folder: ", run_dir)
  return(run_dir)
}


get_phylo_paths <- function(project_name,
                            regions_to_include,
                            run_dir = NULL) {
  if (!exists("base_dir", envir = .GlobalEnv)) {
    stop("`base_dir` is not defined. Run start_project() first.")
  }
  
  regions_to_include <- sort_regions(regions_to_include)
  region_set_name <- paste(regions_to_include, collapse = ".")
  
  region_root_dir <- file.path(base_dir, "phylogenies", region_set_name)
  
  if (is.null(run_dir)) {
    analysis_dir <- region_root_dir
  } else {
    analysis_dir <- normalizePath(run_dir, mustWork = FALSE)
  }
  
  list(
    region_set_name = region_set_name,
    region_root_dir = region_root_dir,
    analysis_dir = analysis_dir,
    prep_dir = file.path(analysis_dir, "prep"),
    single_gene_dir = file.path(analysis_dir, "single_gene_trees"),
    multi_gene_dir = file.path(analysis_dir, "multi_gene_trees")
  )
}



# creating multifastas
create_multifastas <- function(project_name,
                               regions_to_include,
                               run_dir = NULL,
                               use_tree_filter = TRUE,
                               include_col = "include_in_tree",
                               strain_col = "strain.standard.type") {
  
  if (!exists("base_dir", envir = .GlobalEnv)) {
    stop("`base_dir` is not defined. Run start_project() first.")
  }
  
  paths <- get_phylo_paths(
    project_name = project_name,
    regions_to_include = regions_to_include,
    run_dir = run_dir
  )
  
  region_set_name <- paths$region_set_name
  
  # Shared dataset folder from select_regions()
  region_root_dir <- paths$region_root_dir
  
  # THIS run's output folder
  analysis_dir <- paths$analysis_dir
  prep_dir <- paths$prep_dir
  
  if (!dir.exists(region_root_dir)) {
    stop("Expected region-set folder not found: ", region_root_dir,
         "\nDid you run select_regions() for this region set?")
  }
  
  if (!dir.exists(prep_dir)) {
    dir.create(prep_dir, recursive = TRUE)
  }
  
  # One subfolder per region
  for (rg in regions_to_include) {
    rg_dir <- file.path(prep_dir, rg)
    if (!dir.exists(rg_dir)) dir.create(rg_dir, recursive = TRUE)
  }
  
  # ------------------------------------------------------------------
  # READ SHARED INPUT FILES FROM REGION ROOT
  # ------------------------------------------------------------------
  
  attendance_path <- file.path(
    region_root_dir,
    paste0(
      "Region_attendance_sheet_",
      project_name,
      ".",
      region_set_name,
      ".csv"
    )
  )
  
  long_filtered_path <- file.path(
    region_root_dir,
    paste0(
      "selected_accessions_metadata_",
      project_name,
      ".",
      region_set_name,
      ".csv"
    )
  )
  
  if (!file.exists(attendance_path)) {
    stop("Region attendance sheet not found: ", attendance_path)
  }
  
  if (!file.exists(long_filtered_path)) {
    stop("Filtered metadata (long) not found: ", long_filtered_path)
  }
  
  region_attendance <- read.csv(
    attendance_path,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  filtered_long <- read.csv(
    long_filtered_path,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # ------------------------------------------------------------------
  # OPTIONAL TREE FILTERING
  # ------------------------------------------------------------------
  
  if (isTRUE(use_tree_filter)) {
    
    if (include_col %in% colnames(region_attendance)) {
      
      if (!strain_col %in% colnames(region_attendance)) {
        stop(
          "Tree filter column found, but strain column is missing from attendance sheet: ",
          strain_col
        )
      }
      
      if (!strain_col %in% colnames(filtered_long)) {
        stop(
          "Tree filter column found, but strain column is missing from selected metadata: ",
          strain_col
        )
      }
      
      keep_values <- c(
        "TRUE", "true", "True",
        "T", "t",
        "yes", "YES", "Yes",
        "keep", "KEEP", "Keep",
        "1"
      )
      
      keep_rows <- as.character(region_attendance[[include_col]]) %in% keep_values
      
      strains_to_keep <- unique(region_attendance[[strain_col]][keep_rows])
      
      strains_to_keep <- strains_to_keep[
        !is.na(strains_to_keep) &
          strains_to_keep != ""
      ]
      
      message("Tree filter detected: ", include_col)
      message("Strains marked for inclusion: ", length(strains_to_keep))
      
      original_attendance_n <- nrow(region_attendance)
      original_long_n <- nrow(filtered_long)
      
      region_attendance <- region_attendance[
        region_attendance[[strain_col]] %in% strains_to_keep,
        ,
        drop = FALSE
      ]
      
      filtered_long <- filtered_long[
        filtered_long[[strain_col]] %in% strains_to_keep,
        ,
        drop = FALSE
      ]
      
      message(
        "Attendance rows retained: ",
        nrow(region_attendance),
        " / ",
        original_attendance_n
      )
      
      message(
        "Metadata rows retained: ",
        nrow(filtered_long),
        " / ",
        original_long_n
      )
      
      # Write filtered copies INSIDE RUN FOLDER
      filtered_attendance_path <- file.path(
        analysis_dir,
        paste0(
          "Region_attendance_sheet_",
          project_name,
          ".",
          region_set_name,
          "_treefiltered.csv"
        )
      )
      
      filtered_long_path <- file.path(
        analysis_dir,
        paste0(
          "selected_accessions_metadata_",
          project_name,
          ".",
          region_set_name,
          "_treefiltered.csv"
        )
      )
      
      write.csv(region_attendance, filtered_attendance_path, row.names = FALSE)
      write.csv(filtered_long, filtered_long_path, row.names = FALSE)
      
      message("Wrote tree-filtered attendance sheet: ", filtered_attendance_path)
      message("Wrote tree-filtered selected metadata: ", filtered_long_path)
      
    } else {
      
      stop(
        "Tree filtering requested, but no column named '",
        include_col,
        "' was found in:\n",
        attendance_path,
        "\n\nAdd an include_in_tree column to this exact file, then rerun create_multifastas()."
      )
    }
  }
  
  # ------------------------------------------------------------------
  # SANITY CHECKS
  # ------------------------------------------------------------------
  
  needed_cols <- c(
    "Accession",
    "region.standard",
    "fasta.header.type",
    "sequence"
  )
  
  missing_cols <- setdiff(needed_cols, colnames(filtered_long))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing columns in filtered metadata: ",
      paste(missing_cols, collapse = ", "),
      "\nUpstream curation must provide these."
    )
  }
  
  # ------------------------------------------------------------------
  # BUILD ACCESSION VECTORS
  # ------------------------------------------------------------------
  
  region_cols <- intersect(
    regions_to_include,
    colnames(region_attendance)
  )
  
  if (length(region_cols) == 0) {
    stop(
      "None of the requested regions are present as columns in the attendance sheet."
    )
  }
  
  region_accessions <- lapply(region_cols, function(rg) {
    unique(na.omit(region_attendance[[rg]]))
  })
  
  names(region_accessions) <- region_cols
  
  # ------------------------------------------------------------------
  # WRITE RAW FASTAS
  # ------------------------------------------------------------------
  
  manifest <- data.frame(
    region = character(0),
    n_sequences = integer(0),
    fasta_path = character(0),
    stringsAsFactors = FALSE
  )
  
  for (rg in names(region_accessions)) {
    
    acc_vec <- region_accessions[[rg]]
    
    if (length(acc_vec) == 0) {
      message("No accessions found for region: ", rg, " (skipping).")
      next
    }
    
    sub_df <- filtered_long[
      filtered_long$Accession %in% acc_vec &
        filtered_long$region.standard == rg,
    ]
    
    sub_df <- sub_df[
      !is.na(sub_df$sequence) &
        sub_df$sequence != "",
    ]
    
    sub_df <- sub_df[
      order(sub_df$fasta.header.type, decreasing = FALSE),
    ]
    
    headers <- sub_df$fasta.header.type
    
    needs_gt <- !startsWith(headers, ">")
    headers[needs_gt] <- paste0(">", headers[needs_gt])
    
    seqs_fasta <- c(rbind(headers, sub_df$sequence))
    
    rg_dir <- file.path(prep_dir, rg)
    
    fasta_name <- paste0(
      project_name,
      ".",
      region_set_name,
      "_",
      rg,
      ".raw.fasta"
    )
    
    fasta_path <- file.path(rg_dir, fasta_name)
    
    writeLines(seqs_fasta, con = fasta_path)
    
    message("Created multifasta for region ", rg, ": ", fasta_path)
    
    manifest <- rbind(
      manifest,
      data.frame(
        region = rg,
        n_sequences = nrow(sub_df),
        fasta_path = fasta_path,
        stringsAsFactors = FALSE
      )
    )
  }
  
  # ------------------------------------------------------------------
  # MANIFEST
  # ------------------------------------------------------------------
  
  manifest_path <- file.path(
    prep_dir,
    paste0(
      "multifasta_manifest_",
      project_name,
      ".",
      region_set_name,
      ".tsv"
    )
  )
  
  write.table(
    manifest,
    manifest_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  message("Wrote manifest: ", manifest_path)
  
  invisible(manifest)
}




align_regions_mafft <- function(project_name,
                                regions_to_include,
                                run_dir = NULL,
                                threads = max(1, parallel::detectCores() - 1),
                                mafft_args = c("--auto", "--reorder"),
                                force = FALSE) {
  
  if (!exists("base_dir", envir = .GlobalEnv)) {
    stop("`base_dir` is not defined. Run start_project() first.")
  }
  
  paths <- get_phylo_paths(
    project_name = project_name,
    regions_to_include = regions_to_include,
    run_dir = run_dir
  )
  
  regions_to_include <- sort_regions(regions_to_include)
  region_set_name <- paths$region_set_name
  prep_dir <- paths$prep_dir
  
  mafft_path <- Sys.getenv("MAFFT_PATH", unset = "mafft")
  
  check_result <- suppressWarnings(
    system2(mafft_path, "--version", stdout = TRUE, stderr = TRUE)
  )
  
  if (
    length(check_result) == 0 ||
    grepl("not found|No such file", check_result[1], ignore.case = TRUE)
  ) {
    stop(
      "MAFFT not found. Please install it or set MAFFT_PATH in your .Renviron file.\n",
      "Example:  MAFFT_PATH=/usr/local/bin/mafft\n",
      "Then restart R and rerun this command."
    )
  } else {
    message("Using MAFFT executable: ", mafft_path)
  }
  
  if (!dir.exists(prep_dir)) {
    stop(
      "Prep directory not found: ",
      prep_dir,
      "\nDid you run create_multifastas() for this run?"
    )
  }
  
  manifest <- data.frame(
    region = character(0),
    raw_fasta = character(0),
    aligned_fasta = character(0),
    log_path = character(0),
    status = character(0),
    stringsAsFactors = FALSE
  )
  
  for (rg in regions_to_include) {
    
    rg_dir <- file.path(prep_dir, rg)
    
    if (!dir.exists(rg_dir)) {
      warning("Region prep folder missing (skipping): ", rg_dir)
      next
    }
    
    raw_fa <- file.path(
      rg_dir,
      paste0(project_name, ".", region_set_name, "_", rg, ".raw.fasta")
    )
    
    aln_fa <- file.path(
      rg_dir,
      paste0(project_name, ".", region_set_name, "_", rg, ".aligned.fasta")
    )
    
    log_fp <- file.path(
      rg_dir,
      paste0(project_name, ".", region_set_name, "_", rg, ".mafft.log")
    )
    
    if (!file.exists(raw_fa)) {
      warning("Raw FASTA not found for region ", rg, ": ", raw_fa)
      
      manifest <- rbind(
        manifest,
        data.frame(
          region = rg,
          raw_fasta = raw_fa,
          aligned_fasta = NA,
          log_path = log_fp,
          status = "missing_raw",
          stringsAsFactors = FALSE
        )
      )
      
      next
    }
    
    if (file.exists(aln_fa) && !force) {
      message("Aligned FASTA already exists; use force=TRUE to overwrite: ", aln_fa)
      
      manifest <- rbind(
        manifest,
        data.frame(
          region = rg,
          raw_fasta = raw_fa,
          aligned_fasta = aln_fa,
          log_path = log_fp,
          status = "skipped_exists",
          stringsAsFactors = FALSE
        )
      )
      
      next
    }
    
    message("Running MAFFT for region ", rg, " ...")
    
    mafft_args_full <- c(
      "--thread",
      as.character(threads),
      mafft_args,
      raw_fa
    )
    
    exit_code <- tryCatch(
      {
        system2(
          command = mafft_path,
          args = mafft_args_full,
          stdout = aln_fa,
          stderr = log_fp
        )
      },
      error = function(e) {
        warning("MAFFT invocation failed for ", rg, ": ", conditionMessage(e))
        return(1L)
      }
    )
    
    status <- if (
      !is.null(exit_code) &&
      exit_code == 0L &&
      file.exists(aln_fa)
    ) {
      "ok"
    } else {
      "failed"
    }
    
    manifest <- rbind(
      manifest,
      data.frame(
        region = rg,
        raw_fasta = raw_fa,
        aligned_fasta = if (file.exists(aln_fa)) aln_fa else NA,
        log_path = log_fp,
        status = status,
        stringsAsFactors = FALSE
      )
    )
    
    if (status != "ok") {
      warning("MAFFT failed for region ", rg, ". See log: ", log_fp)
    } else {
      message("Aligned FASTA written: ", aln_fa)
    }
  }
  
  align_manifest <- file.path(
    prep_dir,
    paste0("alignment_manifest_", project_name, ".", region_set_name, ".tsv")
  )
  
  write.table(
    manifest,
    align_manifest,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  message("Alignment manifest: ", align_manifest)
  
  invisible(manifest)
}



trim_regions_trimal <- function(project_name,
                                regions_to_include,
                                run_dir = NULL,
                                trimal_args = c("-automated1"),
                                force = FALSE) {
  
  if (!exists("base_dir", envir = .GlobalEnv)) {
    stop("`base_dir` is not defined. Run start_project() first.")
  }
  
  paths <- get_phylo_paths(
    project_name = project_name,
    regions_to_include = regions_to_include,
    run_dir = run_dir
  )
  
  regions_to_include <- sort_regions(regions_to_include)
  region_set_name <- paths$region_set_name
  prep_dir <- paths$prep_dir
  
  trimal_path <- Sys.getenv("TRIMAL_PATH", unset = "trimal")
  
  check_result <- suppressWarnings(
    system2(trimal_path, "--version", stdout = TRUE, stderr = TRUE)
  )
  
  if (
    length(check_result) == 0 ||
    grepl("not found|No such file", check_result[1], ignore.case = TRUE)
  ) {
    stop(
      "trimAl not found. Please install it or set TRIMAL_PATH in your .Renviron file.\n",
      "Example:  TRIMAL_PATH=/usr/local/bin/trimal\n",
      "Then restart R and rerun this command."
    )
  } else {
    message("Using trimAl executable: ", trimal_path)
  }
  
  if (!dir.exists(prep_dir)) {
    stop(
      "Prep directory not found: ",
      prep_dir,
      "\nDid you run align_regions_mafft() for this run?"
    )
  }
  
  manifest <- data.frame(
    region = character(0),
    aligned_fasta = character(0),
    trimmed_fasta = character(0),
    log_path = character(0),
    status = character(0),
    stringsAsFactors = FALSE
  )
  
  for (rg in regions_to_include) {
    
    rg_dir <- file.path(prep_dir, rg)
    
    aln_fa <- file.path(
      rg_dir,
      paste0(project_name, ".", region_set_name, "_", rg, ".aligned.fasta")
    )
    
    trimmed_fa <- file.path(
      rg_dir,
      paste0(project_name, ".", region_set_name, "_", rg, ".trimmed.fasta")
    )
    
    log_fp <- file.path(
      rg_dir,
      paste0(project_name, ".", region_set_name, "_", rg, ".trimal.log")
    )
    
    if (!file.exists(aln_fa)) {
      warning("Aligned FASTA not found for region ", rg, ": ", aln_fa)
      
      manifest <- rbind(
        manifest,
        data.frame(
          region = rg,
          aligned_fasta = aln_fa,
          trimmed_fasta = NA,
          log_path = log_fp,
          status = "missing_aligned",
          stringsAsFactors = FALSE
        )
      )
      
      next
    }
    
    if (file.exists(trimmed_fa) && !force) {
      message("Trimmed FASTA already exists; use force=TRUE to overwrite: ", trimmed_fa)
      
      manifest <- rbind(
        manifest,
        data.frame(
          region = rg,
          aligned_fasta = aln_fa,
          trimmed_fasta = trimmed_fa,
          log_path = log_fp,
          status = "skipped_exists",
          stringsAsFactors = FALSE
        )
      )
      
      next
    }
    
    message("Running trimAl for region ", rg, " ...")
    
    trimal_args_full <- c(
      trimal_args,
      "-in",
      aln_fa,
      "-out",
      trimmed_fa
    )
    
    exit_code <- tryCatch(
      {
        system2(
          command = trimal_path,
          args = trimal_args_full,
          stdout = log_fp,
          stderr = log_fp
        )
      },
      error = function(e) {
        warning("trimAl invocation failed for ", rg, ": ", conditionMessage(e))
        return(1L)
      }
    )
    
    status <- if (
      !is.null(exit_code) &&
      exit_code == 0L &&
      file.exists(trimmed_fa)
    ) {
      "ok"
    } else {
      "failed"
    }
    
    manifest <- rbind(
      manifest,
      data.frame(
        region = rg,
        aligned_fasta = aln_fa,
        trimmed_fasta = if (file.exists(trimmed_fa)) trimmed_fa else NA,
        log_path = log_fp,
        status = status,
        stringsAsFactors = FALSE
      )
    )
    
    if (status != "ok") {
      warning("trimAl failed for region ", rg, ". See log: ", log_fp)
    } else {
      message("Trimmed FASTA written: ", trimmed_fa)
    }
  }
  
  trim_manifest <- file.path(
    prep_dir,
    paste0("trim_manifest_", project_name, ".", region_set_name, ".tsv")
  )
  
  write.table(
    manifest,
    trim_manifest,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  message("Trim manifest: ", trim_manifest)
  
  invisible(manifest)
}


# write FINAL attendance sheet - only includes accessions that made it past the trimming process
# (some accessions/sequences may be removed automatically by trimal, if certain parameters are used)
write_final_region_attendance_sheet <- function(project_name,
                                                regions_to_include,
                                                run_dir = NULL,
                                                strain_col = "strain.standard.type") {
  
  if (!exists("base_dir", envir = .GlobalEnv)) {
    stop("`base_dir` is not defined. Run start_project() first.")
  }
  
  if (is.null(run_dir)) {
    stop("write_final_region_attendance_sheet() requires a run_dir.")
  }
  
  paths <- get_phylo_paths(
    project_name = project_name,
    regions_to_include = regions_to_include,
    run_dir = run_dir
  )
  
  regions_to_include <- sort_regions(regions_to_include)
  region_set_name <- paths$region_set_name
  region_root_dir <- paths$region_root_dir
  analysis_dir <- paths$analysis_dir
  prep_dir <- paths$prep_dir
  
  if (!dir.exists(prep_dir)) {
    stop("Prep directory not found: ", prep_dir)
  }
  
  intended_attendance_path <- file.path(
    analysis_dir,
    paste0(
      "Region_attendance_sheet_",
      project_name,
      ".",
      region_set_name,
      "_treefiltered.csv"
    )
  )
  
  if (!file.exists(intended_attendance_path)) {
    intended_attendance_path <- file.path(
      region_root_dir,
      paste0(
        "Region_attendance_sheet_",
        project_name,
        ".",
        region_set_name,
        ".csv"
      )
    )
  }
  
  if (!file.exists(intended_attendance_path)) {
    stop("Could not find intended/input attendance sheet.")
  }
  
  metadata_path <- file.path(
    analysis_dir,
    paste0(
      "selected_accessions_metadata_",
      project_name,
      ".",
      region_set_name,
      "_treefiltered.csv"
    )
  )
  
  if (!file.exists(metadata_path)) {
    metadata_path <- file.path(
      region_root_dir,
      paste0(
        "selected_accessions_metadata_",
        project_name,
        ".",
        region_set_name,
        ".csv"
      )
    )
  }
  
  if (!file.exists(metadata_path)) {
    stop("Could not find selected accessions metadata.")
  }
  
  intended_attendance <- read.csv(
    intended_attendance_path,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  selected_metadata <- read.csv(
    metadata_path,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  required_meta_cols <- c("Accession", "region.standard", "fasta.header.type")
  missing_meta_cols <- setdiff(required_meta_cols, names(selected_metadata))
  
  if (length(missing_meta_cols) > 0) {
    stop(
      "Selected metadata is missing required columns: ",
      paste(missing_meta_cols, collapse = ", ")
    )
  }
  
  if (!strain_col %in% names(intended_attendance)) {
    stop("Strain column not found in attendance sheet: ", strain_col)
  }
  
  final_attendance <- intended_attendance
  
  retention_report <- data.frame(
    region = character(0),
    accession = character(0),
    strain = character(0),
    fasta_header = character(0),
    status = character(0),
    stringsAsFactors = FALSE
  )
  
  selected_metadata$fasta.header.type <- sub(
    "^>",
    "",
    selected_metadata$fasta.header.type
  )
  
  for (rg in regions_to_include) {
    
    if (!rg %in% names(final_attendance)) {
      warning("Region column not found in attendance sheet: ", rg)
      next
    }
    
    rg_dir <- file.path(prep_dir, rg)
    
    aligned_fa <- file.path(
      rg_dir,
      paste0(project_name, ".", region_set_name, "_", rg, ".aligned.fasta")
    )
    
    trimmed_fa <- file.path(
      rg_dir,
      paste0(project_name, ".", region_set_name, "_", rg, ".trimmed.fasta")
    )
    
    if (!file.exists(aligned_fa)) {
      warning("Aligned FASTA not found for region ", rg, ": ", aligned_fa)
      next
    }
    
    if (!file.exists(trimmed_fa)) {
      warning("Trimmed FASTA not found for region ", rg, ": ", trimmed_fa)
      next
    }
    
    aligned_names <- names(Biostrings::readDNAStringSet(aligned_fa))
    trimmed_names <- names(Biostrings::readDNAStringSet(trimmed_fa))
    
    region_metadata <- selected_metadata[
      selected_metadata$region.standard == rg,
      ,
      drop = FALSE
    ]
    
    header_lookup <- region_metadata[, c("Accession", "fasta.header.type")]
    
    intended_acc <- final_attendance[[rg]]
    intended_acc_clean <- intended_acc[
      !is.na(intended_acc) &
        intended_acc != ""
    ]
    
    retained_acc <- character(0)
    
    for (acc in intended_acc_clean) {
      
      possible_headers <- header_lookup$fasta.header.type[
        header_lookup$Accession == acc
      ]
      
      possible_headers <- possible_headers[
        !is.na(possible_headers) &
          possible_headers != ""
      ]
      
      in_aligned <- any(possible_headers %in% aligned_names)
      in_trimmed <- any(possible_headers %in% trimmed_names)
      
      if (in_trimmed) {
        retained_acc <- c(retained_acc, acc)
      }
      
      strain_value <- final_attendance[[strain_col]][
        which(final_attendance[[rg]] == acc)[1]
      ]
      
      status <- if (in_trimmed) {
        "retained"
      } else if (in_aligned) {
        "removed_by_trimal"
      } else if (length(possible_headers) == 0) {
        "no_header_found_in_metadata"
      } else {
        "not_found_in_aligned_fasta"
      }
      
      retention_report <- rbind(
        retention_report,
        data.frame(
          region = rg,
          accession = acc,
          strain = strain_value,
          fasta_header = paste(possible_headers, collapse = ";"),
          status = status,
          stringsAsFactors = FALSE
        )
      )
    }
    
    removed_acc <- setdiff(intended_acc_clean, retained_acc)
    
    final_attendance[[rg]][
      final_attendance[[rg]] %in% removed_acc
    ] <- NA
  }
  
  region_cols <- intersect(regions_to_include, names(final_attendance))
  
  keep_rows <- rowSums(
    !is.na(final_attendance[, region_cols, drop = FALSE]) &
      final_attendance[, region_cols, drop = FALSE] != ""
  ) > 0
  
  final_attendance <- final_attendance[keep_rows, , drop = FALSE]
  
  intended_copy_path <- file.path(
    analysis_dir,
    paste0(
      "intended_region_attendance_sheet_",
      project_name,
      ".",
      region_set_name,
      ".csv"
    )
  )
  
  final_attendance_path <- file.path(
    analysis_dir,
    paste0(
      "Region_attendance_sheet_",
      project_name,
      ".",
      region_set_name,
      ".csv"
    )
  )
  
  retention_report_path <- file.path(
    analysis_dir,
    paste0(
      "trim_retention_report_",
      project_name,
      ".",
      region_set_name,
      ".tsv"
    )
  )
  
  write.csv(
    intended_attendance,
    intended_copy_path,
    row.names = FALSE
  )
  
  write.csv(
    final_attendance,
    final_attendance_path,
    row.names = FALSE
  )
  
  write.table(
    retention_report,
    retention_report_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  message("Wrote intended attendance sheet: ", intended_copy_path)
  message("Wrote final region attendance sheet: ", final_attendance_path)
  message("Wrote trim retention report: ", retention_report_path)
  
  invisible(
    list(
      intended_attendance = intended_attendance,
      final_attendance = final_attendance,
      retention_report = retention_report,
      intended_attendance_path = intended_copy_path,
      final_attendance_path = final_attendance_path,
      retention_report_path = retention_report_path
    )
  )
}



# running IQTREE modelfinder step for single region 
iqtree_modelfinder_per_region <- function(project_name,
                                          regions_to_include,
                                          run_dir = NULL,
                                          threads = max(1, parallel::detectCores() - 1),
                                          iqtree_args = c("-m", "MFP+MERGE", "-nt", "AUTO", "-quiet"),
                                          single_gene_bootstraps = 1000,
                                          force = FALSE) {
  
  if (!exists("base_dir", envir = .GlobalEnv)) {
    stop("`base_dir` is not defined. Run start_project() first.")
  }
  
  iqtree_bin <- Sys.getenv("IQTREE_PATH")
  
  if (!nzchar(iqtree_bin)) {
    stop("IQTREE_PATH is not set in .Renviron.")
  }
  
  suppressWarnings(system2(iqtree_bin, "-version"))
  
  paths <- get_phylo_paths(
    project_name = project_name,
    regions_to_include = regions_to_include,
    run_dir = run_dir
  )
  
  regions_to_include <- sort_regions(regions_to_include)
  region_set_name <- paths$region_set_name
  prep_dir <- paths$prep_dir
  single_gene_dir <- paths$single_gene_dir
  multi_gene_dir <- paths$multi_gene_dir
  
  if (!dir.exists(prep_dir)) {
    stop(
      "Prep directory not found: ",
      prep_dir,
      "\nDid you run trim_regions_trimal() for this run?"
    )
  }
  
  if (!dir.exists(single_gene_dir)) {
    dir.create(single_gene_dir, recursive = TRUE)
  }
  
  if (!dir.exists(multi_gene_dir)) {
    dir.create(multi_gene_dir, recursive = TRUE)
  }
  
  results <- list()
  
  for (rg in regions_to_include) {
    
    rg_prep_dir <- file.path(prep_dir, rg)
    
    trimmed_fa <- file.path(
      rg_prep_dir,
      paste0(project_name, ".", region_set_name, "_", rg, ".trimmed.fasta")
    )
    
    if (!file.exists(trimmed_fa)) {
      warning("Missing trimmed alignment for region ", rg, ": ", trimmed_fa)
      next
    }
    
    rg_sg_dir <- file.path(single_gene_dir, rg)
    
    if (!dir.exists(rg_sg_dir)) {
      dir.create(rg_sg_dir, recursive = TRUE)
    }
    
    prefix_base <- paste0(
      project_name,
      ".",
      region_set_name,
      "_",
      rg,
      ".modeltest"
    )
    
    prefix <- file.path(rg_sg_dir, prefix_base)
    iqtreefile <- paste0(prefix, ".iqtree")
    
    if (!file.exists(iqtreefile) || force) {
      
      args <- c(
        "-s",
        trimmed_fa,
        "-pre",
        prefix,
        "-bb",
        as.character(single_gene_bootstraps)
      )
      
      if (!any(iqtree_args == "-nt")) {
        args <- c(args, "-nt", as.character(threads))
      }
      
      args <- c(args, iqtree_args)
      
      message("Running ModelFinder with UF bootstraps for region ", rg, " ...")
      message(iqtree_bin, " ", paste(shQuote(args), collapse = " "))
      
      iqtree_log <- paste0(prefix, ".run.log")
      
      exit_status <- system2(
        command = iqtree_bin,
        args = args,
        stdout = iqtree_log,
        stderr = iqtree_log
      )
      
      if (!identical(exit_status, 0L)) {
        warning(
          "IQ-TREE ModelFinder finished with non-zero exit status for region ",
          rg,
          ": ",
          exit_status,
          "\nSee log: ",
          iqtree_log
        )
      }
    } else {
      message("IQ-TREE result already exists; use force=TRUE to overwrite: ", iqtreefile)
    }
    
    if (!file.exists(iqtreefile)) {
      warning("Expected IQ-TREE output not found for region ", rg, ": ", iqtreefile)
      next
    }
    
    iqtxt <- readLines(iqtreefile, warn = FALSE)
    
    model_idx <- grep(
      "Best-fit model according to BIC:",
      iqtxt,
      fixed = TRUE
    )
    
    if (length(model_idx) == 0L) {
      best_model <- NA_character_
      warning("Could not find model line in ", iqtreefile)
    } else {
      best_line <- iqtxt[model_idx[1]]
      
      best_model <- stringr::str_trim(
        sub(
          ".*Best-fit model according to BIC:\\s*",
          "",
          best_line
        )
      )
    }
    
    aln <- Biostrings::readDNAStringSet(trimmed_fa)
    aln_length <- unique(Biostrings::width(aln))[1]
    
    results[[rg]] <- data.frame(
      region = rg,
      best_model = best_model,
      aln_length = aln_length,
      iqtree_file = iqtreefile,
      stringsAsFactors = FALSE
    )
  }
  
  model_fits <- dplyr::bind_rows(results)
  
  if (nrow(model_fits) == 0 || !"region" %in% names(model_fits)) {
    stop(
      "No successful IQ-TREE ModelFinder results were generated.\n",
      "This usually means the trimmed FASTA files were missing, IQ-TREE failed, ",
      "or the expected .iqtree files were not created.\n\n",
      "Check this folder:\n",
      single_gene_dir
    )
  }
  
  out_tsv <- file.path(
    multi_gene_dir,
    paste0("model_fits_", project_name, ".", region_set_name, ".tsv")
  )
  
  readr::write_tsv(model_fits, out_tsv)
  
  message("ModelFinder summary written to: ", out_tsv)
  
  invisible(model_fits)
}


# create input files for iqtree concatenated analysis
concatenate_and_write_partitions <- function(project_name,
                                             regions_to_include,
                                             run_dir = NULL) {
  
  if (!exists("base_dir", envir = .GlobalEnv)) {
    stop("`base_dir` is not defined. Run start_project() first.")
  }
  
  paths <- get_phylo_paths(
    project_name = project_name,
    regions_to_include = regions_to_include,
    run_dir = run_dir
  )
  
  regions_to_include <- sort_regions(regions_to_include)
  region_set_name <- paths$region_set_name
  prep_dir <- paths$prep_dir
  multi_gene_dir <- paths$multi_gene_dir
  
  if (!dir.exists(prep_dir)) {
    stop(
      "Prep directory not found: ",
      prep_dir,
      "\nDid you run trim_regions_trimal() for this run?"
    )
  }
  
  if (!dir.exists(multi_gene_dir)) {
    dir.create(multi_gene_dir, recursive = TRUE)
  }
  
  model_fits_path <- file.path(
    multi_gene_dir,
    paste0("model_fits_", project_name, ".", region_set_name, ".tsv")
  )
  
  if (!file.exists(model_fits_path)) {
    stop(
      "model_fits TSV not found: ",
      model_fits_path,
      "\nDid you run iqtree_modelfinder_per_region() for this run?"
    )
  }
  
  model_fits <- readr::read_tsv(
    model_fits_path,
    show_col_types = FALSE
  )
  
  if (!all(c("region", "best_model") %in% names(model_fits))) {
    stop(
      "model_fits TSV does not contain the required columns: region, best_model\n",
      "File checked: ",
      model_fits_path
    )
  }
  
  model_lookup <- model_fits |>
    dplyr::select(region, best_model)
  
  if (!all(regions_to_include %in% model_lookup$region)) {
    warning("Some regions in regions_to_include are missing from model_fits TSV.")
  }
  
  aln_per_region <- list()
  
  for (rg in regions_to_include) {
    
    rg_dir <- file.path(prep_dir, rg)
    
    trimmed_fa <- file.path(
      rg_dir,
      paste0(project_name, ".", region_set_name, "_", rg, ".trimmed.fasta")
    )
    
    if (!file.exists(trimmed_fa)) {
      stop("Missing trimmed alignment for region ", rg, ": ", trimmed_fa)
    }
    
    aln_per_region[[rg]] <- Biostrings::readDNAStringSet(trimmed_fa)
  }
  
  all_taxa <- sort(
    unique(
      unlist(
        lapply(aln_per_region, names)
      )
    )
  )
  
  # ------------------------------------------------------------
  # Write final taxon-by-gene coverage report
  # ------------------------------------------------------------
  
  coverage_df <- data.frame(
    taxon = all_taxa,
    stringsAsFactors = FALSE
  )
  
  for (rg in regions_to_include) {
    coverage_df[[rg]] <- as.integer(all_taxa %in% names(aln_per_region[[rg]]))
  }
  
  coverage_df$genes_present <- rowSums(
    coverage_df[, regions_to_include, drop = FALSE]
  )
  
  coverage_df$genes_missing <- length(regions_to_include) - coverage_df$genes_present
  
  coverage_path <- file.path(
    multi_gene_dir,
    paste0("final_taxon_gene_coverage_", project_name, ".", region_set_name, ".tsv")
  )
  
  write.table(
    coverage_df,
    coverage_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  message("Final taxon gene coverage written to: ", coverage_path)
  
  # ------------------------------------------------------------
  # Concatenate sequences
  # ------------------------------------------------------------
  
  concat_vec <- vapply(
    all_taxa,
    function(taxon) {
      paste0(
        vapply(
          regions_to_include,
          function(rg) {
            s <- aln_per_region[[rg]]
            
            if (taxon %in% names(s)) {
              as.character(s[[taxon]])
            } else {
              width_rg <- Biostrings::width(s)[1]
              paste(rep("-", width_rg), collapse = "")
            }
          },
          FUN.VALUE = character(1)
        ),
        collapse = ""
      )
    },
    FUN.VALUE = character(1)
  )
  
  concat_dna <- Biostrings::DNAStringSet(concat_vec)
  names(concat_dna) <- all_taxa
  
  region_lengths <- vapply(
    regions_to_include,
    function(rg) {
      unique(Biostrings::width(aln_per_region[[rg]]))[1]
    },
    FUN.VALUE = integer(1)
  )
  
  starts <- cumsum(c(1, head(region_lengths, -1)))
  ends <- cumsum(region_lengths)
  
  if (unique(Biostrings::width(concat_dna))[1] != tail(ends, 1)) {
    warning("Concatenated alignment length does not match sum of region lengths.")
  }
  
  concat_path <- file.path(
    multi_gene_dir,
    paste0("concatenated_", project_name, ".", region_set_name, ".fasta")
  )
  
  Biostrings::writeXStringSet(
    concat_dna,
    filepath = concat_path,
    format = "fasta"
  )
  
  message("Concatenated alignment written to: ", concat_path)
  
  models_ordered <- vapply(
    regions_to_include,
    function(rg) {
      
      row <- model_lookup[
        model_lookup$region == rg,
        ,
        drop = FALSE
      ]
      
      if (nrow(row) == 0L || is.na(row$best_model[1])) {
        warning(
          sprintf(
            "No best-fit model found for region '%s'. Using fallback model 'GTR+G'.",
            rg
          )
        )
        
        return("GTR+G")
      }
      
      row$best_model[1]
    },
    FUN.VALUE = character(1)
  )
  
  partition_lines <- c(
    "#nexus",
    "begin sets;"
  )
  
  part_names <- paste0("part", seq_along(regions_to_include))
  
  for (i in seq_along(regions_to_include)) {
    line <- sprintf(
      "\tcharset %s = %d-%d;",
      part_names[i],
      starts[i],
      ends[i]
    )
    
    partition_lines <- c(partition_lines, line)
  }
  
  part_specs <- paste(
    sprintf("%s:%s", models_ordered, part_names),
    collapse = ", "
  )
  
  partition_lines <- c(
    partition_lines,
    sprintf("\tcharpartition mine = %s;", part_specs),
    "end;"
  )
  
  part_path <- file.path(
    multi_gene_dir,
    paste0("partitions_", project_name, ".", region_set_name, ".nex")
  )
  
  writeLines(partition_lines, part_path)
  
  message("Partition NEXUS file written to: ", part_path)
  
  invisible(
    list(
      concat_fasta = concat_path,
      partitions_nex = part_path,
      coverage_tsv = coverage_path,
      coverage = coverage_df,
      regions = regions_to_include,
      starts = starts,
      ends = ends,
      models = models_ordered
    )
  )
}



iqtree_multigene_partitioned <- function(
    project_name,
    regions_to_include,
    run_dir = NULL,
    threads = 8,
    multigene_bootstraps = 1000,
    iqtree_args = c("-m", "MFP+MERGE"),
    force = TRUE
) {
  
  if (!exists("base_dir", envir = .GlobalEnv)) {
    stop("`base_dir` is not defined. Run start_project() first.")
  }
  
  iqtree_bin <- Sys.getenv("IQTREE_PATH")
  
  if (!nzchar(iqtree_bin)) {
    stop("IQTREE_PATH is not set in .Renviron.")
  }
  
  suppressWarnings(system2(iqtree_bin, "-version"))
  
  paths <- get_phylo_paths(
    project_name = project_name,
    regions_to_include = regions_to_include,
    run_dir = run_dir
  )
  
  region_set_name <- paths$region_set_name
  analysis_dir <- paths$analysis_dir
  multi_gene_dir <- paths$multi_gene_dir
  
  if (!dir.exists(analysis_dir)) {
    stop("Phylogeny run directory does not exist: ", analysis_dir)
  }
  
  if (!dir.exists(multi_gene_dir)) {
    stop("Multi-gene tree directory does not exist: ", multi_gene_dir)
  }
  
  concatenated_fasta <- file.path(
    multi_gene_dir,
    paste0("concatenated_", project_name, ".", region_set_name, ".fasta")
  )
  
  partition_nexus <- file.path(
    multi_gene_dir,
    paste0("partitions_", project_name, ".", region_set_name, ".nex")
  )
  
  iqtree_prefix <- file.path(
    multi_gene_dir,
    paste0("iqtree_", project_name, ".", region_set_name)
  )
  
  expected_treefile <- paste0(iqtree_prefix, ".treefile")
  
  if (!file.exists(concatenated_fasta)) {
    stop("Concatenated alignment not found: ", concatenated_fasta)
  }
  
  if (!file.exists(partition_nexus)) {
    stop("Partition Nexus file not found: ", partition_nexus)
  }
  
  if (file.exists(expected_treefile) && !force) {
    stop(
      "IQ-TREE treefile already exists and force = FALSE:\n  ",
      expected_treefile,
      "\nSet force = TRUE to overwrite or create a new run_dir."
    )
  }
  
  args <- c(
    "-s",
    concatenated_fasta,
    "-p",
    partition_nexus,
    "--ufboot",
    as.character(multigene_bootstraps)
  )
  
  if (!any(iqtree_args == "-nt")) {
    args <- c(args, "-nt", as.character(threads))
  }
  
  args <- c(
    args,
    iqtree_args,
    "--prefix",
    iqtree_prefix
  )
  
  message(
    "Running IQ-TREE with command:\n",
    iqtree_bin,
    " ",
    paste(shQuote(args), collapse = " ")
  )
  
  exit_status <- system2(
    command = iqtree_bin,
    args = args
  )
  
  if (exit_status != 0) {
    warning("IQ-TREE finished with non-zero exit status: ", exit_status)
  } else {
    message("IQ-TREE multigene partitioned run completed successfully.")
    
    if (file.exists(expected_treefile)) {
      message("Treefile: ", expected_treefile)
    }
  }
  
  invisible(
    list(
      status = exit_status,
      cmd = paste(iqtree_bin, paste(args, collapse = " ")),
      output_prefix = iqtree_prefix,
      treefile = expected_treefile
    )
  )
}






# ============================================================
# arborist extra wrappers
# ============================================================

# tree tip annotation

make_tree_host_annotation_table <- function(
    tree_file,
    metadata_file,
    host_ranks = c("kingdom"),
    output_dir = dirname(tree_file),
    output_prefix = NULL
) {
  if (!requireNamespace("ape", quietly = TRUE)) stop("Package 'ape' is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (!file.exists(tree_file)) stop("Tree file not found: ", tree_file)
  if (!file.exists(metadata_file)) stop("Metadata file not found: ", metadata_file)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  if (is.null(output_prefix)) {
    output_prefix <- tools::file_path_sans_ext(basename(tree_file))
  }
  
  host_cols <- paste0("Host.", host_ranks)
  
  meta <- read.csv(metadata_file, stringsAsFactors = FALSE)
  
  required_cols <- c(
    "Accession",
    "org_name",
    "strain.standard",
    "strain.standard.type",
    host_cols
  )
  
  missing_cols <- setdiff(required_cols, names(meta))
  if (length(missing_cols) > 0) {
    stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
  }
  
  clean_value <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    trimws(x)
  }
  
  normalize_tip_text <- function(x) {
    x <- as.character(x)
    x <- gsub("^>|^'", "", x)
    x <- gsub("'$", "", x)
    x <- gsub("[[:space:]]+", "_", x)
    x <- gsub("[^A-Za-z0-9_.-]+", "", x)
    x
  }
  
  collapse_unique_nonempty <- function(x) {
    x <- clean_value(x)
    x <- sort(unique(x[x != ""]))
    paste(x, collapse = ";")
  }
  
  count_unique_nonempty <- function(x) {
    x <- clean_value(x)
    length(unique(x[x != ""]))
  }
  
  meta$org_name <- clean_value(meta$org_name)
  meta$strain.standard <- clean_value(meta$strain.standard)
  meta$strain.standard.type <- clean_value(meta$strain.standard.type)
  
  for (col in host_cols) {
    meta[[col]] <- clean_value(meta[[col]])
  }
  
  meta$strain_label_for_tree <- ifelse(
    meta$strain.standard.type != "",
    meta$strain.standard.type,
    meta$strain.standard
  )
  
  meta$expected_tip_label <- paste0(
    gsub(" ", ".", meta$org_name),
    "_",
    meta$strain_label_for_tree
  )
  
  meta$expected_tip_label_norm <- normalize_tip_text(meta$expected_tip_label)
  
  strain_host_summary <- meta %>%
    dplyr::filter(
      org_name != "",
      strain_label_for_tree != ""
    ) %>%
    dplyr::group_by(
      org_name,
      strain_label_for_tree,
      expected_tip_label,
      expected_tip_label_norm
    ) %>%
    dplyr::summarise(
      accession_count = dplyr::n(),
      accessions = paste(sort(unique(Accession)), collapse = ";"),
      dplyr::across(
        dplyr::all_of(host_cols),
        collapse_unique_nonempty,
        .names = "{.col}_values"
      ),
      dplyr::across(
        dplyr::all_of(host_cols),
        count_unique_nonempty,
        .names = "{.col}_n_values"
      ),
      .groups = "drop"
    )
  
  for (rank in host_ranks) {
    values_col <- paste0("Host.", rank, "_values")
    n_col <- paste0("Host.", rank, "_n_values")
    assignment_col <- paste0("host_", rank)
    conflict_col <- paste0("host_", rank, "_conflict")
    
    strain_host_summary[[assignment_col]] <- dplyr::case_when(
      strain_host_summary[[n_col]] == 0 ~ NA_character_,
      strain_host_summary[[n_col]] == 1 ~ strain_host_summary[[values_col]],
      strain_host_summary[[n_col]] > 1 ~ "CONFLICT"
    )
    
    strain_host_summary[[conflict_col]] <- strain_host_summary[[n_col]] > 1
  }
  
  conflict_cols <- paste0("host_", host_ranks, "_conflict")
  
  conflict_report <- strain_host_summary %>%
    dplyr::filter(rowSums(dplyr::across(dplyr::all_of(conflict_cols))) > 0)
  
  tree <- ape::read.tree(tree_file)
  
  tip_df <- data.frame(
    tip_label = tree$tip.label,
    tip_label_norm = normalize_tip_text(tree$tip.label),
    stringsAsFactors = FALSE
  )
  
  annotation_rows <- lapply(seq_len(nrow(tip_df)), function(i) {
    tip <- tip_df$tip_label[i]
    tip_norm <- tip_df$tip_label_norm[i]
    
    exact_hits <- strain_host_summary[
      strain_host_summary$expected_tip_label_norm == tip_norm,
      ,
      drop = FALSE
    ]
    
    if (nrow(exact_hits) == 1) {
      hit <- exact_hits[1, ]
      match_type <- "exact"
    } else if (nrow(exact_hits) > 1) {
      hit <- exact_hits[1, ]
      match_type <- "multiple_exact_matches"
    } else {
      contained_hits <- strain_host_summary[
        vapply(
          strain_host_summary$expected_tip_label_norm,
          function(x) {
            grepl(x, tip_norm, fixed = TRUE) || grepl(tip_norm, x, fixed = TRUE)
          },
          logical(1)
        ),
        ,
        drop = FALSE
      ]
      
      if (nrow(contained_hits) == 1) {
        hit <- contained_hits[1, ]
        match_type <- "partial"
      } else if (nrow(contained_hits) > 1) {
        hit <- contained_hits[1, ]
        match_type <- "multiple_partial_matches"
      } else {
        out <- data.frame(
          tip_label = tip,
          matched_metadata_label = NA_character_,
          org_name = NA_character_,
          strain_label_for_tree = NA_character_,
          accession_count = NA_integer_,
          accessions = NA_character_,
          match_type = "unmatched",
          stringsAsFactors = FALSE
        )
        
        for (rank in host_ranks) {
          out[[paste0("host_", rank)]] <- NA_character_
          out[[paste0("host_", rank, "_conflict")]] <- NA
        }
        
        return(out)
      }
    }
    
    out <- data.frame(
      tip_label = tip,
      matched_metadata_label = hit$expected_tip_label,
      org_name = hit$org_name,
      strain_label_for_tree = hit$strain_label_for_tree,
      accession_count = hit$accession_count,
      accessions = hit$accessions,
      match_type = match_type,
      stringsAsFactors = FALSE
    )
    
    for (rank in host_ranks) {
      out[[paste0("host_", rank)]] <- hit[[paste0("host_", rank)]]
      out[[paste0("host_", rank, "_conflict")]] <- hit[[paste0("host_", rank, "_conflict")]]
    }
    
    out
  })
  
  annotation_table <- dplyr::bind_rows(annotation_rows)
  
  rank_tag <- paste(host_ranks, collapse = ".")
  
  annotation_path <- file.path(
    output_dir,
    paste0(output_prefix, "_tip_host_annotation_", rank_tag, ".csv")
  )
  
  conflict_path <- file.path(
    output_dir,
    paste0(output_prefix, "_host_conflicts_", rank_tag, ".csv")
  )
  
  unmatched_path <- file.path(
    output_dir,
    paste0(output_prefix, "_unmatched_tips_", rank_tag, ".csv")
  )
  
  write.csv(annotation_table, annotation_path, row.names = FALSE)
  write.csv(conflict_report, conflict_path, row.names = FALSE)
  write.csv(
    annotation_table[annotation_table$match_type == "unmatched", , drop = FALSE],
    unmatched_path,
    row.names = FALSE
  )
  
  message("Tip host annotation table written to: ", annotation_path)
  message("Host conflict report written to: ", conflict_path)
  message("Unmatched tip report written to: ", unmatched_path)
  message("Tips in tree: ", nrow(tip_df))
  message("Matched tips: ", sum(annotation_table$match_type != "unmatched"))
  message("Unmatched tips: ", sum(annotation_table$match_type == "unmatched"))
  message("Host conflicts: ", nrow(conflict_report))
  
  invisible(list(
    annotation_table = annotation_table,
    conflict_report = conflict_report,
    unmatched_tips = annotation_table[annotation_table$match_type == "unmatched", , drop = FALSE]
  ))
}

make_rank_palette <- function(x, palette_name = "Dark 3") {
  vals <- sort(unique(na.omit(as.character(x))))
  vals <- vals[vals != ""]
  
  cols <- grDevices::hcl.colors(
    n = length(vals),
    palette = palette_name
  )
  
  stats::setNames(cols, vals)
}

# aquiring accessions + data from NCBI
ncbi_data_fetch <- function(
    taxa_list,
    max_acc_per_taxa = "max",
    organism_scope = NULL,
    include_filters = NULL,
    exclude_filters = NULL,
    ncbi_database = "nucleotide",
    project_name = get0(
      "project_name",
      envir = .GlobalEnv
    ),
    accession_checkpoint_every = 500,
    accession_fetch_batch_size = 500,
    accession_page_max_retries = 3,
    accession_retry_wait = 5,
    accession_max_history_refreshes = 3,
    metadata_checkpoint_every = 500,
    metadata_progress_every = 5,
    resume = TRUE,
    overwrite_accession_checkpoints = FALSE,
    metadata_batch_size = 250
) {
  
  if (
    is.null(project_name) ||
    !nzchar(project_name)
  ) {
    
    stop(
      "project_name is not set. ",
      "Run start_project(project_name) first, ",
      "or pass project_name explicitly."
    )
  }
  
  
  # ============================================================
  # Resolve and remember NCBI database
  # ============================================================
  
  ncbi_database <- normalize_ncbi_database(
    ncbi_database
  )
  
  
  assign(
    "ncbi_database",
    ncbi_database,
    envir = .GlobalEnv
  )
  
  
  message(
    "\n============================================================"
  )
  
  message(
    "NCBI DATA SOURCE: ",
    toupper(ncbi_database)
  )
  
  message(
    "============================================================"
  )
  
  
  # Automatically locate and register an NCBI API key
  configure_entrez_key()
  
  
  # ============================================================
  # If the user explicitly wants accession checkpoints
  # overwritten, also remove any old completion marker.
  # ============================================================
  
  if (
    overwrite_accession_checkpoints &&
    file.exists(
      accession_completion_file()
    )
  ) {
    
    message(
      "Removing existing accession completion marker because ",
      "overwrite_accession_checkpoints = TRUE."
    )
    
    
    file.remove(
      accession_completion_file()
    )
  }
  
  
  # ============================================================
  # 1. Check whether accession retrieval is really complete
  # ============================================================
  
  retrieval_complete <- (
    resume &&
      !overwrite_accession_checkpoints &&
      accession_retrieval_is_complete(
        taxa_list = taxa_list,
        ncbi_database = ncbi_database
      )
  )
  
  
  # ============================================================
  # Reuse completed accession manifest
  # ============================================================
  
  if (retrieval_complete) {
    
    message(
      "\nAccession retrieval has already been completed for all ",
      length(taxa_list),
      " search group(s) from ",
      ncbi_database,
      "."
    )
    
    
    message(
      "Skipping individual accession checkpoint checks."
    )
    
    
    accession_path <-
      "./intermediate_files/all_pulled_accessions.csv"
    
    
    accession_list <- read.csv(
      accession_path,
      stringsAsFactors = FALSE,
      colClasses = "character"
    )
    
    
    if (
      !"Accession" %in%
      names(accession_list)
    ) {
      
      stop(
        "Existing accession manifest is missing the Accession column: ",
        accession_path
      )
    }
    
    
    if (
      "ncbi_database" %in%
      names(accession_list)
    ) {
      
      manifest_databases <- unique(
        accession_list$ncbi_database[
          !is.na(
            accession_list$ncbi_database
          ) &
            nzchar(
              accession_list$ncbi_database
            )
        ]
      )
      
      
      if (
        length(manifest_databases) > 0 &&
        any(
          manifest_databases !=
          ncbi_database
        )
      ) {
        
        stop(
          "Existing accession manifest belongs to a different NCBI database."
        )
      }
    }
    
    
    message(
      "Loaded ",
      nrow(accession_list),
      " unique accession(s) from existing manifest."
    )
    
    
  } else {
    
    # ==========================================================
    # The completion marker either does not exist or failed
    # validation.
    #
    # Remove only the stale marker.
    #
    # IMPORTANT:
    # Do NOT delete per-search-group accession checkpoints here.
    # They may contain valid partial progress.
    # ==========================================================
    
    if (
      file.exists(
        accession_completion_file()
      ) &&
      !overwrite_accession_checkpoints
    ) {
      
      message(
        "Removing stale or invalid accession completion marker."
      )
      
      
      file.remove(
        accession_completion_file()
      )
    }
    
    
    # ==========================================================
    # Retrieve / resume accessions
    # ==========================================================
    
    accession_list <-
      get_accessions_for_all_taxa(
        taxa_list =
          taxa_list,
        max_acc_per_taxa =
          max_acc_per_taxa,
        organism_scope =
          organism_scope,
        include_filters =
          include_filters,
        exclude_filters =
          exclude_filters,
        checkpoint_every =
          accession_checkpoint_every,
        resume =
          resume,
        overwrite =
          overwrite_accession_checkpoints,
        accession_fetch_batch_size =
          accession_fetch_batch_size,
        page_max_retries =
          accession_page_max_retries,
        retry_wait =
          accession_retry_wait,
        max_history_refreshes =
          accession_max_history_refreshes,
        ncbi_database =
          ncbi_database
      )
    
    
    # ==========================================================
    # Sanity check before declaring accession retrieval complete
    # ==========================================================
    
    if (
      is.null(accession_list) ||
      nrow(accession_list) == 0
    ) {
      
      stop(
        "Accession retrieval returned zero accessions. ",
        "A completion marker will NOT be written."
      )
    }
    
    
    # ==========================================================
    # Only write the completion marker after all search groups
    # have successfully completed.
    # ==========================================================
    
    write_accession_completion_marker(
      taxa_list =
        taxa_list,
      accession_file =
        "./intermediate_files/all_pulled_accessions.csv",
      ncbi_database =
        ncbi_database
    )
  }
  
  
  # ============================================================
  # Final accession sanity check before starting metadata
  # ============================================================
  
  if (
    is.null(accession_list) ||
    nrow(accession_list) == 0
  ) {
    
    stop(
      "Accession manifest contains zero accessions. ",
      "Metadata retrieval will not be started."
    )
  }
  
  
  # ============================================================
  # 2. Fetch metadata
  # ============================================================
  
  retrieve_ncbi_metadata(
    project_name =
      project_name,
    resume =
      resume,
    checkpoint_every =
      metadata_checkpoint_every,
    progress_every =
      metadata_progress_every,
    metadata_batch_size =
      metadata_batch_size,
    ncbi_database =
      ncbi_database
  )
  
  
  invisible(
    accession_list
  )
}


data_curate <- function(project_name = get0("project_name", envir = .GlobalEnv, ifnotfound = NULL),
                        taxa_of_interest = get0("taxa_of_interest", envir = .GlobalEnv, ifnotfound = NULL),
                        my_lab_sequences = get0("my_lab_sequences", envir = .GlobalEnv, ifnotfound = ""),
                        add_strain_taxonomy = TRUE,
                        overwrite_strain_taxonomy = TRUE) {
  
  if (is.null(project_name) || !nzchar(project_name)) {
    stop("project_name is not set. Run start_project(project_name) first, or pass project_name explicitly.")
  }
  
  message("Starting metadata curation for project: ", project_name)
  
  # 1) Merge custom/lab-generated sequences and metadata before curation
  message("\nStep 1: Checking for custom sequences/data...")
  
  merge_metadata_with_custom_file(
    project_name = project_name,
    my_lab_sequences = my_lab_sequences
  )
  
  # 2) Basic metadata curation
  # This creates:
  #   strain.standard
  #   strain.standard.type
  #   org_name
  #
  # If add_strain_taxonomy = TRUE, curate_metadata_basic() also adds:
  #   Strain.taxonomy
  #   Strain.phylum
  #   Strain.class
  #   Strain.order
  #   Strain.family
  #   Strain.genus
  #   Strain.species
  message("\nStep 2: Running basic metadata curation...")
  
  curated_data <- curate_metadata_basic(
    project_name = project_name,
    taxa_of_interest = taxa_of_interest,
    add_strain_taxonomy = add_strain_taxonomy,
    overwrite_strain_taxonomy = overwrite_strain_taxonomy
  )
  
  final_path <- file.path(
    "./metadata_files",
    paste0("all_accessions_pulled_metadata_", project_name, "_curated.csv")
  )
  
  if (!file.exists(final_path)) {
    warning("Expected final curated metadata file was not found: ", final_path)
  } else {
    message("\nFinal curated metadata written to: ", final_path)
  }
  
  message("\nMetadata curation complete.")
  message("Region curation remains separate. Run curate_metadata_regions(project_name) when needed.")
  
  invisible(curated_data)
}


