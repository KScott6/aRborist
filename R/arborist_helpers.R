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
compose_entrez_term <- function(taxon,
                                organism_scope  = NULL,
                                include_filters = NULL,
                                exclude_filters = NULL) {
  # 1) Resolve arguments in priority:
  #    explicit arg > .GlobalEnv override > package default
  if (is.null(organism_scope)) {
    organism_scope <- get0("organism_scope", envir = .GlobalEnv,
                           ifnotfound = default_organism_scope)
  }
  if (is.null(include_filters)) {
    include_filters <- get0("search_include", envir = .GlobalEnv,
                            ifnotfound = default_search_include)
  }
  if (is.null(exclude_filters)) {
    exclude_filters <- get0("search_exclude", envir = .GlobalEnv,
                            ifnotfound = default_search_exclude)
  }
  
  # 2) Base organism term
  q <- sprintf('"%s"[Organism]', taxon)
  
  # 3) Add organism_scope (e.g. txid4751[Organism:exp])
  if (!is.null(organism_scope) && nzchar(organism_scope)) {
    q <- paste(q, organism_scope, sep = " AND ")
  }
  
  # 4) Add positive filters with AND
  include_filters <- include_filters[nzchar(include_filters)]
  if (length(include_filters) > 0) {
    include_str <- paste(include_filters, collapse = " AND ")
    q <- paste(q, include_str, sep = " AND ")
  }
  
  # 5) Add negative filters with NOT (no leading AND)
  exclude_filters <- exclude_filters[nzchar(exclude_filters)]
  if (length(exclude_filters) > 0) {
    not_str <- paste(paste("NOT", exclude_filters), collapse = " ")
    q <- paste(q, not_str)
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

# Sleep helper (uses ncbi_api_key from global env)
# this is not using the "10 requests/sec with API, 3 requests/sec without" timing because I noticed the requests were being bunched up and sent in groups, resulting in a noticable percentage of my requests getting denied. 
# you can mess with the timings if you want, but watch out for denied requests
get_sleep_duration <- function() {
  if (!is.null(ncbi_api_key) && nzchar(ncbi_api_key)) 0.2 else 0.5
}
if (exists("ncbi_api_key") && !is.null(ncbi_api_key) && nzchar(ncbi_api_key)) {
  rentrez::set_entrez_key(ncbi_api_key)
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

# Fetch accessions from NCBI - searching using entrez
fetch_accessions_for_taxon <- function(taxon,
                                       max_acc        = max_acc_per_taxa,
                                       organism_scope = NULL,
                                       include_filters = NULL,
                                       exclude_filters = NULL) {
  cat("Searching term:", taxon, "\n")
  
  # Resolve defaults if not provided
  if (is.null(organism_scope) && exists("organism_scope", .GlobalEnv))
    organism_scope <- get("organism_scope", .GlobalEnv)
  if (is.null(include_filters) && exists("search_include", .GlobalEnv))
    include_filters <- get("search_include", .GlobalEnv)
  if (is.null(exclude_filters) && exists("search_exclude", .GlobalEnv))
    exclude_filters <- get("search_exclude", .GlobalEnv)
  
  # Build full query string
  filters <- compose_entrez_term(
    taxon            = taxon,
    organism_scope   = organism_scope,
    include_filters  = include_filters,
    exclude_filters  = exclude_filters
  )
  
  search <- rentrez::entrez_search(db = "nucleotide", term = filters, retmax = 9999)
  
  if (length(search$ids) == 0) {
    cat("No accessions found for:", taxon, "\n")
    return(data.frame(Accession = character(0), genus = character(0), stringsAsFactors = FALSE))
  }
  
  max_n <- if (identical(max_acc, "max")) Inf else as.numeric(max_acc)
  
  if (length(search$ids) == 9999) {
    cat(taxon, "has ≥ 10,000 NCBI accessions. Using webhistory.\n")
    large_search <- rentrez::entrez_search(db = "nucleotide", term = filters, use_history = TRUE)
    total_accession_count <- as.integer(large_search[["count"]])
    pull_n <- min(total_accession_count, max_n)
    cat(total_accession_count, "accessions available for", taxon, "- pulling a maximum of", pull_n, "\n")
    
    tmp <- paste0("./intermediate_files/temp_file_accessions_from_", taxon, ".txt")
    if (file.exists(tmp)) file.remove(tmp)
    
    for (seq_start in seq(0, pull_n - 1, by = 50)) {
      recs <- rentrez::entrez_fetch(
        db = "nuccore",
        web_history = large_search$web_history,
        rettype = "acc",
        retmax = min(50, pull_n - seq_start),
        retstart = seq_start
      )
      cat(recs, file = tmp, append = TRUE)
      Sys.sleep(get_sleep_duration())
    }
    
    df <- read.table(tmp, stringsAsFactors = FALSE)
    colnames(df) <- "Accession"
    df$genus <- taxon
    file.remove(tmp)
    cat("Accession retrieval for", taxon, "successful\n\n")
    return(df)
  }
  
  ids <- search$ids
  if (is.finite(max_n)) {
    ids <- utils::head(ids, max_n)
  }
  
  if (length(ids) <= 300) {
    summary <- rentrez::entrez_summary(db = "nuccore", id = ids)
    Sys.sleep(get_sleep_duration())
  } else {
    summary <- list()
    idx <- split(seq_along(ids), ceiling(seq_along(ids) / 300))
    for (p in idx) {
      summary[p] <- rentrez::entrez_summary(db = "nuccore", id = ids[p])
      Sys.sleep(get_sleep_duration())
    }
    class(summary) <- c("esummary_list", "list")
  }
  
  tempdf <- data.frame(
    Accession = unname(rentrez::extract_from_esummary(summary, "caption")),
    stringsAsFactors = FALSE
  )
  tempdf$genus <- taxon
  cat("Search complete for", taxon, " (", nrow(tempdf), " accessions)\n", sep = "")
  tempdf
}

# pulling accessions using provided taxa names
get_accessions_for_all_taxa <- function(taxa_list,
                                        max_acc_per_taxa,
                                        organism_scope  = NULL,
                                        include_filters = NULL,
                                        exclude_filters = NULL,
                                        timing_file = "./intermediate_files/fetch_times_accessions.csv") {
  
  # Resolve defaults if not provided
  if (is.null(organism_scope) && exists("organism_scope", .GlobalEnv))
    organism_scope <- get("organism_scope", envir = .GlobalEnv)
  if (is.null(include_filters) && exists("search_include", .GlobalEnv))
    include_filters <- get("search_include", envir = .GlobalEnv)
  if (is.null(exclude_filters) && exists("search_exclude", .GlobalEnv))
    exclude_filters <- get("search_exclude", envir = .GlobalEnv)
  
  taxa_frame_acc <- vector("list", length(taxa_list))
  
  timing_log <- data.frame(Taxon = character(),
                           Num_accessions = integer(),
                           Start_time = character(),
                           End_time = character(),
                           Elapsed_minutes = numeric(),
                           stringsAsFactors = FALSE)
  
  overall_start <- Sys.time()
  
  for (i in seq_along(taxa_list)) {
    term <- taxa_list[i]
    cat("\n=== Starting", term, "(", i, "of", length(taxa_list), ") ===\n")
    start_time <- Sys.time()
    
    tempdf <- tryCatch({
      fetch_accessions_for_taxon(
        taxon           = term,
        max_acc         = max_acc_per_taxa,
        organism_scope  = organism_scope,
        include_filters = include_filters,
        exclude_filters = exclude_filters
      )
    }, error = function(e) {
      cat("ERROR:", conditionMessage(e), "\n")
      data.frame(Accession = character(0), genus = character(0), stringsAsFactors = FALSE)
    })
    
    end_time <- Sys.time()
    elapsed <- as.numeric(difftime(end_time, start_time, units = "mins"))
    cat(sprintf("Finished %s in %.2f minutes\n", term, elapsed))
    
    timing_log <- rbind(timing_log, data.frame(
      Taxon = term,
      Num_accessions = nrow(tempdf),
      Start_time = format(start_time, "%Y-%m-%d %H:%M:%S"),
      End_time = format(end_time, "%Y-%m-%d %H:%M:%S"),
      Elapsed_minutes = round(elapsed, 2),
      stringsAsFactors = FALSE
    ))
    
    outfile_name <- paste0("./intermediate_files/Accessions_for_", term, ".csv")
    write.csv(tempdf, outfile_name, row.names = FALSE, quote = FALSE)
    taxa_frame_acc[[i]] <- tempdf
  }
  
  total_elapsed <- as.numeric(difftime(Sys.time(), overall_start, units = "mins"))
  cat("\nAll taxa completed in", round(total_elapsed, 2), "minutes.\n")
  
  non_empty <- Filter(function(x) nrow(x) > 0, taxa_frame_acc)
  accession_list <- if (length(non_empty)) {
    do.call(rbind, non_empty)
  } else {
    data.frame(Accession = character(0), genus = character(0), stringsAsFactors = FALSE)
  }
  
  write.csv(accession_list, "./intermediate_files/all_pulled_accessions.csv", row.names = FALSE)
  write.csv(timing_log, timing_file, row.names = FALSE)
  cat("Timing log written to ", timing_file, "\n", sep = "")
  
  return(accession_list)
}


# Metadata retrieval, using list of pulled accessions
fetch_metadata_for_accession <- function(accession) {
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
  
  x <- rentrez::entrez_fetch(db = "nuccore", id = accession, rettype = "xml")
  doc <- XML::xmlParse(x)
  
  # top-level fields
  top_locus     <- XML::xpathSApply(doc, "//GBSeq_locus", xmlValue)
  top_len       <- XML::xpathSApply(doc, "//GBSeq_length", xmlValue)
  top_strand    <- XML::xpathSApply(doc, "//GBSeq_strandedness", xmlValue)
  top_moltype   <- XML::xpathSApply(doc, "//GBSeq_moltype", xmlValue)
  top_upd       <- XML::xpathSApply(doc, "//GBSeq_update-date", xmlValue)
  top_create    <- XML::xpathSApply(doc, "//GBSeq_create-date", xmlValue)
  top_def       <- XML::xpathSApply(doc, "//GBSeq_definition", xmlValue)
  top_accver    <- XML::xpathSApply(doc, "//GBSeq_accession-version", xmlValue)
  top_proj      <- XML::xpathSApply(doc, "//GBSeq_project", xmlValue)
  top_org       <- XML::xpathSApply(doc, "//GBSeq_organism", xmlValue)
  top_tax       <- XML::xpathSApply(doc, "//GBSeq_taxonomy", xmlValue)
  top_seq       <- XML::xpathSApply(doc, "//GBSeq_sequence", xmlValue)
  
  # qualifiers
  q_names  <- XML::xpathSApply(doc, "//GBQualifier/GBQualifier_name",  xmlValue)
  q_values <- XML::xpathSApply(doc, "//GBQualifier/GBQualifier_value", xmlValue)
  quals <- data.frame(name = q_names, value = q_values, stringsAsFactors = FALSE)
  
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

# metadata retrieval
retrieve_ncbi_metadata <- function(project_name,
                                   resume = TRUE,
                                   overwrite_checkpoints = FALSE,
                                   checkpoint_dir = "./metadata_files/metadata_checkpoints") {
  
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
  
  accession_list <- read.csv(accession_path, header = TRUE, stringsAsFactors = FALSE)
  
  if (!"Accession" %in% names(accession_list)) {
    stop("Accession list must contain an 'Accession' column.")
  }
  
  if ("genus" %in% names(accession_list)) {
    taxa_groups <- split(accession_list, accession_list$genus)
  } else {
    taxa_groups <- list(ALL = accession_list)
  }
  
  failed_path <- file.path(checkpoint_dir, "metadata_failed_accessions.csv")
  timing_path <- "./intermediate_files/fetch_times_metadata_by_taxon.csv"
  
  failed_log <- data.frame(
    Taxon = character(),
    Accession = character(),
    Error = character(),
    Time = character(),
    stringsAsFactors = FALSE
  )
  
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
  
  for (tx in names(taxa_groups)) {
    
    safe_tx <- gsub("[^A-Za-z0-9_.-]+", "_", tx)
    
    checkpoint_path <- file.path(
      checkpoint_dir,
      paste0("metadata_", safe_tx, ".csv")
    )
    
    block <- taxa_groups[[tx]]
    
    if (file.exists(checkpoint_path) && resume && !overwrite_checkpoints) {
      message("\n--- Skipping ", tx, ": checkpoint already exists ---")
      next
    }
    
    if (file.exists(checkpoint_path) && overwrite_checkpoints) {
      message("\n--- Overwriting existing checkpoint for ", tx, " ---")
      file.remove(checkpoint_path)
    }
    
    taxon_start <- Sys.time()
    
    message("\n--- ", tx, ": ", nrow(block), " accession(s) ---")
    
    metadata_rows <- list()
    success_count <- 0L
    fail_count <- 0L
    
    for (i in seq_len(nrow(block))) {
      
      acc <- block$Accession[i]
      
      message("[", i, "/", nrow(block), "] ", acc, " ...")
      
      entry <- tryCatch(
        {
          fetch_metadata_for_accession(acc)
        },
        error = function(e) {
          fail_count <<- fail_count + 1L
          
          failed_log <<- rbind(
            failed_log,
            data.frame(
              Taxon = tx,
              Accession = acc,
              Error = conditionMessage(e),
              Time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              stringsAsFactors = FALSE
            )
          )
          
          message("  ERROR: ", conditionMessage(e))
          NULL
        }
      )
      
      if (!is.null(entry)) {
        success_count <- success_count + 1L
        metadata_rows[[length(metadata_rows) + 1L]] <- entry
        
        message(
          "  OK | Species: ", entry$organism,
          " | Strain: ",
          dplyr::coalesce(entry$strain, entry$specimen_voucher, entry$isolate, ""),
          " | Host: ", entry$host
        )
      }
      
      Sys.sleep(get_sleep_duration())
    }
    
    if (length(metadata_rows) > 0) {
      taxon_metadata <- plyr::rbind.fill(metadata_rows)
      write.csv(taxon_metadata, checkpoint_path, row.names = FALSE)
      message("Checkpoint written: ", checkpoint_path)
    } else {
      warning("No metadata successfully retrieved for taxon: ", tx)
    }
    
    if (nrow(failed_log) > 0) {
      write.csv(failed_log, failed_path, row.names = FALSE)
      message("Failed accession log written: ", failed_path)
    }
    
    taxon_end <- Sys.time()
    taxon_elapsed <- as.numeric(difftime(taxon_end, taxon_start, units = "mins"))
    
    timing_log <- rbind(
      timing_log,
      data.frame(
        Taxon = tx,
        Num_accessions = nrow(block),
        Num_successful = success_count,
        Num_failed = fail_count,
        Start_time = format(taxon_start, "%Y-%m-%d %H:%M:%S"),
        End_time = format(taxon_end, "%Y-%m-%d %H:%M:%S"),
        Elapsed_minutes = round(taxon_elapsed, 2),
        Checkpoint_file = checkpoint_path,
        stringsAsFactors = FALSE
      )
    )
    
    write.csv(timing_log, timing_path, row.names = FALSE)
    
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
  
  checkpoint_files <- list.files(
    checkpoint_dir,
    pattern = "^metadata_.*\\.csv$",
    full.names = TRUE
  )
  
  checkpoint_files <- checkpoint_files[
    !grepl("metadata_failed_accessions\\.csv$", checkpoint_files)
  ]
  
  if (length(checkpoint_files) == 0) {
    stop("No checkpoint metadata files found in: ", checkpoint_dir)
  }
  
  message("\nCombining ", length(checkpoint_files), " checkpoint file(s).")
  
  metadata_database <- plyr::rbind.fill(
    lapply(checkpoint_files, function(x) {
      read.csv(x, stringsAsFactors = FALSE)
    })
  )
  
  metadata_database <- metadata_database %>%
    dplyr::distinct(Accession, .keep_all = TRUE)
  
  final_path <- paste0(
    "./metadata_files/all_accessions_pulled_metadata_",
    project_name,
    ".csv"
  )
  
  write.csv(metadata_database, final_path, row.names = FALSE)
  
  overall_end <- Sys.time()
  total_elapsed <- as.numeric(difftime(overall_end, overall_start, units = "mins"))
  
  timing_log <- rbind(
    timing_log,
    data.frame(
      Taxon = "TOTAL",
      Num_accessions = nrow(accession_list),
      Num_successful = nrow(metadata_database),
      Num_failed = if (file.exists(failed_path)) nrow(read.csv(failed_path)) else 0L,
      Start_time = format(overall_start, "%Y-%m-%d %H:%M:%S"),
      End_time = format(overall_end, "%Y-%m-%d %H:%M:%S"),
      Elapsed_minutes = round(total_elapsed, 2),
      Checkpoint_file = final_path,
      stringsAsFactors = FALSE
    )
  )
  
  write.csv(timing_log, timing_path, row.names = FALSE)
  
  message("\nMetadata retrieval complete.")
  message("Final metadata written to: ", final_path)
  message("Timing log written to: ", timing_path)
  
  if (file.exists(failed_path)) {
    message("Failed accession log written to: ", failed_path)
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
  
  if (is.null(my_lab_sequences) || !nzchar(my_lab_sequences)) {
    message("No custom sequences file provided; skipping custom sequence merge.")
    return(invisible(NULL))
  }
  
  if (!file.exists(my_lab_sequences)) {
    stop("Custom sequences file not found: ", my_lab_sequences)
  }
  
  metadata_database <- read.csv(metadata_file_path, stringsAsFactors = FALSE)
  custom_sequences  <- read.csv(my_lab_sequences, stringsAsFactors = FALSE, fill = TRUE)
  
  required_cols <- c("Accession", "strain", "sequence", "organism", "gene")
  missing_required <- setdiff(required_cols, names(custom_sequences))
  
  if (length(missing_required) > 0) {
    stop(
      "Custom file is missing required column(s): ",
      paste(missing_required, collapse = ", "),
      "\nRequired columns are: ",
      paste(required_cols, collapse = ", ")
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
  
  merged_data <- plyr::rbind.fill(metadata_database, custom_sequences)
  
  merged_data <- merged_data %>%
    dplyr::distinct(Accession, .keep_all = TRUE)
  
  write.csv(merged_data, metadata_file_path, row.names = FALSE)
  
  message("Merged custom sequences into: ", metadata_file_path)
  message("Custom rows added from: ", my_lab_sequences)
  
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
  if (!is.null(taxa_of_interest) && length(taxa_of_interest) > 0) {
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
add_failed_host_term_counts <- function(failed_df,
                                        project_name,
                                        metadata_dir = "./metadata_files") {
  metadata_path <- file.path(
    metadata_dir,
    paste0("all_accessions_pulled_metadata_", project_name, "_curated.csv")
  )
  
  if (!file.exists(metadata_path)) {
    warning("Cannot add failed host term counts; metadata file not found: ", metadata_path)
    failed_df$accession_count <- NA_integer_
    return(failed_df)
  }
  
  meta <- read.csv(metadata_path, stringsAsFactors = FALSE)
  
  if (!"host.standardized" %in% names(meta)) {
    warning("Cannot add failed host term counts; metadata is missing host.standardized.")
    failed_df$accession_count <- NA_integer_
    return(failed_df)
  }
  
  # If Accession exists, count unique accessions.
  # Otherwise, fall back to row counts.
  if ("Accession" %in% names(meta)) {
    count_df <- meta %>%
      dplyr::filter(!is.na(host.standardized) & host.standardized != "") %>%
      dplyr::group_by(host.standardized) %>%
      dplyr::summarise(
        accession_count = dplyr::n_distinct(Accession),
        .groups = "drop"
      )
  } else {
    count_df <- meta %>%
      dplyr::filter(!is.na(host.standardized) & host.standardized != "") %>%
      dplyr::count(host.standardized, name = "accession_count")
  }
  
  failed_df <- failed_df %>%
    dplyr::select(-dplyr::any_of("accession_count")) %>%
    dplyr::left_join(
      count_df,
      by = c("original_term" = "host.standardized")
    ) %>%
    dplyr::mutate(
      accession_count = dplyr::if_else(
        is.na(accession_count),
        0L,
        as.integer(accession_count)
      )
    ) %>%
    dplyr::arrange(
      dplyr::desc(accession_count),
      original_term
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
  
  fungal_col <- paste0("Strain.", fungal_rank)
  host_col   <- paste0("Host.", host_rank)
  
  if (!fungal_col %in% names(meta)) {
    stop(
      "Column '", fungal_col, "' not found in metadata.\n",
      "Available Strain.* columns: ",
      paste(grep("^Strain\\.", names(meta), value = TRUE), collapse = ", ")
    )
  }
  
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
  
  names(counts_wide)[1] <- fungal_col
  host_cols <- setdiff(names(counts_wide), fungal_col)
  
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
      failed_df[[col]] <- NA_character_
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

.read_host_taxonomy <- function(project_name, host_dir = "./host_assessment") {
  taxonomy_path <- .host_taxonomy_path(project_name, host_dir)
  taxonomy_cols <- .host_taxonomy_cols()
  
  if (file.exists(taxonomy_path)) {
    host_taxonomy <- read.csv(taxonomy_path, stringsAsFactors = FALSE)
    
    if (!"Host.standardized" %in% names(host_taxonomy)) {
      if ("Host.standard" %in% names(host_taxonomy)) {
        message("Detected legacy taxonomy file. Renaming 'Host.standard' to 'Host.standardized'.")
        names(host_taxonomy)[names(host_taxonomy) == "Host.standard"] <- "Host.standardized"
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
    
    host_taxonomy <- host_taxonomy[, taxonomy_cols, drop = FALSE]
  } else {
    host_taxonomy <- .empty_host_taxonomy()
  }
  
  host_taxonomy$Host.standardized <- .host_clean_term(host_taxonomy$Host.standardized)
  
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
  
  terms <- .host_display_term(meta$host.standardized.original)
  
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

.lookup_one_host_taxonomy <- function(term, db = "ncbi") {
  ranks_of_interest <- c(
    "kingdom", "phylum", "class",
    "order", "family", "genus", "species"
  )
  
  res_list <- tryCatch(
    {
      taxize::classification(term, db = db)
    },
    error = function(e) {
      NULL
    }
  )
  
  if (is.null(res_list) || length(res_list) == 0 || is.atomic(res_list)) {
    return(NULL)
  }
  
  res <- res_list[[1]]
  
  if (!.is_valid_tax_table(res)) {
    return(NULL)
  }
  
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
  
  df_row[, .host_taxonomy_cols(), drop = FALSE]
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
    sleep_sec = 0.1,
    overwrite = FALSE,
    term_type = "initial_lookup",
    parent_map = NULL
) {
  if (!dir.exists(host_dir)) dir.create(host_dir, recursive = TRUE)
  
  host_taxonomy <- .read_host_taxonomy(project_name, host_dir)
  failed_df <- .read_failed_terms(project_name, host_dir)
  
  terms <- unique(.host_clean_term(terms))
  terms <- terms[!is.na(terms)]
  
  invalid_terms <- terms[.host_is_invalid_tax_term(terms)]
  valid_terms <- terms[!.host_is_invalid_tax_term(terms)]
  
  if (length(invalid_terms) > 0) {
    failed_df <- .append_failed_terms(
      failed_df,
      invalid_terms,
      term_type = term_type,
      lookup_status = "not_searched_invalid"
    )
  }
  
  already_done <- unique(.host_clean_term(host_taxonomy$Host.standardized))
  
  if (!overwrite) {
    valid_terms <- setdiff(valid_terms, already_done)
  }
  
  if (length(valid_terms) == 0) {
    message("No new valid host terms to query.")
    .write_failed_terms(
      .add_failed_host_term_counts(failed_df, project_name),
      project_name,
      host_dir
    )
    return(invisible(list(
      taxonomy = host_taxonomy,
      failed_table = failed_df,
      newly_failed = character(0),
      newly_successful = character(0)
    )))
  }
  
  message("Host taxonomy lookup starting for ", length(valid_terms), " terms.")
  
  successful_rows <- list()
  failed_terms <- character(0)
  
  for (term in valid_terms) {
    message("  Querying: ", term)
    
    df_row <- .lookup_one_host_taxonomy(term, db = db)
    
    if (is.null(df_row)) {
      message("    FAILED")
      failed_terms <- c(failed_terms, term)
    } else {
      message("    OK")
      successful_rows[[length(successful_rows) + 1L]] <- df_row
    }
    
    Sys.sleep(sleep_sec)
  }
  
  if (length(successful_rows) > 0) {
    new_tax_rows <- dplyr::bind_rows(successful_rows)
    
    if (overwrite) {
      host_taxonomy <- host_taxonomy[
        !host_taxonomy$Host.standardized %in% new_tax_rows$Host.standardized,
        ,
        drop = FALSE
      ]
    }
    
    host_taxonomy <- dplyr::bind_rows(host_taxonomy, new_tax_rows)
    host_taxonomy <- host_taxonomy[!duplicated(host_taxonomy$Host.standardized), ]
  }
  
  if (length(failed_terms) > 0) {
    failed_df <- .append_failed_terms(
      failed_df,
      failed_terms,
      term_type = term_type,
      lookup_status = "failed"
    )
  }
  
  failed_df <- .add_failed_host_term_counts(failed_df, project_name)
  
  .write_host_taxonomy(host_taxonomy, project_name, host_dir)
  .write_failed_terms(failed_df, project_name, host_dir)
  
  invisible(list(
    taxonomy = host_taxonomy,
    failed_table = failed_df,
    newly_failed = unique(failed_terms),
    newly_successful = if (length(successful_rows) > 0) {
      unique(dplyr::bind_rows(successful_rows)$Host.standardized)
    } else {
      character(0)
    }
  ))
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
    sleep_sec = 0.1,
    overwrite = FALSE,
    metadata_file = NULL,
    use_isolation_source = FALSE,
    overwrite_host_standardized = FALSE
) {
  if (!dir.exists(host_dir)) dir.create(host_dir, recursive = TRUE)
  
  if (is.null(metadata_file)) {
    metadata_file <- .host_metadata_path(project_name)
  }
  
  if (!file.exists(metadata_file)) {
    stop("Metadata file not found: ", metadata_file)
  }
  
  message("Running host assessment initial pass...")
  
  meta <- read.csv(metadata_file, stringsAsFactors = FALSE)
  
  meta <- .ensure_host_standardized(
    meta,
    use_isolation_source = use_isolation_source,
    overwrite_host_standardized = overwrite_host_standardized
  )
  
  write.csv(meta, metadata_file, row.names = FALSE)
  
  host_terms <- unique(.host_clean_term(meta$host.standardized))
  
  terms_path <- file.path(
    host_dir,
    paste0("host_terms_for_taxonomy_", project_name, ".csv")
  )
  
  write.csv(
    data.frame(host = host_terms, stringsAsFactors = FALSE),
    terms_path,
    row.names = FALSE
  )
  
  message("Host terms file written to: ", terms_path)
  
  lookup_result <- .search_host_terms(
    terms = host_terms,
    project_name = project_name,
    host_dir = host_dir,
    db = db,
    sleep_sec = sleep_sec,
    overwrite = overwrite,
    term_type = "initial_lookup"
  )
  
  merge_host_taxonomy_into_metadata(
    project_name = project_name,
    host_dir = host_dir,
    metadata_file = metadata_file
  )
  
  message("Initial host assessment pass completed.")
  
  invisible(lookup_result)
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
    sleep_sec = 0.1,
    overwrite = FALSE,
    metadata_file = NULL
) {
  if (!dir.exists(host_dir)) dir.create(host_dir, recursive = TRUE)
  
  if (is.null(metadata_file)) {
    metadata_file <- .host_metadata_path(project_name)
  }
  
  if (!file.exists(metadata_file)) {
    stop("Metadata file not found: ", metadata_file)
  }
  
  failed_path <- .host_failed_path(project_name, host_dir)
  
  if (!file.exists(failed_path)) {
    stop(
      "Failed terms file not found: ", failed_path,
      "\nRun run_host_assessment_initial_pass() first."
    )
  }
  
  message("Running host assessment refinement pass...")
  
  failed_df <- .read_failed_terms(project_name, host_dir)
  host_taxonomy <- .read_host_taxonomy(project_name, host_dir)
  
  replacement_ok <- !.host_is_invalid_tax_term(failed_df$replacement_term)
  replacement_rows <- failed_df[replacement_ok, , drop = FALSE]
  
  if (nrow(replacement_rows) == 0) {
    message("No replacement terms supplied yet. Nothing to refine.")
    
    failed_df <- .add_failed_host_term_counts(failed_df, project_name, metadata_file)
    .write_failed_terms(failed_df, project_name, host_dir)
    
    merge_host_taxonomy_into_metadata(
      project_name = project_name,
      host_dir = host_dir,
      metadata_file = metadata_file
    )
    
    return(invisible(list(
      taxonomy = host_taxonomy,
      failed_table = failed_df
    )))
  }
  
  replacement_terms <- unique(.host_clean_term(replacement_rows$replacement_term))
  replacement_terms <- replacement_terms[!.host_is_invalid_tax_term(replacement_terms)]
  
  usable_taxonomy_rows <- .has_usable_host_taxonomy(host_taxonomy)
  
  already_have_taxonomy <- unique(
    .host_clean_term(host_taxonomy$Host.standardized[usable_taxonomy_rows])
  )
  
  if (overwrite) {
    terms_to_query <- replacement_terms
  } else {
    terms_to_query <- setdiff(replacement_terms, already_have_taxonomy)
  }
  
  message("Replacement terms supplied: ", length(replacement_terms))
  message("Replacement terms already have usable taxonomy: ", length(intersect(replacement_terms, already_have_taxonomy)))
  message("Replacement terms to query this pass: ", length(terms_to_query))
  
  if (length(terms_to_query) > 0) {
    lookup_result <- .search_host_terms(
      terms = terms_to_query,
      project_name = project_name,
      host_dir = host_dir,
      db = db,
      sleep_sec = sleep_sec,
      overwrite = overwrite,
      term_type = "replacement_lookup"
    )
  } else {
    lookup_result <- list(
      taxonomy = host_taxonomy,
      failed_table = failed_df,
      newly_failed = character(0),
      newly_successful = character(0)
    )
  }
  
  host_taxonomy <- .read_host_taxonomy(project_name, host_dir)
  
  usable_taxonomy_rows <- .has_usable_host_taxonomy(host_taxonomy)
  
  successful_taxonomy_terms <- unique(
    .host_clean_term(host_taxonomy$Host.standardized[usable_taxonomy_rows])
  )
  
  meta <- read.csv(metadata_file, stringsAsFactors = FALSE)
  meta <- .ensure_host_standardized(meta)
  
  if (!"host.standardized.before_refinement" %in% names(meta)) {
    meta$host.standardized.before_refinement <- meta$host.standardized
  }
  
  message("Applying successful replacement terms to metadata...")
  
  for (i in seq_len(nrow(replacement_rows))) {
    original <- .host_display_term(replacement_rows$original_term[i])
    replacement <- .host_clean_term(replacement_rows$replacement_term[i])
    
    if (.host_is_invalid_tax_term(replacement)) {
      next
    }
    
    replacement_has_taxonomy <- replacement %in% successful_taxonomy_terms
    
    if (!replacement_has_taxonomy) {
      next
    }
    
    if (original == "NA") {
      idx <- .host_is_invalid_tax_term(meta$host.standardized)
    } else {
      idx <- .host_display_term(meta$host.standardized) == original |
        .host_display_term(meta$host.standardized.original) == original
    }
    
    if (any(idx, na.rm = TRUE)) {
      meta$host.standardized[idx] <- replacement
    }
  }
  
  write.csv(meta, metadata_file, row.names = FALSE)
  
  failed_df <- .read_failed_terms(project_name, host_dir)
  host_taxonomy <- .read_host_taxonomy(project_name, host_dir)
  
  usable_taxonomy_rows <- .has_usable_host_taxonomy(host_taxonomy)
  
  successful_taxonomy_terms <- unique(
    .host_clean_term(host_taxonomy$Host.standardized[usable_taxonomy_rows])
  )
  
  replacement_terms <- unique(.host_clean_term(failed_df$replacement_term))
  replacement_terms <- replacement_terms[!.host_is_invalid_tax_term(replacement_terms)]
  
  for (i in seq_len(nrow(failed_df))) {
    replacement <- .host_clean_term(failed_df$replacement_term[i])
    
    if (.host_is_invalid_tax_term(replacement)) {
      next
    }
    
    if (replacement %in% successful_taxonomy_terms) {
      failed_df$replacement_lookup_status[i] <- "replacement_successful"
      failed_df$lookup_status[i] <- "resolved_by_replacement"
    } else {
      failed_df$replacement_lookup_status[i] <- "replacement_failed"
    }
  }
  
  failed_replacements <- replacement_terms[
    !replacement_terms %in% successful_taxonomy_terms
  ]
  
  if (length(failed_replacements) > 0) {
    message("Replacement terms that still failed taxonomy lookup:")
    message("  - ", paste(failed_replacements, collapse = "\n  - "))
    
    failed_df <- .append_failed_terms(
      failed_df,
      terms = failed_replacements,
      term_type = "replacement_lookup",
      lookup_status = "failed"
    )
  }
  
  failed_df <- .add_failed_host_term_counts(failed_df, project_name, metadata_file)
  
  .write_failed_terms(failed_df, project_name, host_dir)
  
  merge_host_taxonomy_into_metadata(
    project_name = project_name,
    host_dir = host_dir,
    metadata_file = metadata_file
  )
  
  message("Refinement pass completed.")
  
  invisible(list(
    taxonomy = .read_host_taxonomy(project_name, host_dir),
    failed_table = failed_df,
    lookup_result = lookup_result
  ))
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




# final data merge and summary creation
run_host_assessment_summary <- function(
    project_name,
    fungal_rank = "genus",
    host_rank   = "phylum",
    keep_NAs    = FALSE,
    host_dir    = "./host_assessment"
) {
  
  summary <- summarize_host_usage(
    project_name = project_name,
    fungal_rank = fungal_rank,
    host_rank = host_rank,
    keep_NAs = keep_NAs,
    host_dir = host_dir
    # metadata_file intentionally NOT exposed
  )
  
  invisible(summary)
}




# ============================================================
# Phylogeny creation functions
# ============================================================

curate_metadata_regions <- function(project_name,
                                    mapping_file = NULL) {
  
  # Resolve default mapping file location
  if (is.null(mapping_file) || !nzchar(mapping_file)) {
    arborist_root <- get0("arborist_repo",
                          envir      = .GlobalEnv,
                          ifnotfound = normalizePath("~/github/aRborist", mustWork = FALSE))
    
    mapping_file <- file.path(arborist_root, "example_data", "region_replacement_patterns.csv")
  }
  
  # 0) read the input from the basic curation step
  infile  <- paste0("./metadata_files/all_accessions_pulled_metadata_",
                    project_name, "_curated.csv")
  acc_df  <- read.csv(infile, header = TRUE, stringsAsFactors = FALSE)
  
  # 1) set up the component columns we actually want to keep
  acc_df$gene.region.components      <- NA_character_
  acc_df$product.region.components   <- NA_character_
  acc_df$acc_title.region.components <- NA_character_
  
  # 2) load user/packaged mapping file, if file isn't found there are basic backup patterns
  if (file.exists(mapping_file)) {
    message("Using region replacement patterns from: ", mapping_file)
    map_df <- read.csv(mapping_file, stringsAsFactors = FALSE)
  } else {
    warning("Mapping file not found at: ", mapping_file,
            "\nFalling back to built-in defaults. To customize mappings, copy and edit ",
            "'example_data/region_replacement_patterns.csv' in your aRborist repo.")
    map_df <- data.frame(
      pattern = c(
        "tef-*\\d*", "EF1-alpha", "ef1a",
        "b-tub", "TUB2",
        "RBP2", "RPB2",
        "RBP1", "RPB1",
        "elongation factor 1",
        "internal transcribed",
        "26S", "28S",
        "small subunit ribosomal",
        "large[st]* subunit ribosomal",
        "beta-tubulin",
        "actin beta",
        "licensing\\D*7\\D*",
        "polymerase II larg[est]*",
        "polymerase II second largest",
        "large subunit ribosomal"
      ),
      standard = c(
        "TEF", "TEF", "TEF",
        "BTUB", "BTUB",
        "RPB2", "RPB2",
        "RPB1", "RPB1",
        "TEF",
        "ITS",
        "LSU", "LSU",
        "SSU",
        "LSU",
        "BTUB",
        "actin",
        "MCM7",
        "RPB1",
        "RPB2",
        "LSU"
      ),
      stringsAsFactors = FALSE
    )
  }
  
  # helper to append a new component to an existing semicolon list
  append_component <- function(current, add) {
    if (is.na(current) || current == "") {
      add
    } else {
      # avoid duplicates like ITS;ITS
      pattern <- paste0("(^|;)", add, "($|;)")
      if (grepl(pattern, current)) {
        current
      } else {
        paste(current, add, sep = ";")
      }
    }
  }
  
  # 3) CSV-based replacements: accumulate matches instead of overwriting
  for (i in seq_len(nrow(map_df))) {
    pat <- map_df$pattern[i]
    std <- map_df$standard[i]
    
    # gene column
    hit_gene <- stringr::str_detect(acc_df$gene,
                                    stringr::regex(pat, ignore_case = TRUE))
    if (any(hit_gene)) {
      current_vals <- acc_df$gene.region.components[hit_gene]
      acc_df$gene.region.components[hit_gene] <- mapply(
        append_component, current_vals, std, USE.NAMES = FALSE
      )
    }
    
    # product column
    hit_prod <- stringr::str_detect(acc_df$product,
                                    stringr::regex(pat, ignore_case = TRUE))
    if (any(hit_prod)) {
      current_vals <- acc_df$product.region.components[hit_prod]
      acc_df$product.region.components[hit_prod] <- mapply(
        append_component, current_vals, std, USE.NAMES = FALSE
      )
    }
    
    # accession title column
    hit_title <- stringr::str_detect(acc_df$accession_title,
                                     stringr::regex(pat, ignore_case = TRUE))
    if (any(hit_title)) {
      current_vals <- acc_df$acc_title.region.components[hit_title]
      acc_df$acc_title.region.components[hit_title] <- mapply(
        append_component, current_vals, std, USE.NAMES = FALSE
      )
    }
  }
  
  # 4) backup default component search (fills only blanks, super basic)
  detect_components <- function(txt) {
    if (is.null(txt) || is.na(txt) || txt == "") return(NA_character_)
    
    patterns <- list(
      ITS = stringr::regex("internal transcribed spacer|\\bITS\\b", ignore_case = TRUE),
      SSU = stringr::regex("18S|small subunit ribosomal", ignore_case = TRUE),
      LSU = stringr::regex("28S|large subunit ribosomal|26S", ignore_case = TRUE)
    )
    
    found <- character(0)
    for (nm in names(patterns)) {
      if (stringr::str_detect(txt, patterns[[nm]])) {
        found <- c(found, nm)
      }
    }
    found <- unique(found)
    if (length(found) == 0) {
      return(NA_character_)
    } else {
      paste(found, collapse = ";")
    }
  }
  
  for (row_i in seq_len(nrow(acc_df))) {
    # gene
    if (is.na(acc_df$gene.region.components[row_i]) || acc_df$gene.region.components[row_i] == "") {
      comp <- detect_components(acc_df$gene[row_i])
      if (!is.na(comp)) acc_df$gene.region.components[row_i] <- comp
    }
    
    # product
    if (is.na(acc_df$product.region.components[row_i]) || acc_df$product.region.components[row_i] == "") {
      comp <- detect_components(acc_df$product[row_i])
      if (!is.na(comp)) acc_df$product.region.components[row_i] <- comp
    }
    
    # accession title
    if (is.na(acc_df$acc_title.region.components[row_i]) || acc_df$acc_title.region.components[row_i] == "") {
      comp <- detect_components(acc_df$accession_title[row_i])
      if (!is.na(comp)) acc_df$acc_title.region.components[row_i] <- comp
    }
  }
  
  # 5) final region assignment from components, in priority order
  acc_df <- acc_df %>%
    dplyr::mutate(
      region.standard = dplyr::coalesce(
        gene.region.components,
        product.region.components,
        acc_title.region.components
      )
    )
  
  # 6) fasta headers
  acc_df$fasta.header      <- paste0(">", acc_df$org_name, "_", acc_df$strain.standard)
  acc_df$fasta.header.type <- paste0(">", acc_df$org_name, "_", acc_df$strain.standard.type)
  
  # 7) log unmatched rows
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
    
    unmatched_file <- paste0("./metadata_files/unmatched_regions_",
                             project_name, ".csv")
    write.csv(unmatched_df, unmatched_file, row.names = FALSE)
    
    message("Some records did not match any region pattern. ",
            "These were written to: ", unmatched_file)
  }
  
  # 8) write final curated metadata
  outfile <- paste0("./metadata_files/all_accessions_pulled_metadata_",
                    project_name, "_curated.csv")
  write.csv(acc_df, outfile, row.names = FALSE)
  cat("Wrote region-curated metadata to:", outfile, "\n")
}

# filtering metadata to only selected regions
select_regions <- function(project_name,
                           regions_to_include,
                           acc_to_exclude = character(0),
                           min_region_requirement = length(regions_to_include)) {
  
  if (!exists("base_dir", envir = .GlobalEnv)) {
    stop("`base_dir` is not defined. Run start_project() first.")
  }
  
  # Canonicalize regions + region-set name
  regions_to_include <- sort_regions(regions_to_include)
  region_set_name    <- paste(regions_to_include, collapse = ".")
  
  # 1) Define input path to curated metadata
  input_path <- file.path(
    base_dir,
    "metadata_files",
    paste0("all_accessions_pulled_metadata_", project_name, "_curated.csv")
  )
  if (!file.exists(input_path)) {
    stop("Curated metadata not found at: ", input_path)
  }
  
  accession_list <- read.csv(input_path, header = TRUE, stringsAsFactors = FALSE)
  
  # 2) Define region-set name and output directory
  phylo_dir <- file.path(base_dir, "phylogenies", region_set_name)
  if (!dir.exists(phylo_dir)) dir.create(phylo_dir, recursive = TRUE)
  
  # 3) Optional exclusion of specific accessions
  if (!is.null(acc_to_exclude) && length(acc_to_exclude) > 0 && any(acc_to_exclude != "")) {
    accession_list <- accession_list[!accession_list$Accession %in% acc_to_exclude, ]
  }
  
  # 4) Filter: keep only desired regions
  multifasta_prep <- accession_list[accession_list$region.standard %in% regions_to_include, ]
  
  # 5) Write "long" filtered file (per accession)
  output_long <- file.path(
    phylo_dir,
    paste0("selected_accessions_metadata_", project_name, ".", region_set_name, ".csv")
  )
  write.csv(multifasta_prep, output_long, row.names = FALSE)
  
  # 6) Create "wide" region attendance table
  multifasta_prep_complete <- subset(
    multifasta_prep,
    select = c(strain.standard.type, organism, Accession, region.standard)
  )
  
  multifasta_prep_select <- dplyr::distinct(
    multifasta_prep_complete,
    strain.standard.type, region.standard,
    .keep_all = TRUE
  )
  
  select_region_attendance <- tidyr::pivot_wider(
    multifasta_prep_select,
    names_from = "region.standard",
    values_from = "Accession"
  )
  
  # 7) Apply minimum region inclusion threshold
  select_region_attendance_filtered <- select_region_attendance %>%
    dplyr::mutate(
      total = rowSums(!is.na(dplyr::select(., tidyselect::any_of(regions_to_include))))
    ) %>%
    dplyr::filter(total >= min_region_requirement) %>%
    dplyr::select(-total)
  
  output_wide <- file.path(
    phylo_dir,
    paste0("Region_attendance_sheet_", project_name, ".", region_set_name, ".csv")
  )
  write.csv(select_region_attendance_filtered, output_wide, row.names = FALSE)
  
  # 8) Messages for user
  message("✓ Filtered metadata written to: ", output_long)
  message("✓ Region attendance sheet written to: ", output_wide)
  message("✓ Region set: ", region_set_name)
}

create_multifastas <- function(project_name,
                               regions_to_include) {
  if (!exists("base_dir", envir = .GlobalEnv)) {
    stop("`base_dir` is not defined. Run start_project() first.")
  }
  
  # Canonicalize regions + region-set name
  regions_to_include <- sort_regions(regions_to_include)
  region_set_name    <- paste(regions_to_include, collapse = ".")
  
  # 0) get region set and folders
  phylo_dir <- file.path(base_dir, "phylogenies", region_set_name)
  if (!dir.exists(phylo_dir)) {
    stop("Expected region-set folder not found: ", phylo_dir,
         "\nDid you run select_regions() for this region set?")
  }
  
  prep_dir <- file.path(phylo_dir, "prep")
  if (!dir.exists(prep_dir)) dir.create(prep_dir, recursive = TRUE)
  
  # One subfolder per region (downstream: raw, aligned, trimmed files live here)
  for (rg in regions_to_include) {
    rg_dir <- file.path(prep_dir, rg)
    if (!dir.exists(rg_dir)) dir.create(rg_dir, recursive = TRUE)
  }
  
  # 1) Load inputs created by select_regions()
  attendance_path <- file.path(
    phylo_dir,
    paste0("Region_attendance_sheet_", project_name, ".", region_set_name, ".csv")
  )
  long_filtered_path <- file.path(
    phylo_dir,
    paste0("selected_accessions_metadata_", project_name, ".", region_set_name, ".csv")
  )
  
  if (!file.exists(attendance_path)) {
    stop("Region attendance sheet not found: ", attendance_path)
  }
  if (!file.exists(long_filtered_path)) {
    stop("Filtered metadata (long) not found: ", long_filtered_path)
  }
  
  region_attendance <- read.csv(attendance_path, header = TRUE, stringsAsFactors = FALSE)
  filtered_long <- read.csv(long_filtered_path, header = TRUE, stringsAsFactors = FALSE)
  
  # Sanity checks on needed columns
  needed_cols <- c("Accession", "region.standard", "fasta.header.type", "sequence")
  missing_cols <- setdiff(needed_cols, colnames(filtered_long))
  if (length(missing_cols) > 0) {
    stop("Missing columns in filtered metadata: ", paste(missing_cols, collapse = ", "),
         "\nUpstream curation must provide these.")
  }
  
  # 2) Build per-region accession vectors from the wide attendance sheet
  region_cols <- intersect(regions_to_include, colnames(region_attendance))
  if (length(region_cols) == 0) {
    stop("None of the requested regions are present as columns in the attendance sheet.")
  }
  
  # Collect per-region accession IDs (drop NAs), preserving strain/region mapping done upstream
  region_accessions <- lapply(region_cols, function(rg) {
    unique(na.omit(region_attendance[[rg]]))
  })
  names(region_accessions) <- region_cols
  
  # 3) For each region, subset sequences and write a multifasta into prep/<region>/
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
    
    sub_df <- filtered_long[filtered_long$Accession %in% acc_vec &
                              filtered_long$region.standard == rg, ]
    
    # Drop rows with missing sequences
    sub_df <- sub_df[!is.na(sub_df$sequence) & sub_df$sequence != "", ]
    
    # Order by header (stable, human-friendly)
    sub_df <- sub_df[order(sub_df$fasta.header.type, decreasing = FALSE), ]
    
    # Ensure headers begin with '>'
    headers <- sub_df$fasta.header.type
    needs_gt <- !startsWith(headers, ">")
    headers[needs_gt] <- paste0(">", headers[needs_gt])
    
    # Interleave header/sequence to write a simple FASTA
    seqs_fasta <- c(rbind(headers, sub_df$sequence))
    
    rg_dir <- file.path(prep_dir, rg)
    fasta_name <- paste0(project_name, ".", region_set_name, "_", rg, ".raw.fasta")
    fasta_path <- file.path(rg_dir, fasta_name)
    
    writeLines(seqs_fasta, con = fasta_path)
    message("Created multifasta for region ", rg, ": ", fasta_path)
    
    manifest <- rbind(
      manifest,
      data.frame(region = rg,
                 n_sequences = nrow(sub_df),
                 fasta_path = fasta_path,
                 stringsAsFactors = FALSE)
    )
  }
  
  # 4) Write a small manifest for bookkeeping/QC
  manifest_path <- file.path(prep_dir, paste0("multifasta_manifest_", project_name, ".", region_set_name, ".tsv"))
  write.table(manifest, manifest_path, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote manifest: ", manifest_path)
  
  invisible(manifest)
}

align_regions_mafft <- function(project_name,
                                regions_to_include,
                                threads = max(1, parallel::detectCores() - 1),
                                mafft_args = c("--auto", "--reorder"),
                                force = FALSE) {
  if (!exists("base_dir", envir = .GlobalEnv)) {
    stop("`base_dir` is not defined. Run start_project() first.")
  }
  
  # Canonicalize regions + region-set name
  regions_to_include <- sort_regions(regions_to_include)
  region_set_name    <- paste(regions_to_include, collapse = ".")
  
  # 0) Determine MAFFT executable path
  mafft_path <- Sys.getenv("MAFFT_PATH", unset = "mafft")
  
  # Verify MAFFT works before starting
  check_result <- suppressWarnings(system2(mafft_path, "--version", stdout = TRUE, stderr = TRUE))
  if (length(check_result) == 0 || grepl("not found|No such file", check_result[1], ignore.case = TRUE)) {
    stop(
      "MAFFT not found. Please install it or set MAFFT_PATH in your .Renviron file.\n",
      "Example:  MAFFT_PATH=/usr/local/bin/mafft\n",
      "Then restart R and rerun this command."
    )
  } else {
    message("Using MAFFT executable: ", mafft_path)
  }
  
  # 1) Folder setup
  phylo_dir <- file.path(base_dir, "phylogenies", region_set_name)
  prep_dir  <- file.path(phylo_dir, "prep")
  if (!dir.exists(prep_dir)) {
    stop("Prep directory not found: ", prep_dir,
         "\nDid you run create_multifastas() for this region set?")
  }
  
  manifest <- data.frame(
    region     = character(0),
    raw_fasta  = character(0),
    aligned_fasta = character(0),
    log_path   = character(0),
    status     = character(0),
    stringsAsFactors = FALSE
  )
  
  # 2) Loop through regions
  for (rg in regions_to_include) {
    rg_dir <- file.path(prep_dir, rg)
    if (!dir.exists(rg_dir)) {
      warning("Region prep folder missing (skipping): ", rg_dir)
      next
    }
    
    raw_fa <- file.path(rg_dir, paste0(project_name, ".", region_set_name, "_", rg, ".raw.fasta"))
    aln_fa <- file.path(rg_dir, paste0(project_name, ".", region_set_name, "_", rg, ".aligned.fasta"))
    log_fp <- file.path(rg_dir, paste0(project_name, ".", region_set_name, "_", rg, ".mafft.log"))
    
    if (!file.exists(raw_fa)) {
      warning("Raw FASTA not found for region ", rg, ": ", raw_fa)
      manifest <- rbind(manifest, data.frame(
        region = rg, raw_fasta = raw_fa, aligned_fasta = NA, log_path = log_fp, status = "missing_raw",
        stringsAsFactors = FALSE
      ))
      next
    }
    
    if (file.exists(aln_fa) && !force) {
      message("Aligned FASTA already exists (use force=TRUE to overwrite): ", aln_fa)
      manifest <- rbind(manifest, data.frame(
        region = rg, raw_fasta = raw_fa, aligned_fasta = aln_fa, log_path = log_fp, status = "skipped_exists",
        stringsAsFactors = FALSE
      ))
      next
    }
    
    # 3) Run MAFFT alignment
    message("Running MAFFT for region ", rg, " …")
    
    mafft_args_full <- c("--thread", as.character(threads), mafft_args, raw_fa)
    
    exit_code <- tryCatch({
      system2(command = mafft_path,
              args    = mafft_args_full,
              stdout  = aln_fa,
              stderr  = log_fp)
    }, error = function(e) {
      warning("MAFFT invocation failed for ", rg, ": ", conditionMessage(e))
      return(1L)
    })
    
    status <- if (!is.null(exit_code) && exit_code == 0L && file.exists(aln_fa)) "ok" else "failed"
    
    manifest <- rbind(manifest, data.frame(
      region = rg,
      raw_fasta = raw_fa,
      aligned_fasta = if (file.exists(aln_fa)) aln_fa else NA,
      log_path = log_fp,
      status = status,
      stringsAsFactors = FALSE
    ))
    
    if (status != "ok") {
      warning("MAFFT failed for region ", rg, ". See log: ", log_fp)
    } else {
      message("✓ Aligned FASTA written: ", aln_fa)
    }
  }
  
  # 4) Write alignment manifest
  align_manifest <- file.path(prep_dir, paste0("alignment_manifest_", project_name, ".", region_set_name, ".tsv"))
  write.table(manifest, align_manifest, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Alignment manifest: ", align_manifest)
  
  invisible(manifest)
}

trim_regions_trimal <- function(project_name,
                                regions_to_include,
                                trimal_args = c("--automated1"),
                                force = FALSE) {
  if (!exists("base_dir", envir = .GlobalEnv)) {
    stop("`base_dir` is not defined. Run start_project() first.")
  }
  
  # Canonicalize regions + region-set name
  regions_to_include <- sort_regions(regions_to_include)
  region_set_name    <- paste(regions_to_include, collapse = ".")
  
  # 0) Determine trimAl path
  trimal_path <- Sys.getenv("TRIMAL_PATH", unset = "trimal")
  
  # Verify that trimAl runs
  check_result <- suppressWarnings(system2(trimal_path, "--version", stdout = TRUE, stderr = TRUE))
  if (length(check_result) == 0 || grepl("not found|No such file", check_result[1], ignore.case = TRUE)) {
    stop(
      "trimAl not found. Please install it or set TRIMAL_PATH in your .Renviron file.\n",
      "Example:  TRIMAL_PATH=/usr/local/bin/trimal\n",
      "Then restart R and rerun this command."
    )
  } else {
    message("Using trimAl executable: ", trimal_path)
  }
  
  # 1) Folder setup
  phylo_dir <- file.path(base_dir, "phylogenies", region_set_name)
  prep_dir  <- file.path(phylo_dir, "prep")
  if (!dir.exists(prep_dir)) {
    stop("Prep directory not found: ", prep_dir,
         "\nDid you run align_regions_mafft() for this region set?")
  }
  
  manifest <- data.frame(
    region     = character(0),
    aligned_fasta = character(0),
    trimmed_fasta = character(0),
    log_path   = character(0),
    status     = character(0),
    stringsAsFactors = FALSE
  )
  
  # 2) Loop through regions
  for (rg in regions_to_include) {
    rg_dir <- file.path(prep_dir, rg)
    aln_fa <- file.path(rg_dir, paste0(project_name, ".", region_set_name, "_", rg, ".aligned.fasta"))
    trimmed_fa <- file.path(rg_dir, paste0(project_name, ".", region_set_name, "_", rg, ".trimmed.fasta"))
    log_fp <- file.path(rg_dir, paste0(project_name, ".", region_set_name, "_", rg, ".trimal.log"))
    
    if (!file.exists(aln_fa)) {
      warning("Aligned FASTA not found for region ", rg, ": ", aln_fa)
      next
    }
    
    if (file.exists(trimmed_fa) && !force) {
      message("Trimmed FASTA already exists (use force=TRUE to overwrite): ", trimmed_fa)
      manifest <- rbind(manifest, data.frame(
        region = rg,
        aligned_fasta = aln_fa,
        trimmed_fasta = trimmed_fa,
        log_path = log_fp,
        status = "skipped_exists",
        stringsAsFactors = FALSE
      ))
      next
    }
    
    # 3) Run trimAl
    message("Running trimAl for region ", rg, " …")
    
    trimal_args_full <- c(trimal_args, "-in", aln_fa, "-out", trimmed_fa)
    exit_code <- tryCatch({
      system2(command = trimal_path,
              args    = trimal_args_full,
              stdout  = log_fp,
              stderr  = log_fp)
    }, error = function(e) {
      warning("trimAl invocation failed for ", rg, ": ", conditionMessage(e))
      return(1L)
    })
    
    status <- if (!is.null(exit_code) && exit_code == 0L && file.exists(trimmed_fa)) "ok" else "failed"
    
    manifest <- rbind(manifest, data.frame(
      region = rg,
      aligned_fasta = aln_fa,
      trimmed_fasta = if (file.exists(trimmed_fa)) trimmed_fa else NA,
      log_path = log_fp,
      status = status,
      stringsAsFactors = FALSE
    ))
    
    if (status != "ok") {
      warning("trimAl failed for region ", rg, ". See log: ", log_fp)
    } else {
      message("✓ Trimmed FASTA written: ", trimmed_fa)
    }
  }
  
  # 4) Write manifest
  trim_manifest <- file.path(prep_dir, paste0("trim_manifest_", project_name, ".", region_set_name, ".tsv"))
  write.table(manifest, trim_manifest, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Trim manifest: ", trim_manifest)
  
  invisible(manifest)
}

# running IQTREE modelfinder step for single region 
iqtree_modelfinder_per_region <- function(project_name,
                                          regions_to_include,
                                          threads    = max(1, parallel::detectCores() - 1),
                                          iqtree_args = c("-m", "MFP+MERGE", "-nt", "AUTO", "-quiet"),
                                          single_gene_bootstraps = 1000,
                                          force = FALSE) {
  
  if (!exists("base_dir", envir = .GlobalEnv)) {
    stop("`base_dir` is not defined. Run start_project() first.")
  }
  
  # 0) IQ-TREE software check
  iqtree_bin <- Sys.getenv("IQTREE_PATH")
  if (!nzchar(iqtree_bin)) {
    stop("IQTREE_PATH is not set in .Renviron.")
  }
  suppressWarnings(system2(iqtree_bin, "-version"))
  
  # Canonicalize regions + region-set name
  regions_to_include <- sort_regions(regions_to_include)
  region_set_name    <- paste(regions_to_include, collapse = ".")
  
  # 1) Paths
  phylo_dir <- file.path(base_dir, "phylogenies", region_set_name)
  prep_dir  <- file.path(phylo_dir, "prep")
  single_gene_dir <- file.path(phylo_dir, "single_gene_trees")
  multi_gene_dir  <- file.path(phylo_dir, "multi_gene_trees")
  
  if (!dir.exists(prep_dir)) stop("Prep directory not found: ", prep_dir)
  if (!dir.exists(single_gene_dir)) dir.create(single_gene_dir, recursive = TRUE)
  if (!dir.exists(multi_gene_dir)) dir.create(multi_gene_dir, recursive = TRUE)
  
  results <- list()
  
  # 2) Loop over regions
  for (rg in regions_to_include) {
    rg_prep_dir <- file.path(prep_dir, rg)
    trimmed_fa <- file.path(
      rg_prep_dir,
      paste0(project_name, ".", region_set_name, "_", rg, ".trimmed.fasta")
    )
    
    if (!file.exists(trimmed_fa)) {
      warning("Missing trimmed alignment for region ", rg)
      next
    }
    
    # Output location for single-gene IQ-TREE result
    rg_sg_dir <- file.path(single_gene_dir, rg)
    if (!dir.exists(rg_sg_dir)) dir.create(rg_sg_dir, recursive = TRUE)
    
    prefix_base <- paste0(project_name, ".", region_set_name, "_", rg, ".modeltest")
    prefix      <- file.path(rg_sg_dir, prefix_base)
    iqtreefile  <- paste0(prefix, ".iqtree")
    
    # 3) Run IQ-TREE (ModelFinder + bootstraps)
    if (!file.exists(iqtreefile) || force) {
      
      args <- c(
        "-s", trimmed_fa,
        "-pre", prefix,
        "-bb", as.character(single_gene_bootstraps)
      )
      
      # Add -nt <threads> unless user set -nt in iqtree_args
      if (!any(iqtree_args == "-nt")) {
        args <- c(args, "-nt", as.character(threads))
      }
      
      # Add user-specified options
      args <- c(args, iqtree_args)
      
      message("Running ModelFinder with UF bootstraps for region ", rg, " …")
      system2(iqtree_bin, args, stdout = TRUE, stderr = TRUE)
    }
    
    # 4) Extract best model
    iqtxt <- readLines(iqtreefile, warn = FALSE)
    model_idx <- grep("Best-fit model according to BIC:", iqtxt, fixed = TRUE)
    
    if (length(model_idx) == 0L) {
      best_model <- NA_character_
      warning("Could not find model line in ", iqtreefile)
    } else {
      best_line <- iqtxt[model_idx[1]]
      best_model <- stringr::str_trim(sub(".*Best-fit model according to BIC:\\s*", "",
                                          best_line))
    }
    
    # alignment length
    aln <- Biostrings::readDNAStringSet(trimmed_fa)
    aln_length <- unique(Biostrings::width(aln))[1]
    
    results[[rg]] <- data.frame(
      region      = rg,
      best_model  = best_model,
      aln_length  = aln_length,
      iqtree_file = iqtreefile,
      stringsAsFactors = FALSE
    )
  }
  
  # 5) Output TSV to multi_gene_trees folder
  model_fits <- dplyr::bind_rows(results)
  
  out_tsv <- file.path(
    multi_gene_dir,
    paste0("model_fits_", project_name, ".", region_set_name, ".tsv")
  )
  
  readr::write_tsv(model_fits, out_tsv)
  message("ModelFinder summary written to: ", out_tsv)
  
  invisible(model_fits)
}

# Create Projects/<project_name>/ and run your normal structure inside it.
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

concatenate_and_write_partitions <- function(project_name,
                                             regions_to_include) {
  if (!exists("base_dir", envir = .GlobalEnv)) {
    stop("`base_dir` is not defined. Run start_project() first.")
  }
  
  # Canonicalize regions + region-set name
  regions_to_include <- sort_regions(regions_to_include)
  region_set_name    <- paste(regions_to_include, collapse = ".")
  
  # 0) Paths and region-set name
  phylo_dir       <- file.path(base_dir, "phylogenies", region_set_name)
  prep_dir        <- file.path(phylo_dir, "prep")
  multi_gene_dir  <- file.path(phylo_dir, "multi_gene_trees")
  
  if (!dir.exists(prep_dir)) {
    stop("Prep directory not found: ", prep_dir,
         "\nDid you run trim_regions_trimal() for this region set?")
  }
  if (!dir.exists(multi_gene_dir)) dir.create(multi_gene_dir, recursive = TRUE)
  
  # 1) Read model_fits TSV (from previous step)
  model_fits_path <- file.path(
    multi_gene_dir,
    paste0("model_fits_", project_name, ".", region_set_name, ".tsv")
  )
  if (!file.exists(model_fits_path)) {
    stop("model_fits TSV not found: ", model_fits_path,
         "\nDid you run iqtree_modelfinder_per_region() ?")
  }
  
  model_fits <- readr::read_tsv(model_fits_path, show_col_types = FALSE)
  
  model_lookup <- model_fits |>
    dplyr::select(region, best_model)
  
  if (!all(regions_to_include %in% model_lookup$region)) {
    warning("Some regions in regions_to_include are missing from model_fits TSV.")
  }
  
  # 2) Read trimmed alignments for each region
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
  
  # 3) Build union of taxa and concatenate sequences in the specified order
  all_taxa <- sort(unique(unlist(lapply(aln_per_region, names))))
  
  concat_vec <- vapply(all_taxa, function(taxon) {
    paste0(
      vapply(regions_to_include, function(rg) {
        s <- aln_per_region[[rg]]
        if (taxon %in% names(s)) {
          as.character(s[[taxon]])
        } else {
          width_rg <- Biostrings::width(s)[1]
          paste(rep("-", width_rg), collapse = "")
        }
      }, FUN.VALUE = character(1)),
      collapse = ""
    )
  }, FUN.VALUE = character(1))
  
  concat_dna <- Biostrings::DNAStringSet(concat_vec)
  names(concat_dna) <- all_taxa
  
  # 4) Compute partition coordinates (1-based, inclusive)
  region_lengths <- vapply(regions_to_include, function(rg) {
    unique(Biostrings::width(aln_per_region[[rg]]))[1]
  }, FUN.VALUE = integer(1))
  
  starts <- cumsum(c(1, head(region_lengths, -1)))
  ends   <- cumsum(region_lengths)
  
  if (unique(Biostrings::width(concat_dna))[1] != tail(ends, 1)) {
    warning("Concatenated alignment length does not match sum of region lengths.")
  }
  
  # 5) Write concatenated supermatrix FASTA
  concat_path <- file.path(
    multi_gene_dir,
    paste0("concatenated_", project_name, ".", region_set_name, ".fasta")
  )
  Biostrings::writeXStringSet(concat_dna, filepath = concat_path, format = "fasta")
  message("Concatenated alignment written to: ", concat_path)
  
  # 6) Build NEXUS partition file with correct models per region
  models_ordered <- vapply(regions_to_include, function(rg) {
    row <- model_lookup[model_lookup$region == rg, , drop = FALSE]
    
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
  }, FUN.VALUE = character(1))
  
  partition_lines <- c("#nexus",
                       "begin sets;")
  
  part_names <- paste0("part", seq_along(regions_to_include))
  
  for (i in seq_along(regions_to_include)) {
    line <- sprintf("\tcharset %s = %d-%d;",
                    part_names[i], starts[i], ends[i])
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
  
  invisible(list(
    concat_fasta    = concat_path,
    partitions_nex  = part_path,
    regions         = regions_to_include,
    starts          = starts,
    ends            = ends,
    models          = models_ordered
  ))
}

iqtree_multigene_partitioned <- function(
    project_name,
    regions_to_include,
    threads = 8,
    multigene_bootstraps = 1000,
    iqtree_args = c("-m", "MFP+MERGE"),
    force = TRUE
) {
  if (!exists("base_dir", envir = .GlobalEnv)) {
    stop("`base_dir` is not defined. Run start_project() first.")
  }
  
  # 0) IQ-TREE software check
  iqtree_bin <- Sys.getenv("IQTREE_PATH")
  if (!nzchar(iqtree_bin)) {
    stop("IQTREE_PATH is not set in .Renviron.")
  }
  suppressWarnings(system2(iqtree_bin, "-version"))
  
  # Canonicalize regions + region-set name
  regions_to_include <- sort_regions(regions_to_include)
  region_set_name    <- paste(regions_to_include, collapse = ".")
  
  # 1) Core paths
  phylo_dir <- file.path(base_dir, "phylogenies", region_set_name)
  if (!dir.exists(phylo_dir)) {
    stop("Phylogenies directory does not exist for region set: ", phylo_dir)
  }
  
  multi_gene_dir <- file.path(phylo_dir, "multi_gene_trees")
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
  
  # 2) Basic checks
  if (!file.exists(concatenated_fasta)) {
    stop("Concatenated alignment not found: ", concatenated_fasta)
  }
  if (!file.exists(partition_nexus)) {
    stop("Partition Nexus file not found: ", partition_nexus)
  }
  
  if (file.exists(expected_treefile) && !force) {
    stop(
      "IQ-TREE treefile already exists and force = FALSE:\n  ",
      expected_treefile, "\n",
      "Set force = TRUE to overwrite or move/rename the existing files."
    )
  }
  
  # 3) Build IQ-TREE command
  args <- c(
    "-s", concatenated_fasta,
    "-p", partition_nexus,
    "--ufboot", as.character(multigene_bootstraps)
  )
  
  if (!any(iqtree_args == "-nt")) {
    args <- c(args, "-nt", as.character(threads))
  }
  
  args <- c(args, iqtree_args, "--prefix", iqtree_prefix)
  
  message("Running IQ-TREE with command:\n",
          iqtree_bin, " ", paste(shQuote(args), collapse = " "))
  
  # 4) Run IQ-TREE
  exit_status <- system2(command = iqtree_bin, args = args)
  
  if (exit_status != 0) {
    warning("IQ-TREE finished with non-zero exit status: ", exit_status)
  } else {
    message("IQ-TREE multigene partitioned run completed successfully.")
    if (file.exists(expected_treefile)) {
      message("Treefile: ", expected_treefile)
    }
  }
  
  invisible(list(
    status        = exit_status,
    cmd           = paste(iqtree_bin, paste(args, collapse = " ")),
    output_prefix = iqtree_prefix,
    treefile      = expected_treefile
  ))
}





# ============================================================
# arborist extra wrappers
# ============================================================

# aquiring accessions + data from NCBI
ncbi_data_fetch <- function(taxa_list,
                            max_acc_per_taxa = "max",
                            organism_scope   = NULL,
                            include_filters  = NULL,
                            exclude_filters  = NULL,
                            project_name     = get0("project_name", envir = .GlobalEnv)) {
  
  if (is.null(project_name) || !nzchar(project_name)) {
    stop("project_name is not set. Run start_project(project_name) first, or pass project_name explicitly.")
  }
  
  # 1) Fetch accessions
  accession_list <- get_accessions_for_all_taxa(
    taxa_list        = taxa_list,
    max_acc_per_taxa = max_acc_per_taxa,
    organism_scope   = organism_scope,
    include_filters  = include_filters,
    exclude_filters  = exclude_filters
  )
  
  # 2) Fetch metadata for those accessions
  retrieve_ncbi_metadata(project_name)
  
  invisible(accession_list)
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


