# ============================================================
# aRborist helper functions
# ============================================================

# Package setup
required_packages <- c(
  "rentrez","stringr","plyr","dplyr","withr","XML",
  "data.table","tidyr","phylotools","scales",
  "purrr","readr","phytools","RColorBrewer",
  "maps","ggplot2","tidygeocoder",
  "ggrepel","taxize","Biostrings", "yaml"
)

installed_packages <- required_packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(required_packages[!installed_packages])
}

load_required_packages <- function() {
  invisible(lapply(required_packages, library, character.only = TRUE))
}


# Defaults for entrez search term (need to change later so users can override in run script)
if (!exists("organism_scope", .GlobalEnv) ||
    is.null(organism_scope) ||
    !nzchar(organism_scope)) {
  assign("organism_scope", "txid4751[Organism:exp]", .GlobalEnv)
}

if (!exists("search_options", .GlobalEnv) ||
    is.null(search_options) ||
    !nzchar(search_options)) {
  assign("search_options",
         "(biomol_genomic[PROP]) AND (100[SLEN]:5000[SLEN]) AND is_nuccore[filter] AND NOT mitochondrion[filter]",
         .GlobalEnv)
}



compose_entrez_term <- function(taxon,
                                organism_scope = NULL,
                                extra_filters  = NULL) {
  q <- sprintf('("%s"[Organism])', taxon)

  if (!is.null(organism_scope) && nzchar(organism_scope)) {
    q <- paste(q, organism_scope, sep = " AND ")
  }
  if (!is.null(extra_filters) && nzchar(extra_filters)) {
    q <- paste(q, extra_filters, sep = " AND ")
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

# Project structure
setup_project_structure <- function(project_dir,
                                    subdirs = c("intermediate_files", "metadata_files", "results_files",
                                                "temp_files", "multifastas", "multifastas/aligned_fastas",
                                                "multigene_tree", "multigene_tree/prep", "single_gene_trees")) {
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

# Fetch accessions from NCBI
fetch_accessions_for_taxon <- function(taxon,
                                       max_acc = max_acc_per_taxa,
                                       organism_scope = NULL,
                                       extra_filters  = NULL) {
  cat("Searching term:", taxon, "\n")

  if (is.null(organism_scope) && exists("organism_scope", .GlobalEnv))
    organism_scope <- get("organism_scope", .GlobalEnv)
  if (is.null(extra_filters) && exists("search_options", .GlobalEnv))
    extra_filters <- get("search_options", .GlobalEnv)

  filters <- compose_entrez_term(taxon, organism_scope, extra_filters)

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

    tmp <- paste0("./temp_files/temp_file_accessions_from_", taxon, ".txt")
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


get_accessions_for_all_taxa <- function(taxa_list,
                                        max_acc_per_taxa,
                                        organism_scope = NULL,
                                        extra_filters  = NULL,
                                        timing_file = "./intermediate_files/fetch_times_accessions.csv") {
  # Resolve defaults if not provided
  if (is.null(organism_scope) && exists("organism_scope", .GlobalEnv))
    organism_scope <- get("organism_scope", envir = .GlobalEnv)
  if (is.null(extra_filters) && exists("search_options", .GlobalEnv))
    extra_filters <- get("search_options", envir = .GlobalEnv)

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
        taxon = term,
        max_acc = max_acc_per_taxa,
        organism_scope = organism_scope,
        extra_filters = extra_filters
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
    # save per-taxon accessions (even if empty, so users see it was checked)
    outfile_name <- paste0("./intermediate_files/Accessions_for_", term, ".csv")
    write.csv(tempdf, outfile_name, row.names = FALSE, quote = FALSE)
    taxa_frame_acc[[i]] <- tempdf
  }

  total_elapsed <- as.numeric(difftime(Sys.time(), overall_start, units = "mins"))
  cat("\nAll taxa completed in", round(total_elapsed, 2), "minutes.\n")

  # bind all results (handle case where some/all are empty)
  non_empty <- Filter(function(x) nrow(x) > 0, taxa_frame_acc)
  accession_list <- if (length(non_empty)) {
    do.call(rbind, non_empty)
  } else {
    data.frame(Accession = character(0), genus = character(0), stringsAsFactors = FALSE)
  }
# write log with fetch times
  write.csv(accession_list, "./intermediate_files/all_pulled_accessions.csv", row.names = FALSE)
  write.csv(timing_log, timing_file, row.names = FALSE)
  cat("Timing log written to ", timing_file, "\n", sep = "")

  return(accession_list)
}


# Metadata retrieval
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
retrieve_ncbi_metadata <- function(project_name) {
  accession_list <- read.csv("./intermediate_files/all_pulled_accessions.csv", header = TRUE)

  # Split by taxon if available, otherwise treat as a single group ("ALL")
  if ("genus" %in% names(accession_list)) {
    taxa_groups <- split(accession_list, accession_list$genus)
  } else {
    taxa_groups <- list(ALL = accession_list)
  }

  metadata_database_list <- vector("list", length = nrow(accession_list))
  fill_idx <- 1L

  timing_log <- data.frame(
    Taxon = character(),
    Num_accessions = integer(),
    Start_time = character(),
    End_time = character(),
    Elapsed_minutes = numeric(),
    stringsAsFactors = FALSE
  )

  overall_start <- Sys.time()
  cat("\nStarting metadata retrieval for", nrow(accession_list), "accessions across",
      length(taxa_groups), "taxon group(s)...\n")

  # Loop over taxa (genus) blocks,time each block
  for (tx in names(taxa_groups)) {
    block <- taxa_groups[[tx]]
    taxon_start <- Sys.time()
    cat(sprintf("\n--- %s: %d accession(s) ---\n", tx, nrow(block)))

    # Process each accession in this taxon
for (i in seq_len(nrow(block))) {
  acc <- block$Accession[i]
  cat(sprintf("[%d/%d] %s ... ", i, nrow(block), acc))
  tryCatch({
    entry <- fetch_metadata_for_accession(acc)
    metadata_database_list[[fill_idx]] <- entry
    fill_idx <- fill_idx + 1L

    # metadata summary printout
    cat("#", i,
        "| Accession:", entry$Accession,
        "| Species:", entry$organism,
        "| Strain:", dplyr::coalesce(entry$strain, entry$specimen_voucher, entry$isolate, ""),
        "| Isolation source:", entry$isolation_source,
        "| Host:", entry$host, "\n")

    cat("done\n")

  }, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
  })
  Sys.sleep(get_sleep_duration())
}


    taxon_end <- Sys.time()
    taxon_elapsed <- as.numeric(difftime(taxon_end, taxon_start, units = "mins"))
    timing_log <- rbind(
      timing_log,
      data.frame(
        Taxon = tx,
        Num_accessions = nrow(block),
        Start_time = format(taxon_start, "%Y-%m-%d %H:%M:%S"),
        End_time = format(taxon_end, "%Y-%m-%d %H:%M:%S"),
        Elapsed_minutes = round(taxon_elapsed, 2),
        stringsAsFactors = FALSE
      )
    )
    cat(sprintf("%s completed in %.2f minutes\n", tx, taxon_elapsed))
  }

  overall_end <- Sys.time()
  total_elapsed <- as.numeric(difftime(overall_end, overall_start, units = "mins"))
  cat("\nAll metadata retrieved in", round(total_elapsed, 2), "minutes total.\n")

  # Bind all non-NULL entries
  metadata_database <- plyr::rbind.fill(Filter(Negate(is.null), metadata_database_list))
  write.csv(
    metadata_database,
    paste0("./metadata_files/all_accessions_pulled_metadata_", project_name, ".csv"),
    row.names = FALSE
  )
  cat("Metadata saved for project:", project_name, "\n")

  # Add a TOTAL row and write timing CSV
  timing_log <- rbind(
    timing_log,
    data.frame(
      Taxon = "TOTAL",
      Num_accessions = nrow(accession_list),
      Start_time = format(overall_start, "%Y-%m-%d %H:%M:%S"),
      End_time = format(overall_end, "%Y-%m-%d %H:%M:%S"),
      Elapsed_minutes = round(total_elapsed, 2),
      stringsAsFactors = FALSE
    )
  )
  write.csv(timing_log, "./intermediate_files/fetch_times_metadata_by_taxon.csv", row.names = FALSE)
  cat("Timing log written to ./intermediate_files/fetch_times_metadata_by_taxon.csv\n")
}



# Custom sequences merge
merge_metadata_with_custom_file <- function(project_name) {
  metadata_file_path <- paste0("./metadata_files/all_accessions_pulled_metadata_", project_name, ".csv")
  metadata_database <- read.csv(metadata_file_path, header = TRUE)

  # if path is empty or missing, skip merge step
  if (is.null(my_lab_sequences) || !nzchar(my_lab_sequences)) {
    cat("No custom sequences file provided; skipping merge.\n")
    return(invisible(NULL))
  }

  custom_sequences <- read.csv(my_lab_sequences, header = TRUE, fill = TRUE)

  if (!all(c("Accession","strain","sequence","organism","gene") %in% colnames(custom_sequences))) {
    stop("Custom file must contain at least 'Accession', 'strain', 'sequence', 'organism', and 'gene' columns. If your sequences do not have accessions, you can simply use their strain name or other unique identifier.")
  }

  merged_data <- plyr::rbind.fill(metadata_database, custom_sequences)
  write.csv(merged_data, metadata_file_path, row.names = FALSE)
  cat("Merged custom sequences into:", metadata_file_path, "\n")
}


# Metadata curation and region selection
curate_metadata_basic <- function(project_name,
                                  taxa_of_interest = NULL) {
  metadata_file_path <- paste0("./metadata_files/all_accessions_pulled_metadata_", project_name, ".csv")
  accession_list <- read.csv(metadata_file_path, header = TRUE, stringsAsFactors = FALSE)

  if (!"specimen_voucher" %in% names(accession_list)) accession_list$specimen_voucher <- NA_character_
  if (!"strain" %in% names(accession_list))           accession_list$strain           <- NA_character_
  if (!"isolate" %in% names(accession_list))          accession_list$isolate          <- NA_character_
  if (!"type_material" %in% names(accession_list))    accession_list$type_material    <- NA_character_
  if (!"geo_loc_name" %in% names(accession_list))     accession_list$geo_loc_name     <- NA_character_

  if (!is.null(taxa_of_interest) && length(taxa_of_interest) > 0) {
    pat <- paste0("^(", paste(taxa_of_interest, collapse = "|"), ")\\b")
    accession_list <- accession_list[grepl(pat, accession_list$organism), ]
  }

# turn "" into NA
name_cols <- c("specimen_voucher", "strain", "isolate")
accession_list[name_cols] <- lapply(accession_list[name_cols], function(x) {
  x[x == ""] <- NA_character_
  x
})

accession_list <- accession_list %>%
  dplyr::mutate(strain.standard = dplyr::coalesce(specimen_voucher, strain, isolate, Accession))

remove_char_pattern <- "[><\\s:;_\\-\\.()&|#/\\\\,'\"!?\\[\\]{}+=%\\*\\^~@$]"
accession_list$strain.standard <- stringr::str_remove_all(
  accession_list$strain.standard,
  remove_char_pattern
)


  accession_list$strain.standard.type <- ifelse(
    !is.na(accession_list$type_material) & accession_list$type_material != "",
    paste0(accession_list$strain.standard, ".TYPE"),
    accession_list$strain.standard
  )

  accession_list$org_name <- gsub("\\s+", "\\.", accession_list$organism)

  write.csv(
    accession_list,
    paste0("./metadata_files/all_accessions_pulled_metadata_", project_name, "_curated_basic.csv"),
    row.names = FALSE
  )
  cat("Wrote basic curated metadata.\n")
}


curate_metadata_regions <- function(project_name,
                                    mapping_file = "./aRborist/example_data/region_replacement_patterns.csv") {
  # 0) read the input from the basic curation step
  infile  <- paste0("./metadata_files/all_accessions_pulled_metadata_",
                    project_name, "_curated_basic.csv")
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
            "'./aRborist/example_data/region_replacement_patterns.csv'.")
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

    #accession title column (left off for now, but same pattern works)
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


filter_metadata <- function(project_name, taxa_of_interest, acc_to_exclude = NULL) {
  accession_list <- read.csv(paste0("./metadata_files/all_accessions_pulled_metadata_", project_name, "_curated.csv"),
                             header = TRUE)
  accession_list$complex_name <- accession_list$organism
  accession_list <- accession_list %>%
    tidyr::separate(complex_name, c("listed.genus", "V1", "V2", "V3"), sep = " ", extra = "merge") %>%
    dplyr::select(-V1, -V2, -V3)

  accession_list <- accession_list[accession_list$listed.genus %in% taxa_of_interest, ]

  if (!is.null(acc_to_exclude) && length(acc_to_exclude) > 0 && any(acc_to_exclude != "")) {
    accession_list <- accession_list[!accession_list$Accession %in% acc_to_exclude, ]
  }

  accession_list[] <- lapply(accession_list, function(x) iconv(x, from = "UTF-8", to = "UTF-8", sub = "byte"))

  write.csv(accession_list,
            paste0("./metadata_files/all_accessions_pulled_metadata_", project_name, "_curated_filtered.csv"),
            row.names = FALSE)
  cat("Wrote curated + filtered metadata.\n")
}

select_regions <- function(project_name, min_region_requirement = length(regions_to_include)) {
  accession_list <- read.csv(paste0("./metadata_files/all_accessions_pulled_metadata_", project_name, "_curated_filtered.csv"),
                             header = TRUE)

  multifasta_prep <- accession_list[accession_list$region.standard %in% regions_to_include, ]
  write.csv(multifasta_prep, paste0("./metadata_files/selected_accessions_metadata_", project_name, ".csv"),
            row.names = FALSE)

  multifasta_prep_complete <- subset(multifasta_prep, select = c(strain.standard.type, organism, Accession, region.standard))

  complete_region_attendance <- multifasta_prep_complete %>%
    tidyr::pivot_wider(names_from = "region.standard", values_from = "Accession", values_fn = list) %>%
    dplyr::mutate(across(all_of(regions_to_include), ~replace(., lengths(.) == 0, NA)))

  multifasta_prep_select <- dplyr::distinct(multifasta_prep_complete, strain.standard.type, region.standard, .keep_all = TRUE)
  select_region_attendance <- multifasta_prep_select %>%
    tidyr::pivot_wider(names_from = "region.standard", values_from = "Accession")

  select_region_attendance_filtered <- select_region_attendance %>%
    dplyr::mutate(total = rowSums(!is.na(dplyr::select(., -strain.standard.type, -organism)))) %>%
    dplyr::filter(total >= min_region_requirement) %>%
    dplyr::select(-total)

  write.csv(select_region_attendance_filtered,
            paste0("./metadata_files/Region_attendance_sheet_", project_name, ".csv"),
            row.names = FALSE)
  cat("Wrote region attendance sheet.\n")
}

create_multifastas <- function(project_name, output_dir = "./multifastas/") {
  select_region_attendance_filtered <- read.csv(paste0("./metadata_files/Region_attendance_sheet_", project_name, ".csv"),
                                                header = TRUE)
  filtered_accessions <- data.frame(Accession = unlist(select_region_attendance_filtered[, -c(1:2)]),
                                    row.names = NULL)

  accession_list <- read.csv(paste0("./metadata_files/all_accessions_pulled_metadata_", project_name, "_curated_filtered.csv"),
                             header = TRUE)
  filtered_accessions_metadata <- merge(filtered_accessions, accession_list, by = "Accession", all.x = TRUE)

  multifasta_prep_simple <- subset(filtered_accessions_metadata,
                                   select = c(strain.standard, organism, Accession, region.standard, fasta.header.type, sequence))
  region_dfs <- split(multifasta_prep_simple, with(multifasta_prep_simple, region.standard), drop = TRUE)

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  for (i in seq_along(region_dfs)) {
    locus_df <- region_dfs[[i]]
    locus_df <- locus_df[order(locus_df$fasta.header.type, decreasing = FALSE), ]
    region.name <- unique(locus_df$region.standard)
    seqs_fasta <- c(rbind(locus_df$fasta.header.type, locus_df$sequence))
    filename <- paste0(output_dir, region.name, ".fasta")
    write.table(seqs_fasta, filename, row.names = FALSE, quote = FALSE, col.names = FALSE)
    message("Created multifasta for region: ", region.name)
  }
  cat("Multifasta files created in ", output_dir, "\n")
}


# Create Projects/<project_name>/ and run your normal structure inside it.
start_project <- function(projects_dir = file.path(path.expand("~"), "aRborist_Projects"),
                          project_name) {
  if (!dir.exists(projects_dir)) dir.create(projects_dir, recursive = TRUE)

  assign("project_name", project_name, envir = .GlobalEnv)
  assign("base_dir", file.path(projects_dir, project_name), envir = .GlobalEnv)
  if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)
  setup_project_structure(base_dir)
  invisible(normalizePath(base_dir))
}

# Save the run options for this project so it's reproducible later
save_project_config <- function(project_dir = getwd(),
                                project_name,
                                taxa_of_interest,
                                regions_to_include = NULL,
                                max_acc_per_taxa = NULL,           # <— no default
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


# Simple wrappers (used by the runner)
ncbi_data_fetch <- function(project_name, max_acc_per_taxa = "max", taxa_of_interest = NULL) {
  get_accessions_for_all_taxa(taxa_of_interest, max_acc_per_taxa)
  retrieve_ncbi_metadata(project_name)
}

data_curate <- function(project_name, taxa_of_interest = NULL, acc_to_exclude = NULL,
                        min_region_requirement = length(regions_to_include)) {
  merge_metadata_with_custom_file(project_name)
  curate_metadata(project_name)
  filter_metadata(project_name, taxa_of_interest, acc_to_exclude)
  select_regions(project_name, min_region_requirement)
  create_multifastas(project_name)
}
