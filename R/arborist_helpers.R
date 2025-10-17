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


# Defaults for entrez search term (user can override in run script)
if (!exists("organism_scope")) {
  organism_scope <- "txid4751[Organism:exp]"  # default is Fungi + descendants
}
if (!exists("search_options")) {
  search_options <- "(biomol_genomic[PROP] AND (100[SLEN]:5000[SLEN])) NOT Contig[All Fields] NOT scaffold[All Fields] NOT genome[All Fields]" # this is my preferred default
}

compose_entrez_term <- function(taxon, organism_scope = NULL, extra_filters = NULL) {
  if (is.null(organism_scope) && exists("organism_scope", .GlobalEnv))
    organism_scope <- get("organism_scope", envir = .GlobalEnv)
  if (is.null(extra_filters) && exists("search_options", .GlobalEnv))
    extra_filters <- get("search_options", envir = .GlobalEnv)

  parts <- c(sprintf('("%s"[All Fields])', taxon), organism_scope, extra_filters)
  paste(Filter(function(z) !is.null(z) && nzchar(z), parts), collapse = " AND ")
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

  # resolving defaults from globals, if needed
  if (is.null(organism_scope) && exists("organism_scope", .GlobalEnv))
    organism_scope <- get("organism_scope", envir = .GlobalEnv)
  if (is.null(extra_filters) && exists("search_options", .GlobalEnv))
    extra_filters <- get("search_options", envir = .GlobalEnv)

  filters <- compose_entrez_term(taxon, organism_scope, extra_filters)

  search <- rentrez::entrez_search(db = "nucleotide", term = filters, retmax = 9999)

  # No results returns a zero-row df
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

    temp_filename <- paste0("./temp_files/temp_file_accessions_from_", taxon, ".txt")
    if (file.exists(temp_filename)) file.remove(temp_filename)

    for (seq_start in seq(0, pull_n - 1, by = 50)) {
      recs <- rentrez::entrez_fetch(db = "nuccore",
                                    web_history = large_search$web_history,
                                    rettype = "acc",
                                    retmax = min(50, pull_n - seq_start),
                                    retstart = seq_start)
      cat(recs, file = temp_filename, append = TRUE)
      Sys.sleep(get_sleep_duration())
      cat(seq_start + min(49, pull_n - seq_start - 1), "accessions recorded\r")
    }

    large_temp_df <- read.table(temp_filename, stringsAsFactors = FALSE)
    colnames(large_temp_df) <- c("Accession")
    large_temp_df$genus <- taxon
    if (file.exists(temp_filename)) file.remove(temp_filename)
    cat("Accession retrieval for", taxon, "successful\n\n")
    return(large_temp_df)
  }

  # Normal path (≤ 9,999)
  ids <- search$ids
  if (is.finite(max_n)) ids <- utils::head(ids, max_n)

  if (length(ids) <= 300) {
    summary <- rentrez::entrez_summary(db = "nuccore", id = ids)
    Sys.sleep(get_sleep_duration())
  } else {
    summary <- list()
    index <- split(seq_along(ids), ceiling(seq_along(ids) / 300))
    for (p in index) {
      summary[p] <- rentrez::entrez_summary(db = "nuccore", id = ids[p])
      Sys.sleep(get_sleep_duration())
    }
    class(summary) <- c("esummary_list", "list")
  }

  tempdf <- data.frame(Accession = unname(rentrez::extract_from_esummary(summary, "caption")),
                       stringsAsFactors = FALSE)
  tempdf$genus <- taxon
  cat("Search complete for", taxon, " (", nrow(tempdf), " accessions)\n", sep = "")
  return(tempdf)
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
  if (!exists("metadata_categories_keep", .GlobalEnv)) {
    metadata_categories_keep <- c(
      "GBSeq_locus","GBSeq_length","GBSeq_strandedness","GBSeq_moltype",
      "GBSeq_update.date","GBSeq_create.date","GBSeq_definition",
      "GBSeq_accession.version","GBSeq_project","GBSeq_organism","GBSeq_taxonomy",
      "GBSeq_sequence","GBSeq_feature.table","_title","_journal","ref_id","pubmed"
    )
  }

  out.xml <- rentrez::entrez_fetch(db = "nuccore", id = accession, rettype = "xml")
  list.out <- XML::xmlToList(out.xml)
  if (length(list.out) == 0L) {
    return(data.frame(Accession = accession, stringsAsFactors = FALSE))
  }
  accession_dfs <- lapply(list.out, data.frame, stringsAsFactors = FALSE)
  all_metadata_df <- accession_dfs[[1]]as.data.frame(metadata_entry, stringsAsFactors = FALSE)
  }

  keep_cols <- intersect(metadata_categories_keep, names(all_metadata_df))
  select_metadata_df <- all_metadata_df[, keep_cols, drop = FALSE]

  basic_pat   <- "GBSeq_locus|GBSeq_length|GBSeq_strandedness|GBSeq_update\\.date|GBSeq_create\\.date|GBSeq_definition|GBSeq_accession\\.version|GBSeq_project|GBSeq_organism|GBSeq_taxonomy|GBSeq_sequence"
  feature_pat <- "GBSeq_feature\\.table|GBSeq_feature-table"

  basic_cols   <- grep(basic_pat,   names(select_metadata_df), value = TRUE)
  feature_cols <- grep(feature_pat, names(select_metadata_df), value = TRUE)

  basic_info_df    <- select_metadata_df[, basic_cols,   drop = FALSE]
  feature_table_df <- select_metadata_df[, feature_cols, drop = FALSE]

  if (ncol(feature_table_df) > 0L) {
    name_idx <- which(grepl("_name$", colnames(feature_table_df)))
    valid_pairs <- name_idx[name_idx + 1L <= ncol(feature_table_df)]
    if (length(valid_pairs)) {
      feat_trans <- data.frame(feature_table_df[, valid_pairs + 1L, drop = FALSE],
                               check.names = FALSE, stringsAsFactors = FALSE)
      colnames(feat_trans) <- feature_table_df[, valid_pairs, drop = TRUE]
    } else {
      feat_trans <- basic_info_df[, 0, drop = FALSE]
    }
  } else {
    feat_trans <- basic_info_df[, 0, drop = FALSE]
  }

  metadata_entry <- cbind(basic_info_df, feat_trans)
  names(metadata_entry) <- sub("^GBSeq_", "", names(metadata_entry))
  if ("locus" %in% names(metadata_entry)) data.table::setnames(metadata_entry, "locus", "Accession") else metadata_entry$Accession <- accession
  if ("definition" %in% names(metadata_entry)) data.table::setnames(metadata_entry, "definition", "accession_title")

  as.data.frame(metadata_entry, stringsAsFactors = FALSE)
}


  # combine basic + transformed features
  metadata_entry <- cbind(basic_info_df, feat_trans)

  # normalize column names and key fields 
  names(metadata_entry) <- sub("^GBSeq_", "", names(metadata_entry))
  if ("locus" %in% names(metadata_entry)) {
    data.table::setnames(metadata_entry, "locus", "Accession")
  } else {
    # ensure Accession column exists even if locus was absent
    metadata_entry$Accession <- accession
  }
  if ("definition" %in% names(metadata_entry)) {
    data.table::setnames(metadata_entry, "definition", "accession_title")
  }

  as.data.frame(metadata_entry, stringsAsFactors = FALSE)
}


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
curate_metadata <- function(project_name) {
  metadata_file_path <- paste0("./metadata_files/all_accessions_pulled_metadata_", project_name, ".csv")
  accession_list <- read.csv(metadata_file_path, header = TRUE)

  # strain.standard naming priority: voucher > strain > isolate > Accession
  accession_list <- accession_list %>%
    dplyr::mutate(strain.standard = dplyr::coalesce(specimen_voucher, strain, isolate, Accession))

  # strip punctuation/whitespace differences
  remove_char <- c("\\>", "\\<", "\\s+", ":", ";", "_", "-", "\\.", "\\(", "\\)", "&", "\\|")
  accession_list$strain.standard <- stringr::str_remove_all(accession_list$strain.standard, paste(remove_char, collapse = "|"))

  # mark TYPE if ANYTHING in type_material category
  accession_list$strain.standard.type <- ifelse(!is.na(accession_list$type_material),
                                                paste(accession_list$strain.standard, ".TYPE", sep = ""),
                                                accession_list$strain.standard)

  # gene/product/title standardization (add/modify as needed - this is def not an exhaustive list)
  accession_list$gene.standard <- accession_list$gene
  accession_list <- accession_list %>%
    dplyr::mutate(gene.standard = dplyr::case_when(
      stringr::str_detect(gene.standard, stringr::regex("tef-*\\d*", ignore_case = TRUE)) ~ "TEF",
      stringr::str_detect(gene.standard, stringr::regex("EF1-alpha", ignore_case = TRUE)) ~ "TEF",
      stringr::str_detect(gene.standard, stringr::regex("ef1a", ignore_case = TRUE)) ~ "TEF",
      stringr::str_detect(gene.standard, stringr::regex("b-tub", ignore_case = TRUE)) ~ "BTUB",
      stringr::str_detect(gene.standard, stringr::regex("TUB2", ignore_case = TRUE)) ~ "BTUB",
      stringr::str_detect(gene.standard, stringr::regex("RBP2", ignore_case = TRUE)) ~ "RPB2",
      stringr::str_detect(gene.standard, stringr::regex("RPB2", ignore_case = TRUE)) ~ "RPB2",
      stringr::str_detect(gene.standard, stringr::regex("RBP1", ignore_case = TRUE)) ~ "RPB1",
      stringr::str_detect(gene.standard, stringr::regex("RPB1", ignore_case = TRUE)) ~ "RPB1",
      stringr::str_detect(gene.standard, stringr::regex("act[16]", ignore_case = TRUE)) ~ "actin",
      TRUE ~ gene.standard
    ))

  accession_list$product.standard <- accession_list$product
  accession_list <- accession_list %>%
    dplyr::mutate(product.standard = dplyr::case_when(
      stringr::str_detect(product.standard, stringr::regex("elongation factor 1", ignore_case = TRUE)) ~ "TEF",
      stringr::str_detect(product.standard, stringr::regex("18S", ignore_case = TRUE)) ~ "SSU",
      stringr::str_detect(product.standard, stringr::regex("16S", ignore_case = TRUE)) ~ "SSU",
      stringr::str_detect(product.standard, stringr::regex("26S", ignore_case = TRUE)) ~ "LSU",
      stringr::str_detect(product.standard, stringr::regex("28S", ignore_case = TRUE)) ~ "LSU",
      stringr::str_detect(product.standard, stringr::regex("actin beta", ignore_case = TRUE)) ~ "actin",
      stringr::str_detect(product.standard, stringr::regex("beta-tubulin", ignore_case = TRUE)) ~ "BTUB",
      stringr::str_detect(product.standard, stringr::regex("licensing\\D*7\\D*", ignore_case = TRUE)) ~ "MCM7",
      stringr::str_detect(product.standard, stringr::regex("polymerase II larg[est]*", ignore_case = TRUE)) ~ "RPB1",
      stringr::str_detect(product.standard, stringr::regex("polymerase II second largest", ignore_case = TRUE)) ~ "RPB2",
      stringr::str_detect(product.standard, stringr::regex("small subunit ribosomal", ignore_case = TRUE)) ~ "SSU",
      stringr::str_detect(product.standard, stringr::regex("^large[st]* subunit ribosomal", ignore_case = TRUE)) ~ "LSU",
      stringr::str_detect(product.standard, stringr::regex("internal transcribed", ignore_case = TRUE)) ~ "ITS",
      TRUE ~ product.standard
    ))

  accession_list$acc_title.standard <- accession_list$accession_title
  accession_list <- accession_list %>%
    dplyr::mutate(acc_title.standard = dplyr::case_when(
      stringr::str_detect(acc_title.standard, stringr::regex("actin beta", ignore_case = TRUE)) ~ "actin",
      stringr::str_detect(acc_title.standard, stringr::regex("beta-*\\s*tubulin", ignore_case = TRUE)) ~ "BTUB",
      stringr::str_detect(acc_title.standard, stringr::regex("licensing\\D*7\\D*", ignore_case = TRUE)) ~ "MCM7",
      stringr::str_detect(acc_title.standard, stringr::regex("RPB1", ignore_case = TRUE)) ~ "RPB1",
      stringr::str_detect(acc_title.standard, stringr::regex("RPB2", ignore_case = TRUE)) ~ "RPB2",
      stringr::str_detect(acc_title.standard, stringr::regex("internal transcribed.*complete sequence", ignore_case = TRUE)) ~ "ITS",
      stringr::str_detect(acc_title.standard, stringr::regex("translation elongation", ignore_case = TRUE)) ~ "TEF",
      TRUE ~ acc_title.standard
    ))

  accession_list <- accession_list %>%
    dplyr::mutate(region.standard = dplyr::coalesce(gene.standard, product.standard, acc_title.standard))

  accession_list$org_name <- gsub("\\s+", "\\.", accession_list$organism)
  accession_list$fasta.header <- paste(">", accession_list$org_name, "_", accession_list$strain.standard, sep = "")
  accession_list$fasta.header.type <- paste(">", accession_list$org_name, "_", accession_list$strain.standard.type, sep = "")

  write.csv(accession_list,
            paste0("./metadata_files/all_accessions_pulled_metadata_", project_name, "_curated.csv"),
            row.names = FALSE)
  cat("Wrote curated metadata.\n")
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


