# ============================================================
# 04_tss_consensus.R
#
# Aim
# ---
# Identify proximal transcription start sites (TSSs) for the
# harmonised GcvB candidate target genes and derive a consensus
# between two TSS-enriched RNA-seq conditions:
#
#   WT plate
#   WT broth
#
#
# WORKFLOW
# --------
#
# 1. Import harmonised candidate loci.
# 2. Add genomic coordinates and strand from the reference.
# 3. Import WT plate and WT broth TSS-enriched GTF files.
# 4. For each gene, identify same-strand upstream TSS calls.
# 5. Search up to 2000 nt upstream.
# 6. Define proximal TSS calls as <= 500 nt upstream.
# 7. Select the nearest proximal TSS independently for:
#       - WT plate
#       - WT broth
# 8. Compare selected coordinates.
# 9. Calls within +/- 10 nt are considered concordant.
# 10. Assign a representative TSS where the result is
#     interpretable.
#
#
# CONSENSUS RULES
# ---------------
#
# Both_concordant
#   Proximal TSS present in both conditions and selected
#   coordinates differ by <= 10 nt.
#
# Both_discordant
#   Proximal TSS present in both conditions but selected
#   coordinates differ by > 10 nt.
#
# WT_plate_only
#   Proximal TSS detected only in WT plate.
#
# WT_broth_only
#   Proximal TSS detected only in WT broth.
#
# No_proximal_TSS
#   No proximal TSS detected in either condition.
#
#
# REPRESENTATIVE TSS
# ------------------
#
# Both_concordant:
#   retain the WT plate coordinate.
#
# WT_plate_only:
#   retain the WT plate coordinate.
#
# WT_broth_only:
#   retain the WT broth coordinate.
#
# Both_discordant:
#   no representative coordinate is assigned.
#
# No_proximal_TSS:
#   no representative coordinate is assigned.
#
#
# INPUTS
# ------
#
# results/tables/target_candidates_annotated.csv
#
# data/reference/reference_gene_coordinates.csv
#
# data/raw/tss/
#   WT plate TSS-enriched GTF
#   WT broth TSS-enriched GTF
#
#
# OUTPUTS
# -------
#
# results/tables/tss_consensus.csv
# results/tables/tss_candidate_calls.csv
# results/tables/tss_support_summary.csv
# results/tables/tss_consensus_qc.csv
#
# ============================================================



# ============================================================
# USER SETTINGS
# ============================================================

# Set the working directory to the root of the analysis project
# before running this script.
#
# Example:
# setwd("/path/to/sRNA-target-pipeline")

project_dir <- getwd()

candidate_filename <-
  "target_candidates_annotated.csv"


reference_filename <-
  "reference_gene_coordinates.csv"



# Terms used to identify the two GTF files.
#
# Every term in each vector must occur in the filename.

wt_plate_terms <- c(
  "2331WT",
  "TSS",
  "enriched"
)


wt_broth_terms <- c(
  "2331BR",
  "TSS",
  "enriched"
)



# Maximum distance upstream searched for any candidate TSS.

upstream_search_max <- 2000L


# Maximum distance upstream for a TSS to be considered proximal.

proximal_tss_max <- 500L


# Maximum coordinate difference between plate and broth
# selected TSS calls for them to be considered concordant.

tss_concordance_tolerance <- 10L



# ============================================================
# PACKAGES
# ============================================================

required_pkgs <- c(
  "readr",
  "dplyr",
  "stringr",
  "tibble"
)


missing_pkgs <- required_pkgs[
  !vapply(
    required_pkgs,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]


if (length(missing_pkgs) > 0) {
  
  install.packages(
    missing_pkgs
  )
}


suppressPackageStartupMessages({
  
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  
})



# ============================================================
# PATHS
# ============================================================

project_dir <- getwd()


candidate_file <- file.path(
  project_dir,
  "results",
  "tables",
  candidate_filename
)


reference_file <- file.path(
  project_dir,
  "data",
  "reference",
  reference_filename
)


tss_dir <- file.path(
  project_dir,
  "data",
  "raw",
  "tss"
)


results_dir <- file.path(
  project_dir,
  "results",
  "tables"
)


dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)



# ============================================================
# CHECK INPUT PATHS
# ============================================================

if (!file.exists(candidate_file)) {
  
  stop(
    "Candidate file not found:\n",
    candidate_file,
    "\n\nRun 01_harmonise_target_predictions.R first."
  )
}


if (!file.exists(reference_file)) {
  
  stop(
    "Reference coordinate file not found:\n",
    reference_file
  )
}


if (!dir.exists(tss_dir)) {
  
  stop(
    "TSS directory not found:\n",
    tss_dir
  )
}



# ============================================================
# HELPER: CHECK REQUIRED COLUMNS
# ============================================================

check_columns <- function(
    df,
    required,
    label
) {
  
  missing_columns <- setdiff(
    required,
    names(df)
  )
  
  
  if (length(missing_columns) > 0) {
    
    stop(
      label,
      " is missing required column(s):\n",
      paste(
        missing_columns,
        collapse = ", "
      ),
      "\n\nAvailable columns:\n",
      paste(
        names(df),
        collapse = ", "
      )
    )
  }
}



# ============================================================
# HELPER: LOCATE ONE GTF FILE
# ============================================================

find_gtf <- function(
    directory,
    required_terms,
    label
) {
  
  files <- list.files(
    directory,
    pattern = "\\.gtf$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  
  if (length(files) == 0) {
    
    stop(
      "No GTF files found in:\n",
      directory
    )
  }
  
  
  hits <- files
  
  
  for (term in required_terms) {
    
    hits <- hits[
      grepl(
        term,
        basename(hits),
        ignore.case = TRUE
      )
    ]
  }
  
  
  if (length(hits) == 0) {
    
    stop(
      "Could not identify ",
      label,
      ".\n\nTerms searched:\n",
      paste(
        required_terms,
        collapse = ", "
      ),
      "\n\nAvailable GTF files:\n",
      paste(
        basename(files),
        collapse = "\n"
      )
    )
  }
  
  
  # If duplicate/copy files exist, prefer the simplest filename.
  
  hits <- hits[
    order(
      nchar(basename(hits)),
      basename(hits)
    )
  ]
  
  
  if (length(hits) > 1) {
    
    warning(
      "Multiple files matched ",
      label,
      ". Using:\n",
      basename(hits[1])
    )
  }
  
  
  message(
    label,
    ": ",
    basename(hits[1])
  )
  
  
  return(hits[1])
}



# ============================================================
# LOCATE TSS FILES
# ============================================================

wt_plate_gtf <- find_gtf(
  directory = tss_dir,
  required_terms = wt_plate_terms,
  label = "WT plate TSS GTF"
)


wt_broth_gtf <- find_gtf(
  directory = tss_dir,
  required_terms = wt_broth_terms,
  label = "WT broth TSS GTF"
)



# ============================================================
# IMPORT CANDIDATE LOCI
# ============================================================

candidates <- readr::read_csv(
  candidate_file,
  show_col_types = FALSE
)


check_columns(
  candidates,
  c(
    "reference_locus",
    "gene",
    "annotation"
  ),
  "Candidate table"
)


candidates <- candidates |>
  dplyr::transmute(
    
    reference_locus =
      stringr::str_trim(
        as.character(reference_locus)
      ),
    
    gene =
      stringr::str_trim(
        as.character(gene)
      ),
    
    annotation =
      as.character(annotation)
    
  ) |>
  dplyr::filter(
    !is.na(reference_locus),
    reference_locus != ""
  ) |>
  dplyr::distinct(
    reference_locus,
    .keep_all = TRUE
  )



# ============================================================
# IMPORT REFERENCE COORDINATES
# ============================================================

reference <- readr::read_csv(
  reference_file,
  show_col_types = FALSE
)


check_columns(
  reference,
  c(
    "reference_locus",
    "strand",
    "gene_start",
    "gene_end"
  ),
  "Reference coordinate table"
)


reference <- reference |>
  dplyr::transmute(
    
    reference_locus =
      stringr::str_trim(
        as.character(reference_locus)
      ),
    
    strand =
      stringr::str_trim(
        as.character(strand)
      ),
    
    gene_start =
      suppressWarnings(
        as.integer(gene_start)
      ),
    
    gene_end =
      suppressWarnings(
        as.integer(gene_end)
      )
    
  ) |>
  dplyr::filter(
    !is.na(reference_locus),
    reference_locus != ""
  ) |>
  dplyr::distinct(
    reference_locus,
    .keep_all = TRUE
  )



# ============================================================
# ADD GENOMIC COORDINATES TO CANDIDATE GENES
# ============================================================

genes_all <- candidates |>
  dplyr::left_join(
    reference,
    by = "reference_locus"
  )



# ============================================================
# IDENTIFY UNMAPPED CANDIDATES
# ============================================================

unmapped_genes <- genes_all |>
  dplyr::filter(
    is.na(strand) |
      is.na(gene_start) |
      is.na(gene_end)
  )


if (nrow(unmapped_genes) > 0) {
  
  warning(
    nrow(unmapped_genes),
    " candidate locus/loci could not be assigned complete ",
    "reference coordinates."
  )
}



# ============================================================
# RETAIN CANDIDATES WITH VALID COORDINATES
# ============================================================

genes <- genes_all |>
  dplyr::filter(
    !is.na(strand),
    !is.na(gene_start),
    !is.na(gene_end)
  )


if (nrow(genes) == 0) {
  
  stop(
    "No candidate loci could be mapped to the reference ",
    "coordinate table."
  )
}


if (any(!genes$strand %in% c("+", "-"))) {
  
  bad_strands <- unique(
    genes$strand[
      !genes$strand %in% c("+", "-")
    ]
  )
  
  stop(
    "Unexpected strand values in reference table:\n",
    paste(
      bad_strands,
      collapse = ", "
    )
  )
}



# ============================================================
# DEFINE TRANSLATION START
# ============================================================
#
# For a + strand gene:
#
#   translation start = lowest genomic coordinate
#
# For a - strand gene:
#
#   translation start = highest genomic coordinate

genes <- genes |>
  dplyr::mutate(
    
    gene_left =
      pmin(
        gene_start,
        gene_end
      ),
    
    gene_right =
      pmax(
        gene_start,
        gene_end
      ),
    
    translation_start =
      dplyr::case_when(
        
        strand == "+" ~
          gene_left,
        
        strand == "-" ~
          gene_right,
        
        TRUE ~
          NA_integer_
      )
  )



# ============================================================
# READ TSS GTF
# ============================================================

read_tss_gtf <- function(
    path,
    condition_label
) {
  
  x <- readr::read_tsv(
    path,
    comment = "#",
    col_names = c(
      "seqname",
      "source",
      "feature",
      "start",
      "end",
      "score",
      "strand",
      "frame",
      "attribute"
    ),
    col_types = readr::cols(
      seqname = readr::col_character(),
      source = readr::col_character(),
      feature = readr::col_character(),
      start = readr::col_integer(),
      end = readr::col_integer(),
      score = readr::col_character(),
      strand = readr::col_character(),
      frame = readr::col_character(),
      attribute = readr::col_character()
    ),
    progress = FALSE,
    show_col_types = FALSE
  )
  
  
  x <- x |>
    dplyr::mutate(
      
      condition =
        condition_label,
      
      # A TSS is expected to be a single genomic position.
      #
      # If a feature spans >1 nt, use its transcript-facing
      # boundary.
      
      TSS_coord =
        dplyr::case_when(
          
          strand == "+" ~
            start,
          
          strand == "-" ~
            end,
          
          TRUE ~
            NA_integer_
        ),
      
      score_numeric =
        suppressWarnings(
          as.numeric(score)
        )
    ) |>
    dplyr::filter(
      strand %in% c("+", "-"),
      !is.na(TSS_coord)
    )
  
  
  return(x)
}



# ============================================================
# IMPORT TSS DATA
# ============================================================

wt_plate <- read_tss_gtf(
  path = wt_plate_gtf,
  condition_label = "WT_plate"
)


wt_broth <- read_tss_gtf(
  path = wt_broth_gtf,
  condition_label = "WT_broth"
)



cat(
  "\nRaw WT plate TSS calls: ",
  nrow(wt_plate),
  "\n",
  sep = ""
)


cat(
  "Raw WT broth TSS calls: ",
  nrow(wt_broth),
  "\n",
  sep = ""
)



# ============================================================
# HELPER: FIND UPSTREAM TSS CANDIDATES
# ============================================================
#
# Positive upstream_distance means the TSS is upstream of the
# translation start in transcriptional orientation.
#
#
# + strand:
#
#      TSS -----> gene
#
#      distance =
#          translation_start - TSS_coord
#
#
# - strand:
#
#      gene <----- TSS
#
#      distance =
#          TSS_coord - translation_start
#
#
# The many-to-many join is intentional:
#
# each gene is initially compared with all same-strand TSS
# calls, after which genomic distance is used to retain only
# plausible upstream calls.

add_tss_distance <- function(
    gene_table,
    tss_table
) {
  
  gene_min <- gene_table |>
    dplyr::select(
      reference_locus,
      gene,
      annotation,
      strand,
      gene_start,
      gene_end,
      translation_start
    )
  
  
  tss_min <- tss_table |>
    dplyr::transmute(
      
      condition = condition,
      
      TSS_coord = TSS_coord,
      
      TSS_strand = strand,
      
      score_numeric = score_numeric,
      
      tss_source = source,
      
      tss_feature = feature,
      
      tss_attribute = attribute
    )
  
  
  joined <- dplyr::inner_join(
    gene_min,
    tss_min,
    by = c(
      "strand" = "TSS_strand"
    ),
    relationship = "many-to-many"
  )
  
  
  joined <- joined |>
    dplyr::mutate(
      
      upstream_distance =
        dplyr::case_when(
          
          strand == "+" ~
            translation_start - TSS_coord,
          
          strand == "-" ~
            TSS_coord - translation_start,
          
          TRUE ~
            NA_integer_
        )
    ) |>
    dplyr::filter(
      !is.na(upstream_distance),
      upstream_distance >= 0L,
      upstream_distance <= upstream_search_max
    )
  
  
  return(joined)
}



# ============================================================
# GENERATE TSS CANDIDATE TABLES
# ============================================================

plate_candidates <- add_tss_distance(
  gene_table = genes,
  tss_table = wt_plate
)


broth_candidates <- add_tss_distance(
  gene_table = genes,
  tss_table = wt_broth
)



# Combined table is retained for inspection/export only.

all_tss_candidates <- dplyr::bind_rows(
  plate_candidates,
  broth_candidates
) |>
  dplyr::arrange(
    reference_locus,
    condition,
    upstream_distance,
    dplyr::desc(score_numeric)
  )



# ============================================================
# HELPER: SELECT NEAREST PROXIMAL TSS
# ============================================================
#
# This avoids slice(), which previously caused method dispatch
# problems in the R environment.
#
# Ties are resolved by:
#
#   1. smallest upstream distance;
#   2. highest numeric TSS score;
#   3. smallest genomic coordinate.
#
# The third criterion simply guarantees deterministic output
# in the unlikely event of a complete tie.

select_nearest_proximal <- function(
    candidate_table
) {
  
  proximal <- candidate_table |>
    dplyr::filter(
      upstream_distance <= proximal_tss_max
    )
  
  
  if (nrow(proximal) == 0) {
    
    return(
      tibble::tibble(
        reference_locus = character(),
        selected_TSS = integer(),
        selected_distance = integer(),
        selected_score = double()
      )
    )
  }
  
  
  proximal <- proximal |>
    dplyr::group_by(
      reference_locus
    ) |>
    dplyr::arrange(
      upstream_distance,
      dplyr::desc(score_numeric),
      TSS_coord,
      .by_group = TRUE
    ) |>
    dplyr::mutate(
      selection_rank =
        dplyr::row_number()
    ) |>
    dplyr::filter(
      selection_rank == 1L
    ) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      
      reference_locus,
      
      selected_TSS =
        as.integer(TSS_coord),
      
      selected_distance =
        as.integer(upstream_distance),
      
      selected_score =
        as.numeric(score_numeric)
    )
  
  
  return(proximal)
}



# ============================================================
# SELECT WT PLATE TSS
# ============================================================

plate_selected <- select_nearest_proximal(
  plate_candidates
)


names(plate_selected)[
  names(plate_selected) == "selected_TSS"
] <- "WT_plate_TSS"


names(plate_selected)[
  names(plate_selected) == "selected_distance"
] <- "WT_plate_distance"


names(plate_selected)[
  names(plate_selected) == "selected_score"
] <- "WT_plate_score"



# ============================================================
# SELECT WT BROTH TSS
# ============================================================

broth_selected <- select_nearest_proximal(
  broth_candidates
)


names(broth_selected)[
  names(broth_selected) == "selected_TSS"
] <- "WT_broth_TSS"


names(broth_selected)[
  names(broth_selected) == "selected_distance"
] <- "WT_broth_distance"


names(broth_selected)[
  names(broth_selected) == "selected_score"
] <- "WT_broth_score"



# ============================================================
# COUNT PROXIMAL TSS CALLS — WT PLATE
# ============================================================
#
# Plate and broth counts are deliberately generated separately.
#
# This removes the previous pivot_wider() dependency and avoids
# objects/columns named WT_plate or WT_broth being created and
# then renamed later.

plate_counts <- plate_candidates |>
  dplyr::filter(
    upstream_distance <= proximal_tss_max
  ) |>
  dplyr::group_by(
    reference_locus
  ) |>
  dplyr::summarise(
    n_proximal_WT_plate =
      dplyr::n(),
    .groups = "drop"
  )



# ============================================================
# COUNT PROXIMAL TSS CALLS — WT BROTH
# ============================================================

broth_counts <- broth_candidates |>
  dplyr::filter(
    upstream_distance <= proximal_tss_max
  ) |>
  dplyr::group_by(
    reference_locus
  ) |>
  dplyr::summarise(
    n_proximal_WT_broth =
      dplyr::n(),
    .groups = "drop"
  )



# ============================================================
# BUILD CONSENSUS TABLE
# ============================================================

tss_consensus <- genes |>
  dplyr::left_join(
    plate_selected,
    by = "reference_locus"
  ) |>
  dplyr::left_join(
    broth_selected,
    by = "reference_locus"
  ) |>
  dplyr::left_join(
    plate_counts,
    by = "reference_locus"
  ) |>
  dplyr::left_join(
    broth_counts,
    by = "reference_locus"
  )



# ============================================================
# REPLACE MISSING COUNTS WITH ZERO
# ============================================================

tss_consensus <- tss_consensus |>
  dplyr::mutate(
    
    n_proximal_WT_plate =
      dplyr::coalesce(
        n_proximal_WT_plate,
        0L
      ),
    
    n_proximal_WT_broth =
      dplyr::coalesce(
        n_proximal_WT_broth,
        0L
      )
  )



# ============================================================
# COMPARE PLATE AND BROTH TSS CALLS
# ============================================================

tss_consensus <- tss_consensus |>
  dplyr::mutate(
    
    TSS_difference =
      dplyr::case_when(
        
        !is.na(WT_plate_TSS) &
          !is.na(WT_broth_TSS) ~
          
          as.integer(
            abs(
              WT_plate_TSS -
                WT_broth_TSS
            )
          ),
        
        TRUE ~
          NA_integer_
      )
  )



# ============================================================
# ASSIGN SUPPORT CLASS
# ============================================================

tss_consensus <- tss_consensus |>
  dplyr::mutate(
    
    TSS_support =
      dplyr::case_when(
        
        !is.na(WT_plate_TSS) &
          !is.na(WT_broth_TSS) &
          TSS_difference <=
          tss_concordance_tolerance ~
          
          "Both_concordant",
        
        
        !is.na(WT_plate_TSS) &
          !is.na(WT_broth_TSS) &
          TSS_difference >
          tss_concordance_tolerance ~
          
          "Both_discordant",
        
        
        !is.na(WT_plate_TSS) &
          is.na(WT_broth_TSS) ~
          
          "WT_plate_only",
        
        
        is.na(WT_plate_TSS) &
          !is.na(WT_broth_TSS) ~
          
          "WT_broth_only",
        
        
        TRUE ~
          
          "No_proximal_TSS"
      )
  )



# ============================================================
# ASSIGN REPRESENTATIVE TSS
# ============================================================
#
# Coordinates are NOT averaged.
#
# For concordant calls, WT plate is retained exactly.

tss_consensus <- tss_consensus |>
  dplyr::mutate(
    
    representative_TSS =
      dplyr::case_when(
        
        TSS_support == "Both_concordant" ~
          WT_plate_TSS,
        
        TSS_support == "WT_plate_only" ~
          WT_plate_TSS,
        
        TSS_support == "WT_broth_only" ~
          WT_broth_TSS,
        
        TRUE ~
          NA_integer_
      )
  )



# ============================================================
# CALCULATE REPRESENTATIVE 5' UTR LENGTH
# ============================================================
#
# This is the distance between the representative TSS and the
# annotated translation start.
#
# It is therefore always expressed as a positive nucleotide
# distance.

tss_consensus <- tss_consensus |>
  dplyr::mutate(
    
    UTR5_length =
      dplyr::case_when(
        
        strand == "+" &
          !is.na(representative_TSS) ~
          
          translation_start -
          representative_TSS,
        
        
        strand == "-" &
          !is.na(representative_TSS) ~
          
          representative_TSS -
          translation_start,
        
        
        TRUE ~
          NA_integer_
      )
  )



# ============================================================
# TRANSCRIPT-RELATIVE TSS POSITION
# ============================================================
#
# Translation start = 0
#
# Upstream nucleotide positions are negative.
#
# Therefore:
#
#   a TSS 100 nt upstream -> -100

tss_consensus <- tss_consensus |>
  dplyr::mutate(
    
    representative_TSS_rel_start =
      dplyr::case_when(
        
        !is.na(UTR5_length) ~
          
          -as.integer(
            UTR5_length
          ),
        
        TRUE ~
          NA_integer_
      )
  )



# ============================================================
# MULTIPLE PROXIMAL TSS INFORMATION
# ============================================================
#
# This is descriptive QC information rather than a manual
# exclusion flag.

tss_consensus <- tss_consensus |>
  dplyr::mutate(
    
    multiple_proximal_TSS =
      (
        n_proximal_WT_plate > 1L |
          n_proximal_WT_broth > 1L
      )
  )



# ============================================================
# ORDER FINAL CONSENSUS TABLE
# ============================================================

tss_consensus <- tss_consensus |>
  dplyr::select(
    
    reference_locus,
    gene,
    annotation,
    
    strand,
    gene_start,
    gene_end,
    translation_start,
    
    WT_plate_TSS,
    WT_plate_distance,
    WT_plate_score,
    n_proximal_WT_plate,
    
    WT_broth_TSS,
    WT_broth_distance,
    WT_broth_score,
    n_proximal_WT_broth,
    
    TSS_difference,
    TSS_support,
    
    representative_TSS,
    UTR5_length,
    representative_TSS_rel_start,
    
    multiple_proximal_TSS
  ) |>
  dplyr::arrange(
    reference_locus
  )



# ============================================================
# SUPPORT SUMMARY
# ============================================================

tss_summary <- tss_consensus |>
  dplyr::group_by(
    TSS_support
  ) |>
  dplyr::summarise(
    n_genes =
      dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    dplyr::desc(n_genes)
  )



# ============================================================
# QC SUMMARY
# ============================================================

qc_summary <- tibble::tibble(
  
  metric = c(
    
    "candidate_loci_input",
    
    "candidate_loci_with_reference_coordinates",
    
    "candidate_loci_without_reference_coordinates",
    
    "raw_WT_plate_TSS_calls",
    
    "raw_WT_broth_TSS_calls",
    
    "candidate_gene_TSS_pairs_within_2000nt",
    
    "proximal_plate_gene_TSS_pairs_within_500nt",
    
    "proximal_broth_gene_TSS_pairs_within_500nt",
    
    "loci_with_multiple_proximal_TSS",
    
    "loci_with_representative_TSS"
  ),
  
  
  value = c(
    
    nrow(candidates),
    
    nrow(genes),
    
    nrow(unmapped_genes),
    
    nrow(wt_plate),
    
    nrow(wt_broth),
    
    nrow(all_tss_candidates),
    
    sum(
      plate_candidates$upstream_distance <=
        proximal_tss_max,
      na.rm = TRUE
    ),
    
    sum(
      broth_candidates$upstream_distance <=
        proximal_tss_max,
      na.rm = TRUE
    ),
    
    sum(
      tss_consensus$multiple_proximal_TSS,
      na.rm = TRUE
    ),
    
    sum(
      !is.na(
        tss_consensus$representative_TSS
      )
    )
  )
)



# ============================================================
# OUTPUT FILES
# ============================================================

consensus_file <- file.path(
  results_dir,
  "tss_consensus.csv"
)


candidate_calls_file <- file.path(
  results_dir,
  "tss_candidate_calls.csv"
)


summary_file <- file.path(
  results_dir,
  "tss_support_summary.csv"
)


qc_file <- file.path(
  results_dir,
  "tss_consensus_qc.csv"
)



# ============================================================
# EXPORT RESULTS
# ============================================================

readr::write_csv(
  tss_consensus,
  consensus_file,
  na = ""
)


readr::write_csv(
  all_tss_candidates,
  candidate_calls_file,
  na = ""
)


readr::write_csv(
  tss_summary,
  summary_file,
  na = ""
)


readr::write_csv(
  qc_summary,
  qc_file,
  na = ""
)



# ============================================================
# CONSOLE SUMMARY
# ============================================================

cat(
  "\n",
  "============================================\n",
  "TSS consensus analysis complete\n",
  "============================================\n\n",
  sep = ""
)


cat(
  "Candidate file:\n",
  basename(candidate_file),
  "\n\n",
  sep = ""
)


cat(
  "WT plate GTF:\n",
  basename(wt_plate_gtf),
  "\n\n",
  sep = ""
)


cat(
  "WT broth GTF:\n",
  basename(wt_broth_gtf),
  "\n\n",
  sep = ""
)


cat(
  "Candidate loci imported: ",
  nrow(candidates),
  "\n",
  sep = ""
)


cat(
  "Candidate loci with reference coordinates: ",
  nrow(genes),
  "\n",
  sep = ""
)


cat(
  "Candidate loci without reference coordinates: ",
  nrow(unmapped_genes),
  "\n\n",
  sep = ""
)


cat(
  "Raw WT plate TSS calls: ",
  nrow(wt_plate),
  "\n",
  sep = ""
)


cat(
  "Raw WT broth TSS calls: ",
  nrow(wt_broth),
  "\n\n",
  sep = ""
)


cat(
  "TSS support classes:\n"
)


print(
  tss_summary
)


cat(
  "\nLoci with multiple proximal TSS calls: ",
  sum(
    tss_consensus$multiple_proximal_TSS,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


cat(
  "Loci with representative TSS: ",
  sum(
    !is.na(
      tss_consensus$representative_TSS
    )
  ),
  "\n\n",
  sep = ""
)


cat(
  "Outputs:\n",
  consensus_file,
  "\n",
  candidate_calls_file,
  "\n",
  summary_file,
  "\n",
  qc_file,
  "\n\n",
  sep = ""
)

