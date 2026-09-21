# ============================================================
# 01_harmonise_target_predictions.R
#
# Aim
# ---
# Import CSV outputs from TargetRNA2, IntaRNA and CopraRNA,
# standardise them into a common format, and generate:
#
#   1) a complete harmonised candidate set containing all loci;
#   2) an annotated candidate set containing only loci with an
#      assigned gene name;
#   3) method-level prediction tables;
#   4) overlap/UpSet-ready tables for downstream analysis.
#
# The complete candidate set preserves all harmonised predictions.
#
# The annotated subset contains loci for which a gene name could
# be assigned from the prediction output, reference annotation,
# or another prediction method.
#
# No candidates are removed from the complete set on the basis
# of annotation quality.
#
#
# REPOSITORY LAYOUT
# -----------------
#
# data/
#   raw/
#     target_predictions/
#       TargetRNA2/
#         *.csv
#       IntaRNA/
#         *.csv
#       CopraRNA/
#         *.csv
#
#   reference/
#       reference_annotation.csv
#
# results/
#   tables/
#
#
# OUTPUTS
# -------
#
# Complete harmonised set:
#
#   target_predictions_method_level_all.csv
#   target_candidates_all.csv
#   target_prediction_overlap_all.csv
#
# Annotated subset:
#
#   target_predictions_method_level_annotated.csv
#   target_candidates_annotated.csv
#   target_prediction_overlap_annotated.csv
#
# QC:
#
#   target_predictions_qc.csv
#   target_prediction_coordinate_qc.csv
#
#
# RUNNING THE SCRIPT
# ------------------
#
# 1. Save/export each prediction-tool result as CSV.
# 2. Place each file in the appropriate input folder.
# 3. Run 00_extract_gbk_annotation.R.
# 4. Check the USER SETTINGS below.
# 5. Run this script from the repository root.
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


# ------------------------------------------------------------
# Reference annotation
# ------------------------------------------------------------

reference_annotation_file <- "reference_annotation.csv"

# Please check the fields in the output file for each prediction tool and
# update the following as required. 

# ------------------------------------------------------------
# TargetRNA2
# ------------------------------------------------------------
#
# Expected columns:
#
# Rank
# Locus Tag 1
# Gene
# Locus Tag APP5 NC_009053.1
# Synonym
# In silico
# Pvalue
# Energy
# mRNA_start
# mRNA_stop
# sRNA_start
# sRNA_stop

targetrna_rank_column <- "Rank"

targetrna_locus_column <- "Locus Tag 1"

targetrna_gene_column <- "Gene"

targetrna_source_locus_column <- "Locus Tag APP5 NC_009053.1"

targetrna_synonym_column <- "Synonym"

targetrna_source_column <- "In silico"

targetrna_p_column <- "Pvalue"

targetrna_energy_column <- "Energy"

targetrna_mrna_start_column <- "mRNA_start"

targetrna_mrna_stop_column <- "mRNA_stop"

targetrna_srna_start_column <- "sRNA_start"

targetrna_srna_stop_column <- "sRNA_stop"



# ------------------------------------------------------------
# IntaRNA
# ------------------------------------------------------------
#
# Expected columns include:
#
# Item
# IntaRNA Rank
# Database
# Species
# Locus 1
# Gene Name Consensus
# Gene Function
# In silico Source
# Locus 2
# Gene
# P-value
# Energy
# Postion mRNA
# Position sRNA
# fdr

intarna_rank_column <- "IntaRNA Rank"

intarna_database_column <- "Database"

intarna_species_column <- "Species"

intarna_locus_column <- "Locus 1"

intarna_gene_column <- "Gene Name Consensus"

intarna_annotation_column <- "Gene Function"

intarna_source_column <- "In silico Source"

intarna_source_locus_column <- "Locus 2"

intarna_source_gene_column <- "Gene"

intarna_p_column <- "P-value"

intarna_energy_column <- "Energy"

intarna_mrna_position_column <- "Postion mRNA"

intarna_srna_position_column <- "Position sRNA"

intarna_fdr_column <- "fdr"



# ------------------------------------------------------------
# CopraRNA
# ------------------------------------------------------------
#
# CopraRNA stores the organism-specific prediction inside the
# corresponding genome/accession column.
#
# For A. pleuropneumoniae MIDG2331:
#
#   NZ_LN908249

coprarna_genome_column <- "NZ_LN908249"

coprarna_fdr_column <- "fdr"

coprarna_combined_p_column <- "p-value"

coprarna_annotation_column <- "Annotation"


# ------------------------------------------------------------
# CopraRNA mRNA coordinate system
# ------------------------------------------------------------
#
# CopraRNA reports organism-specific interaction coordinates in
# the submitted target window rather than directly relative to
# the CDS start.
#
# For the APP analysis, each target window contains 200 nt
# upstream of the CDS start. Because CopraRNA uses 1-based
# positions within that window, window position 201 corresponds
# to CDS-relative position 0.
#
# Therefore:
#
#   CDS-relative coordinate = raw CopraRNA coordinate - 201
#
# Change this setting if a different upstream window is used.

coprarna_upstream_window_nt <- 200L



# ------------------------------------------------------------
# Optional statistical filtering
# ------------------------------------------------------------
#
# Leave as NA_real_ when the supplied CSV files already contain
# the prediction sets to be analysed.
#
# Set a numeric threshold only when additional filtering is
# required.

targetrna_p_max <- NA_real_

intarna_p_max <- NA_real_

coprarna_p_max <- NA_real_



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
  install.packages(missing_pkgs)
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


raw_dir <- file.path(
  project_dir,
  "data",
  "raw",
  "target_predictions"
)


targetrna_dir <- file.path(
  raw_dir,
  "TargetRNA2"
)


intarna_dir <- file.path(
  raw_dir,
  "IntaRNA"
)


coprarna_dir <- file.path(
  raw_dir,
  "CopraRNA"
)


reference_dir <- file.path(
  project_dir,
  "data",
  "reference"
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
# INPUT FILE DISCOVERY
# ============================================================

find_single_csv <- function(directory, label) {
  
  if (!dir.exists(directory)) {
    stop(
      "Input directory does not exist for ",
      label,
      ":\n",
      directory
    )
  }
  
  
  files <- list.files(
    directory,
    pattern = "\\.csv$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  
  if (length(files) == 0) {
    stop(
      "No CSV file was found for ",
      label,
      " in:\n",
      directory
    )
  }
  
  
  if (length(files) > 1) {
    stop(
      "More than one CSV file was found for ",
      label,
      " in:\n",
      directory,
      "\n\nKeep one input CSV in this folder."
    )
  }
  
  
  files
}



targetrna_file <- find_single_csv(
  targetrna_dir,
  "TargetRNA2"
)


intarna_file <- find_single_csv(
  intarna_dir,
  "IntaRNA"
)


coprarna_file <- find_single_csv(
  coprarna_dir,
  "CopraRNA"
)



# ============================================================
# HELPER FUNCTIONS
# ============================================================

normalise_locus <- function(x) {
  
  x <- str_trim(
    as.character(x)
  )
  
  
  x[
    is.na(x) |
      x == "" |
      str_to_upper(x) == "NA" |
      str_to_upper(x) == "N/A"
  ] <- NA_character_
  
  
  str_to_upper(x)
}



clean_gene <- function(x) {
  
  x <- str_trim(
    as.character(x)
  )
  
  
  x[
    is.na(x) |
      x == "" |
      x == "-" |
      str_to_upper(x) == "NA" |
      str_to_upper(x) == "N/A"
  ] <- NA_character_
  
  
  x
}



clean_annotation <- function(x) {
  
  x <- str_trim(
    as.character(x)
  )
  
  
  x[
    is.na(x) |
      x == "" |
      x == "-" |
      str_to_upper(x) == "NA" |
      str_to_upper(x) == "N/A"
  ] <- NA_character_
  
  
  x
}



clean_text <- function(x) {
  
  x <- str_trim(
    as.character(x)
  )
  
  
  x[
    is.na(x) |
      x == "" |
      str_to_upper(x) == "NA" |
      str_to_upper(x) == "N/A"
  ] <- NA_character_
  
  
  x
}



has_usable_gene_name <- function(x) {
  
  x <- clean_gene(x)
  
  
  !is.na(x) &
    !str_detect(
      x,
      regex(
        "^MIDG2331_RS\\d+$",
        ignore_case = TRUE
      )
    ) &
    !str_detect(
      x,
      regex(
        "^APL_\\d+$",
        ignore_case = TRUE
      )
    )
}



check_columns <- function(df, required, label) {
  
  missing_columns <- setdiff(
    required,
    names(df)
  )
  
  
  if (length(missing_columns) > 0) {
    
    stop(
      label,
      " input is missing required column(s):\n",
      paste(
        missing_columns,
        collapse = ", "
      ),
      "\n\nAvailable columns are:\n",
      paste(
        names(df),
        collapse = ", "
      )
    )
  }
}



apply_max_filter <- function(df, column, threshold) {
  
  if (is.na(threshold)) {
    return(df)
  }
  
  
  df %>%
    filter(
      !is.na(.data[[column]]),
      .data[[column]] <= threshold
    )
}



first_non_missing <- function(x) {
  
  x <- x[
    !is.na(x) &
      x != ""
  ]
  
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  
  x[1]
}



most_common_non_missing <- function(x) {
  
  x <- x[
    !is.na(x) &
      x != ""
  ]
  
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  
  counts <- sort(
    table(x),
    decreasing = TRUE
  )
  
  
  names(counts)[1]
}



# ============================================================
# REFERENCE ANNOTATION
# ============================================================

annotation_path <- file.path(
  reference_dir,
  reference_annotation_file
)


if (!file.exists(annotation_path)) {
  
  stop(
    "Reference annotation file not found:\n",
    annotation_path,
    "\n\nRun 00_extract_gbk_annotation.R first."
  )
}



reference_annotation <- read_csv(
  annotation_path,
  show_col_types = FALSE
)


check_columns(
  reference_annotation,
  c(
    "reference_locus",
    "gene",
    "annotation"
  ),
  "Reference annotation"
)



reference_annotation <- reference_annotation %>%
  transmute(
    
    reference_locus =
      normalise_locus(reference_locus),
    
    reference_gene =
      clean_gene(gene),
    
    reference_annotation =
      clean_annotation(annotation)
    
  ) %>%
  filter(
    !is.na(reference_locus)
  ) %>%
  distinct(
    reference_locus,
    .keep_all = TRUE
  )



# ============================================================
# PARSER: TargetRNA2
# ============================================================

parse_targetrna2 <- function(file) {
  
  raw <- read_csv(
    file,
    show_col_types = FALSE,
    na = c("", "NA", "N/A")
  )
  
  
  required <- c(
    targetrna_rank_column,
    targetrna_locus_column,
    targetrna_gene_column,
    targetrna_source_locus_column,
    targetrna_synonym_column,
    targetrna_source_column,
    targetrna_p_column,
    targetrna_energy_column,
    targetrna_mrna_start_column,
    targetrna_mrna_stop_column,
    targetrna_srna_start_column,
    targetrna_srna_stop_column
  )
  
  
  check_columns(
    raw,
    required,
    "TargetRNA2"
  )
  
  
  result <- raw %>%
    transmute(
      
      method =
        "TargetRNA2",
      
      rank =
        suppressWarnings(
          as.integer(
            .data[[targetrna_rank_column]]
          )
        ),
      
      database =
        "TargetRNA2",
      
      species =
        NA_character_,
      
      reference_locus =
        normalise_locus(
          .data[[targetrna_locus_column]]
        ),
      
      source_locus =
        clean_text(
          .data[[targetrna_source_locus_column]]
        ),
      
      synonym =
        clean_text(
          .data[[targetrna_synonym_column]]
        ),
      
      prediction_gene =
        clean_gene(
          .data[[targetrna_gene_column]]
        ),
      
      annotation =
        NA_character_,
      
      in_silico_source =
        clean_text(
          .data[[targetrna_source_column]]
        ),
      
      probability =
        NA_real_,
      
      p_value =
        suppressWarnings(
          as.numeric(
            .data[[targetrna_p_column]]
          )
        ),
      
      fdr =
        NA_real_,
      
      combined_p_value =
        NA_real_,
      
      energy =
        suppressWarnings(
          as.numeric(
            .data[[targetrna_energy_column]]
          )
        ),
      
      mrna_start =
        suppressWarnings(
          as.integer(
            .data[[targetrna_mrna_start_column]]
          )
        ),
      
      mrna_stop =
        suppressWarnings(
          as.integer(
            .data[[targetrna_mrna_stop_column]]
          )
        ),
      
      srna_start =
        suppressWarnings(
          as.integer(
            .data[[targetrna_srna_start_column]]
          )
        ),
      
      srna_stop =
        suppressWarnings(
          as.integer(
            .data[[targetrna_srna_stop_column]]
          )
        ),
      
      source_file =
        basename(file)
    )
  
  
  apply_max_filter(
    result,
    "p_value",
    targetrna_p_max
  )
}



# ============================================================
# HELPER: PARSE IntaRNA POSITION FIELD
# ============================================================

parse_intarna_position <- function(x) {
  
  matched <- str_match(
    as.character(x),
    "^\\s*(-?\\d+)\\s*--\\s*(-?\\d+)\\s*$"
  )
  
  
  tibble(
    
    start =
      suppressWarnings(
        as.integer(
          matched[, 2]
        )
      ),
    
    stop =
      suppressWarnings(
        as.integer(
          matched[, 3]
        )
      )
  )
}



# ============================================================
# PARSER: IntaRNA
# ============================================================

parse_intarna <- function(file) {
  
  raw <- read_csv(
    file,
    show_col_types = FALSE,
    na = c("", "NA", "N/A")
  )
  
  
  required <- c(
    intarna_rank_column,
    intarna_database_column,
    intarna_species_column,
    intarna_locus_column,
    intarna_gene_column,
    intarna_annotation_column,
    intarna_source_column,
    intarna_p_column,
    intarna_energy_column,
    intarna_mrna_position_column,
    intarna_srna_position_column
  )
  
  
  check_columns(
    raw,
    required,
    "IntaRNA"
  )
  
  
  mrna_positions <- parse_intarna_position(
    raw[[intarna_mrna_position_column]]
  )
  
  
  srna_positions <- parse_intarna_position(
    raw[[intarna_srna_position_column]]
  )
  
  
  if (intarna_fdr_column %in% names(raw)) {
    
    parsed_fdr <- suppressWarnings(
      as.numeric(
        raw[[intarna_fdr_column]]
      )
    )
    
  } else {
    
    parsed_fdr <- rep(
      NA_real_,
      nrow(raw)
    )
  }
  
  
  if (intarna_source_locus_column %in% names(raw)) {
    
    source_locus <- clean_text(
      raw[[intarna_source_locus_column]]
    )
    
  } else {
    
    source_locus <- rep(
      NA_character_,
      nrow(raw)
    )
  }
  
  
  if (intarna_source_gene_column %in% names(raw)) {
    
    synonym <- clean_gene(
      raw[[intarna_source_gene_column]]
    )
    
  } else {
    
    synonym <- rep(
      NA_character_,
      nrow(raw)
    )
  }
  
  
  result <- raw %>%
    mutate(
      
      parsed_mrna_start =
        mrna_positions$start,
      
      parsed_mrna_stop =
        mrna_positions$stop,
      
      parsed_srna_start =
        srna_positions$start,
      
      parsed_srna_stop =
        srna_positions$stop,
      
      parsed_fdr =
        parsed_fdr,
      
      parsed_source_locus =
        source_locus,
      
      parsed_synonym =
        synonym
      
    ) %>%
    transmute(
      
      method =
        "IntaRNA",
      
      rank =
        suppressWarnings(
          as.integer(
            .data[[intarna_rank_column]]
          )
        ),
      
      database =
        clean_text(
          .data[[intarna_database_column]]
        ),
      
      species =
        clean_text(
          .data[[intarna_species_column]]
        ),
      
      reference_locus =
        normalise_locus(
          .data[[intarna_locus_column]]
        ),
      
      source_locus =
        parsed_source_locus,
      
      synonym =
        parsed_synonym,
      
      prediction_gene =
        clean_gene(
          .data[[intarna_gene_column]]
        ),
      
      annotation =
        clean_annotation(
          .data[[intarna_annotation_column]]
        ),
      
      in_silico_source =
        clean_text(
          .data[[intarna_source_column]]
        ),
      
      probability =
        NA_real_,
      
      p_value =
        suppressWarnings(
          as.numeric(
            .data[[intarna_p_column]]
          )
        ),
      
      fdr =
        parsed_fdr,
      
      combined_p_value =
        NA_real_,
      
      energy =
        suppressWarnings(
          as.numeric(
            .data[[intarna_energy_column]]
          )
        ),
      
      mrna_start =
        parsed_mrna_start,
      
      mrna_stop =
        parsed_mrna_stop,
      
      srna_start =
        parsed_srna_start,
      
      srna_stop =
        parsed_srna_stop,
      
      source_file =
        basename(file)
    )
  
  
  apply_max_filter(
    result,
    "p_value",
    intarna_p_max
  )
}



# ============================================================
# HELPER: PARSE CopraRNA REFERENCE FIELD
# ============================================================

parse_coprna_reference_field <- function(x) {
  
  # Expected structure:
  #
  # locus(
  #   gene |
  #   energy |
  #   p-value |
  #   mRNA start |
  #   mRNA stop |
  #   sRNA start |
  #   sRNA stop |
  #   additional identifier
  # )
  
  pattern <- paste0(
    "^\\s*([^\\(]+)\\(",
    "([^|]*)\\|",
    "([^|]*)\\|",
    "([^|]*)\\|",
    "([^|]*)\\|",
    "([^|]*)\\|",
    "([^|]*)\\|",
    "([^|]*)\\|",
    "([^\\)]*)\\)\\s*$"
  )
  
  
  matched <- str_match(
    as.character(x),
    pattern
  )
  
  
  tibble(
    
    reference_locus =
      normalise_locus(
        matched[, 2]
      ),
    
    prediction_gene =
      clean_gene(
        matched[, 3]
      ),
    
    energy =
      suppressWarnings(
        as.numeric(
          matched[, 4]
        )
      ),
    
    p_value =
      suppressWarnings(
        as.numeric(
          matched[, 5]
        )
      ),
    
    mrna_start =
      suppressWarnings(
        as.integer(
          matched[, 6]
        )
      ),
    
    mrna_stop =
      suppressWarnings(
        as.integer(
          matched[, 7]
        )
      ),
    
    srna_start =
      suppressWarnings(
        as.integer(
          matched[, 8]
        )
      ),
    
    srna_stop =
      suppressWarnings(
        as.integer(
          matched[, 9]
        )
      ),
    
    source_locus =
      clean_text(
        matched[, 10]
      )
  )
}



# ============================================================
# PARSER: CopraRNA
# ============================================================

parse_coprna <- function(file) {
  
  raw <- read_csv(
    file,
    show_col_types = FALSE,
    na = c("", "NA", "N/A")
  )
  
  
  required <- c(
    coprarna_genome_column,
    coprarna_fdr_column,
    coprarna_combined_p_column,
    coprarna_annotation_column
  )
  
  
  check_columns(
    raw,
    required,
    "CopraRNA"
  )
  
  
  raw <- raw %>%
    filter(
      !is.na(
        .data[[coprarna_genome_column]]
      )
    )
  
  
  parsed <- parse_coprna_reference_field(
    raw[[coprarna_genome_column]]
  )
  
  
  result <- bind_cols(
    raw,
    parsed
  ) %>%
    transmute(
      
      method =
        "CopraRNA",
      
      rank =
        row_number(),
      
      database =
        "CopraRNA",
      
      species =
        NA_character_,
      
      reference_locus =
        reference_locus,
      
      source_locus =
        source_locus,
      
      synonym =
        NA_character_,
      
      prediction_gene =
        prediction_gene,
      
      annotation =
        clean_annotation(
          .data[[coprarna_annotation_column]]
        ),
      
      in_silico_source =
        "CopraRNA",
      
      probability =
        NA_real_,
      
      p_value =
        p_value,
      
      fdr =
        suppressWarnings(
          as.numeric(
            .data[[coprarna_fdr_column]]
          )
        ),
      
      combined_p_value =
        suppressWarnings(
          as.numeric(
            .data[[coprarna_combined_p_column]]
          )
        ),
      
      energy =
        energy,
      
      # CopraRNA reports mRNA positions as 1-based coordinates
      # within the supplied target window. Convert them here to
      # positions relative to the CDS start so that all three
      # prediction methods use the same coordinate system.
      #
      # With 200 nt upstream:
      #   raw position 201 -> CDS-relative position 0
      #   raw position 200 -> -1
      #   raw position 202 -> +1
      mrna_start =
        mrna_start -
        (coprarna_upstream_window_nt + 1L),
      
      mrna_stop =
        mrna_stop -
        (coprarna_upstream_window_nt + 1L),
      
      srna_start =
        srna_start,
      
      srna_stop =
        srna_stop,
      
      source_file =
        basename(file)
    )
  
  
  apply_max_filter(
    result,
    "p_value",
    coprarna_p_max
  )
}



# ============================================================
# IMPORT PREDICTION DATA
# ============================================================

targetrna <- parse_targetrna2(
  targetrna_file
)


intarna <- parse_intarna(
  intarna_file
)


coprarna <- parse_coprna(
  coprarna_file
)



# ============================================================
# COMBINE METHOD-LEVEL PREDICTIONS
# ============================================================
#
# Repeated loci are deliberately retained here.
#
# A locus predicted by three methods therefore occurs in three
# rows in this table.

method_level_all <- bind_rows(
  targetrna,
  intarna,
  coprarna
) %>%
  filter(
    !is.na(reference_locus)
  )



# ============================================================
# ADD REFERENCE ANNOTATION
# ============================================================

method_level_all <- method_level_all %>%
  left_join(
    reference_annotation,
    by = "reference_locus"
  )



# ============================================================
# RESOLVE GENE NAME AND ANNOTATION BY LOCUS
# ============================================================
#
# Prediction outputs do not always contain the same level of
# annotation.
#
# A gene name available from one prediction method or from the
# reference annotation may therefore be propagated to other rows
# corresponding to the same reference locus.
#
# This changes metadata only.
# It does NOT change prediction-method membership.

locus_metadata <- method_level_all %>%
  group_by(
    reference_locus
  ) %>%
  summarise(
    
    resolved_gene =
      most_common_non_missing(
        c(
          prediction_gene,
          reference_gene
        )
      ),
    
    resolved_annotation =
      most_common_non_missing(
        c(
          annotation,
          reference_annotation
        )
      ),
    
    .groups =
      "drop"
  )



method_level_all <- method_level_all %>%
  left_join(
    locus_metadata,
    by = "reference_locus"
  ) %>%
  mutate(
    
    gene =
      coalesce(
        resolved_gene,
        prediction_gene,
        reference_gene
      ),
    
    annotation =
      coalesce(
        resolved_annotation,
        annotation,
        reference_annotation
      )
    
  ) %>%
  select(
    method,
    rank,
    database,
    species,
    reference_locus,
    source_locus,
    synonym,
    prediction_gene,
    gene,
    annotation,
    in_silico_source,
    probability,
    p_value,
    fdr,
    combined_p_value,
    energy,
    mrna_start,
    mrna_stop,
    srna_start,
    srna_stop,
    source_file
  ) %>%
  arrange(
    reference_locus,
    method,
    rank
  )



# ============================================================
# CREATE ANNOTATED METHOD-LEVEL SUBSET
# ============================================================

method_level_annotated <- method_level_all %>%
  filter(
    has_usable_gene_name(gene)
  )



# ============================================================
# FUNCTION: BUILD UNIQUE CANDIDATE TABLE
# ============================================================

build_candidate_table <- function(method_data) {
  
  method_order <- c(
    "CopraRNA",
    "TargetRNA2",
    "IntaRNA"
  )
  
  
  method_data %>%
    group_by(
      reference_locus
    ) %>%
    summarise(
      
      gene =
        first_non_missing(
          gene
        ),
      
      annotation =
        first_non_missing(
          annotation
        ),
      
      CopraRNA =
        any(
          method == "CopraRNA"
        ),
      
      TargetRNA2 =
        any(
          method == "TargetRNA2"
        ),
      
      IntaRNA =
        any(
          method == "IntaRNA"
        ),
      
      n_methods =
        n_distinct(
          method
        ),
      
      prediction_support =
        paste(
          method_order[
            method_order %in%
              unique(method)
          ],
          collapse = "; "
        ),
      
      .groups =
        "drop"
    ) %>%
    arrange(
      desc(n_methods),
      gene,
      reference_locus
    ) %>%
    mutate(
      candidate_id =
        row_number()
    ) %>%
    relocate(
      candidate_id
    )
}



# ============================================================
# BUILD COMPLETE AND ANNOTATED CANDIDATE TABLES
# ============================================================

candidate_all <- build_candidate_table(
  method_level_all
)


candidate_annotated <- build_candidate_table(
  method_level_annotated
)



# ============================================================
# FUNCTION: BUILD UPSET-READY TABLE
# ============================================================

build_overlap_table <- function(candidate_data) {
  
  candidate_data %>%
    select(
      candidate_id,
      reference_locus,
      gene,
      annotation,
      CopraRNA,
      TargetRNA2,
      IntaRNA,
      n_methods
    )
}



overlap_all <- build_overlap_table(
  candidate_all
)


overlap_annotated <- build_overlap_table(
  candidate_annotated
)



# ============================================================
# QC SUMMARY
# ============================================================

method_counts_all <- method_level_all %>%
  distinct(
    method,
    reference_locus
  ) %>%
  count(
    method,
    name = "n_unique_loci"
  ) %>%
  mutate(
    dataset = "all"
  )



method_counts_annotated <- method_level_annotated %>%
  distinct(
    method,
    reference_locus
  ) %>%
  count(
    method,
    name = "n_unique_loci"
  ) %>%
  mutate(
    dataset = "annotated"
  )



qc_summary <- bind_rows(
  
  method_counts_all,
  
  method_counts_annotated,
  
  tibble(
    dataset = "all",
    method = "Combined",
    n_unique_loci = nrow(candidate_all)
  ),
  
  tibble(
    dataset = "annotated",
    method = "Combined",
    n_unique_loci = nrow(candidate_annotated)
  )
  
) %>%
  select(
    dataset,
    method,
    n_unique_loci
  )



# ============================================================
# COORDINATE-SYSTEM QC
# ============================================================
#
# These summaries provide a quick check that the interaction
# coordinates from all prediction tools are now expressed in a
# common CDS-relative coordinate system.
#
# A negative coordinate lies upstream of the CDS start, zero is
# the CDS start, and a positive coordinate lies within the CDS.

coordinate_qc <- method_level_annotated %>%
  mutate(
    binding_left = pmin(
      mrna_start,
      mrna_stop,
      na.rm = FALSE
    ),
    binding_right = pmax(
      mrna_start,
      mrna_stop,
      na.rm = FALSE
    )
  ) %>%
  group_by(
    method
  ) %>%
  summarise(
    n_interactions = n(),
    min_mrna_coordinate = suppressWarnings(
      min(binding_left, na.rm = TRUE)
    ),
    max_mrna_coordinate = suppressWarnings(
      max(binding_right, na.rm = TRUE)
    ),
    n_entirely_upstream = sum(
      binding_right < 0,
      na.rm = TRUE
    ),
    n_overlapping_cds_start = sum(
      binding_left <= 0 &
        binding_right >= 0,
      na.rm = TRUE
    ),
    n_entirely_downstream = sum(
      binding_left > 0,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(
    method
  )



# ============================================================
# OUTPUT FILE PATHS
# ============================================================

method_all_file <- file.path(
  results_dir,
  "target_predictions_method_level_all.csv"
)


candidate_all_file <- file.path(
  results_dir,
  "target_candidates_all.csv"
)


overlap_all_file <- file.path(
  results_dir,
  "target_prediction_overlap_all.csv"
)



method_annotated_file <- file.path(
  results_dir,
  "target_predictions_method_level_annotated.csv"
)


candidate_annotated_file <- file.path(
  results_dir,
  "target_candidates_annotated.csv"
)


overlap_annotated_file <- file.path(
  results_dir,
  "target_prediction_overlap_annotated.csv"
)


qc_file <- file.path(
  results_dir,
  "target_predictions_qc.csv"
)


coordinate_qc_file <- file.path(
  results_dir,
  "target_prediction_coordinate_qc.csv"
)



# ============================================================
# EXPORT COMPLETE SET
# ============================================================

write_csv(
  method_level_all,
  method_all_file,
  na = ""
)


write_csv(
  candidate_all,
  candidate_all_file,
  na = ""
)


write_csv(
  overlap_all,
  overlap_all_file,
  na = ""
)



# ============================================================
# EXPORT ANNOTATED SUBSET
# ============================================================

write_csv(
  method_level_annotated,
  method_annotated_file,
  na = ""
)


write_csv(
  candidate_annotated,
  candidate_annotated_file,
  na = ""
)


write_csv(
  overlap_annotated,
  overlap_annotated_file,
  na = ""
)



# ============================================================
# EXPORT QC
# ============================================================

write_csv(
  qc_summary,
  qc_file,
  na = ""
)


write_csv(
  coordinate_qc,
  coordinate_qc_file,
  na = ""
)



# ============================================================
# CONSOLE SUMMARY
# ============================================================

cat(
  "\n============================================\n",
  "Target-prediction harmonisation complete\n",
  "============================================\n\n",
  sep = ""
)



cat(
  "Input files\n",
  "-----------\n",
  "TargetRNA2: ",
  basename(targetrna_file),
  "\n",
  "IntaRNA:    ",
  basename(intarna_file),
  "\n",
  "CopraRNA:   ",
  basename(coprarna_file),
  "\n\n",
  sep = ""
)



# ------------------------------------------------------------
# Complete set
# ------------------------------------------------------------

cat(
  "COMPLETE HARMONISED SET\n",
  "-----------------------\n"
)


print(
  method_counts_all
)


cat(
  "\nTotal unique candidate loci: ",
  nrow(candidate_all),
  "\n\n",
  sep = ""
)



cat(
  "Exact overlap combinations:\n"
)


print(
  candidate_all %>%
    count(
      CopraRNA,
      TargetRNA2,
      IntaRNA,
      name = "n_candidates"
    ) %>%
    arrange(
      desc(n_candidates)
    )
)



cat(
  "\nTargets predicted by all three methods:\n"
)


print(
  candidate_all %>%
    filter(
      CopraRNA,
      TargetRNA2,
      IntaRNA
    ) %>%
    select(
      reference_locus,
      gene,
      annotation
    ),
  n = Inf
)



# ------------------------------------------------------------
# Annotated set
# ------------------------------------------------------------

cat(
  "\n\nANNOTATED SUBSET\n",
  "----------------\n"
)


print(
  method_counts_annotated
)


cat(
  "\nTotal unique annotated candidate loci: ",
  nrow(candidate_annotated),
  "\n\n",
  sep = ""
)



cat(
  "Exact overlap combinations:\n"
)


print(
  candidate_annotated %>%
    count(
      CopraRNA,
      TargetRNA2,
      IntaRNA,
      name = "n_candidates"
    ) %>%
    arrange(
      desc(n_candidates)
    )
)



cat(
  "\nTargets predicted by all three methods:\n"
)


print(
  candidate_annotated %>%
    filter(
      CopraRNA,
      TargetRNA2,
      IntaRNA
    ) %>%
    select(
      reference_locus,
      gene,
      annotation
    ),
  n = Inf
)



# ------------------------------------------------------------
# Coordinate-system QC
# ------------------------------------------------------------

cat(
  "\n\nCDS-relative interaction-coordinate QC:\n"
)


print(
  coordinate_qc
)


cat(
  "\nCopraRNA coordinates were converted using ",
  coprarna_upstream_window_nt,
  " nt upstream of the CDS start ",
  "(raw window position ",
  coprarna_upstream_window_nt + 1L,
  " = CDS-relative 0).\n",
  sep = ""
)



# ------------------------------------------------------------
# Completion
# ------------------------------------------------------------

cat(
  "\n\nFiles written to:\n",
  results_dir,
  "\n",
  sep = ""
)

