# ============================================================
# 12_export_appendix_table.R
#
# Aim
# ---
# Generate two candidate-target tables:
#
# 1. A compact manuscript/appendix table containing the
#    information required to interpret and prioritise the
#    96 candidate sRNA targets.
#
# 2. A detailed reproducible prioritisation table for the
#    GitHub repository containing sequence, interaction,
#    transcript, operon, motif, functional and STRING context.
#
#
# IMPORTANT
# ---------
#
# No numerical "priority score" is calculated.
#
# The purpose of the detailed table is to expose the individual
# prioritisation features so that candidate selection remains
# transparent rather than reducing heterogeneous biological
# evidence to an arbitrary score.
#
#
# INPUTS
# ------
#
# results/tables/
#   target_candidates_annotated.csv
#   transcript_context.csv
#   target_interactions_transcript_context.csv
#   motif_occurrences.csv
#   STRING_network_nodes.csv
#
#
# OUTPUTS
# -------
#
# results/tables/
#   candidate_prioritisation_appendix.csv
#   candidate_prioritisation_detailed.csv
#   candidate_prioritisation_qc.csv
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
# Input files
# ------------------------------------------------------------

candidate_file <-
  "results/tables/target_candidates_annotated.csv"


transcript_file <-
  "results/tables/transcript_context.csv"


interaction_file <-
  "results/tables/target_interactions_transcript_context.csv"


motif_file <-
  "results/tables/motif_occurrences.csv"


string_nodes_file <-
  "results/tables/STRING_network_nodes.csv"



# ------------------------------------------------------------
# Output files
# ------------------------------------------------------------

appendix_output_file <-
  "results/tables/candidate_prioritisation_appendix.csv"


detailed_output_file <-
  "results/tables/candidate_prioritisation_detailed.csv"


qc_output_file <-
  "results/tables/candidate_prioritisation_qc.csv"



# ------------------------------------------------------------
# sRNA seed regions
# ------------------------------------------------------------
#
# APP GcvB coordinates used throughout the analysis.
#
# R1 = 78-90
# R2 = 137-148

R1_start <- 78L
R1_end   <- 90L

R2_start <- 137L
R2_end   <- 148L



# ============================================================
# PACKAGES
# ============================================================

required_pkgs <- c(
  "readr",
  "dplyr",
  "tidyr",
  "stringr",
  "purrr",
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
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  
})



# ============================================================
# HELPERS
# ============================================================

check_file <- function(file) {
  
  if (!file.exists(file)) {
    
    stop(
      "Required input file not found:\n",
      file
    )
  }
}



check_columns <- function(data, required, label) {
  
  missing <- setdiff(
    required,
    names(data)
  )
  
  
  if (length(missing) > 0) {
    
    stop(
      label,
      " is missing required column(s):\n",
      paste(
        missing,
        collapse = ", "
      )
    )
  }
}



collapse_unique <- function(x, separator = "; ") {
  
  x <- as.character(x)
  
  x <- x[
    !is.na(x) &
      stringr::str_trim(x) != ""
  ]
  
  x <- unique(x)
  
  
  if (length(x) == 0) {
    
    return(
      NA_character_
    )
  }
  
  
  paste(
    x,
    collapse = separator
  )
}



intervals_overlap <- function(
    left1,
    right1,
    left2,
    right2
) {
  
  !is.na(left1) &
    !is.na(right1) &
    left1 <= right2 &
    right1 >= left2
}



interval_distance <- function(
    left1,
    right1,
    left2,
    right2
) {
  
  if (
    any(
      is.na(
        c(
          left1,
          right1,
          left2,
          right2
        )
      )
    )
  ) {
    
    return(
      NA_integer_
    )
  }
  
  
  if (
    left1 <= right2 &&
    right1 >= left2
  ) {
    
    return(
      0L
    )
  }
  
  
  if (
    right1 < left2
  ) {
    
    return(
      as.integer(
        left2 - right1
      )
    )
  }
  
  
  as.integer(
    left1 - right2
  )
}



# ============================================================
# CHECK INPUT FILES
# ============================================================

check_file(
  candidate_file
)

check_file(
  transcript_file
)

check_file(
  interaction_file
)

check_file(
  motif_file
)

check_file(
  string_nodes_file
)



# ============================================================
# READ INPUTS
# ============================================================

candidates <- readr::read_csv(
  candidate_file,
  show_col_types = FALSE
)


transcript_context <- readr::read_csv(
  transcript_file,
  show_col_types = FALSE
)


interactions <- readr::read_csv(
  interaction_file,
  show_col_types = FALSE
)


motif_occurrences <- readr::read_csv(
  motif_file,
  show_col_types = FALSE
)


string_nodes <- readr::read_csv(
  string_nodes_file,
  show_col_types = FALSE
)



# ============================================================
# CHECK CANDIDATE TABLE
# ============================================================

check_columns(
  
  candidates,
  
  c(
    "reference_locus",
    "gene",
    "annotation",
    "CopraRNA",
    "IntaRNA",
    "TargetRNA2"
  ),
  
  "Candidate table"
)



# Preserve the existing curated candidate order.
#
# This avoids creating an artificial numerical prioritisation
# score purely for table ordering.

candidates <- candidates |>
  dplyr::mutate(
    
    candidate_order =
      dplyr::row_number(),
    
    CopraRNA =
      as.logical(
        CopraRNA
      ),
    
    IntaRNA =
      as.logical(
        IntaRNA
      ),
    
    TargetRNA2 =
      as.logical(
        TargetRNA2
      ),
    
    
    n_prediction_methods =
      as.integer(
        CopraRNA
      ) +
      as.integer(
        IntaRNA
      ) +
      as.integer(
        TargetRNA2
      ),
    
    
    prediction_support =
      purrr::pmap_chr(
        
        list(
          CopraRNA,
          IntaRNA,
          TargetRNA2
        ),
        
        function(C, I, T) {
          
          methods <- character(0)
          
          
          if (isTRUE(C)) {
            methods <- c(
              methods,
              "C"
            )
          }
          
          
          if (isTRUE(I)) {
            methods <- c(
              methods,
              "I"
            )
          }
          
          
          if (isTRUE(T)) {
            methods <- c(
              methods,
              "T"
            )
          }
          
          
          paste(
            methods,
            collapse = ", "
          )
        }
      )
  )



# ============================================================
# CHECK TRANSCRIPT CONTEXT TABLE
# ============================================================

check_columns(
  
  transcript_context,
  
  c(
    "reference_locus",
    "TSS_support",
    "representative_TSS",
    "representative_TSS_rel_start",
    "local_UTR5_length",
    "operon_aware_UTR5_length",
    "operon_position",
    "operon_role",
    "transcript_context_class",
    "UTR_interpretation"
  ),
  
  "Transcript-context table"
)



# ============================================================
# CHECK METHOD-LEVEL INTERACTION TABLE
# ============================================================

check_columns(
  
  interactions,
  
  c(
    "reference_locus",
    "method",
    "binding_left",
    "binding_right",
    "srna_binding_left",
    "srna_binding_right",
    "binding_context"
  ),
  
  "Method-level transcript-context table"
)



# ============================================================
# SEED-REGION CLASSIFICATION
# ============================================================

interactions <- interactions |>
  dplyr::mutate(
    
    overlaps_R1 =
      intervals_overlap(
        srna_binding_left,
        srna_binding_right,
        R1_start,
        R1_end
      ),
    
    
    overlaps_R2 =
      intervals_overlap(
        srna_binding_left,
        srna_binding_right,
        R2_start,
        R2_end
      ),
    
    
    seed_region =
      dplyr::case_when(
        
        overlaps_R1 &
          overlaps_R2 ~
          
          "R1; R2",
        
        
        overlaps_R1 ~
          
          "R1",
        
        
        overlaps_R2 ~
          
          "R2",
        
        
        TRUE ~
          
          "Outside R1/R2"
      )
  )



# ============================================================
# SUMMARISE INTERACTIONS BY CANDIDATE
# ============================================================

interaction_summary <- interactions |>
  dplyr::group_by(
    reference_locus
  ) |>
  
  dplyr::summarise(
    
    prediction_methods =
      collapse_unique(
        method
      ),
    
    
    n_method_level_interactions =
      dplyr::n(),
    
    
    seed_regions =
      collapse_unique(
        seed_region
      ),
    
    
    any_R1_binding =
      any(
        overlaps_R1,
        na.rm = TRUE
      ),
    
    
    any_R2_binding =
      any(
        overlaps_R2,
        na.rm = TRUE
      ),
    
    
    any_conserved_seed_binding =
      any(
        overlaps_R1 |
          overlaps_R2,
        na.rm = TRUE
      ),
    
    
    binding_contexts =
      collapse_unique(
        binding_context
      ),
    
    
    any_5UTR_binding =
      any(
        binding_context ==
          "5UTR",
        na.rm = TRUE
      ),
    
    
    any_candidate_internal_5UTR_binding =
      any(
        binding_context ==
          "candidate_internal_5UTR",
        na.rm = TRUE
      ),
    
    
    any_start_codon_overlap =
      any(
        binding_context ==
          "overlaps_start_codon",
        na.rm = TRUE
      ),
    
    
    any_CDS_binding =
      any(
        binding_context ==
          "CDS",
        na.rm = TRUE
      ),
    
    
    minimum_binding_start_rel_CDS =
      suppressWarnings(
        min(
          binding_left,
          na.rm = TRUE
        )
      ),
    
    
    maximum_binding_end_rel_CDS =
      suppressWarnings(
        max(
          binding_right,
          na.rm = TRUE
        )
      ),
    
    
    .groups =
      "drop"
  ) |>
  
  dplyr::mutate(
    
    minimum_binding_start_rel_CDS =
      dplyr::if_else(
        is.infinite(
          minimum_binding_start_rel_CDS
        ),
        NA_real_,
        minimum_binding_start_rel_CDS
      ),
    
    
    maximum_binding_end_rel_CDS =
      dplyr::if_else(
        is.infinite(
          maximum_binding_end_rel_CDS
        ),
        NA_real_,
        maximum_binding_end_rel_CDS
      )
  )



# ============================================================
# APPENDIX INTERACTION-CONTEXT LABEL
# ============================================================
#
# The appendix intentionally uses broad, readable categories.
#
# Detailed classifications remain in the repository table.

appendix_interaction_context <- interactions |>
  dplyr::mutate(
    
    appendix_context =
      dplyr::case_when(
        
        binding_context %in%
          c(
            "5UTR",
            "candidate_internal_5UTR",
            "upstream_unresolved",
            "upstream_of_selected_TSS"
          ) ~
          
          "5' leader / upstream",
        
        
        binding_context ==
          "overlaps_start_codon" ~
          
          "Overlaps start codon",
        
        
        binding_context ==
          "CDS" ~
          
          "CDS",
        
        
        stringr::str_detect(
          binding_context,
          "TSS"
        ) ~
          
          "5' leader / upstream",
        
        
        TRUE ~
          
          binding_context
      )
  ) |>
  
  dplyr::group_by(
    reference_locus
  ) |>
  
  dplyr::summarise(
    
    interaction_context =
      collapse_unique(
        appendix_context
      ),
    
    .groups =
      "drop"
  )



# ============================================================
# FORMAT TSS EVIDENCE FOR APPENDIX
# ============================================================

transcript_context <- transcript_context |>
  dplyr::mutate(
    
    operon_role_lower =
      stringr::str_to_lower(
        dplyr::coalesce(
          operon_role,
          ""
        )
      ),
    
    
    internal_or_last =
      stringr::str_detect(
        operon_role_lower,
        "internal|last"
      ),
    
    
    appendix_TSS_evidence =
      dplyr::case_when(
        
        TSS_support ==
          "Both_concordant" &
          internal_or_last ~
          
          paste0(
            "Concordant proximal signal ",
            "(internal/last operon gene)"
          ),
        
        
        TSS_support ==
          "WT_plate_only" &
          internal_or_last ~
          
          paste0(
            "WT plate proximal signal ",
            "(internal/last operon gene)"
          ),
        
        
        TSS_support ==
          "WT_broth_only" &
          internal_or_last ~
          
          paste0(
            "WT broth proximal signal ",
            "(internal/last operon gene)"
          ),
        
        
        TSS_support ==
          "Both_discordant" &
          internal_or_last ~
          
          paste0(
            "WT plate + broth discordant signal ",
            "(internal/last operon gene)"
          ),
        
        
        TSS_support ==
          "No_proximal_TSS" &
          internal_or_last ~
          
          "No independent proximal TSS detected",
        
        
        TSS_support ==
          "Both_concordant" ~
          
          "WT plate + broth (concordant)",
        
        
        TSS_support ==
          "Both_discordant" ~
          
          "WT plate + broth (discordant)",
        
        
        TSS_support ==
          "WT_plate_only" ~
          
          "WT plate only",
        
        
        TSS_support ==
          "WT_broth_only" ~
          
          "WT broth only",
        
        
        TSS_support ==
          "No_proximal_TSS" ~
          
          "No proximal TSS",
        
        
        TRUE ~
          
          TSS_support
      )
  )



# ============================================================
# FORMAT OPERON POSITION
# ============================================================

transcript_context <- transcript_context |>
  dplyr::mutate(
    
    appendix_operon_position =
      dplyr::case_when(
        
        stringr::str_detect(
          operon_role_lower,
          "single"
        ) ~
          
          "Single-gene",
        
        
        stringr::str_detect(
          operon_role_lower,
          "first"
        ) ~
          
          "First in operon",
        
        
        stringr::str_detect(
          operon_role_lower,
          "internal"
        ) ~
          
          "Internal in operon",
        
        
        stringr::str_detect(
          operon_role_lower,
          "last"
        ) ~
          
          "Last in operon",
        
        
        TRUE ~
          
          as.character(
            operon_position
          )
      ),
    
    
    operon_favourable_for_gene_level_UTR =
      stringr::str_detect(
        operon_role_lower,
        "single|first"
      ),
    
    
    supported_gene_level_UTR =
      !is.na(
        operon_aware_UTR5_length
      )
  )



# ============================================================
# MOTIF OCCURRENCES
# ============================================================

check_columns(
  
  motif_occurrences,
  
  c(
    "dataset",
    "reference_locus",
    "matched_sequence",
    "n_mismatches",
    "motif_start_rel_CDS",
    "motif_end_rel_CDS"
  ),
  
  "Motif occurrence table"
)



candidate_motifs <- motif_occurrences |>
  dplyr::filter(
    dataset ==
      "candidate"
  ) |>
  
  dplyr::mutate(
    
    exact_CA_motif =
      n_mismatches ==
      0L,
    
    
    one_mismatch_CA_motif =
      n_mismatches ==
      1L
  )



motif_summary <- candidate_motifs |>
  dplyr::group_by(
    reference_locus
  ) |>
  
  dplyr::summarise(
    
    CA_motif_present_le1 =
      any(
        n_mismatches <= 1,
        na.rm = TRUE
      ),
    
    
    exact_CA_motif_present =
      any(
        n_mismatches == 0,
        na.rm = TRUE
      ),
    
    
    one_mismatch_CA_motif_present =
      any(
        n_mismatches == 1,
        na.rm = TRUE
      ),
    
    
    best_CA_motif_mismatches =
      min(
        n_mismatches,
        na.rm = TRUE
      ),
    
    
    CA_motif_sequences =
      collapse_unique(
        matched_sequence
      ),
    
    
    CA_motif_start_rel_CDS =
      collapse_unique(
        motif_start_rel_CDS
      ),
    
    
    CA_motif_end_rel_CDS =
      collapse_unique(
        motif_end_rel_CDS
      ),
    
    
    .groups =
      "drop"
  )



# ============================================================
# DISTANCE BETWEEN MOTIF AND PREDICTED mRNA INTERACTION
# ============================================================
#
# For motif-containing candidates, calculate the smallest
# interval distance between any <=1-mismatch motif occurrence
# and any method-specific predicted mRNA interaction.
#
# 0 = intervals overlap.

motif_interaction_pairs <- candidate_motifs |>
  dplyr::filter(
    n_mismatches <=
      1
  ) |>
  
  dplyr::select(
    
    reference_locus,
    
    matched_sequence,
    
    n_mismatches,
    
    motif_start_rel_CDS,
    
    motif_end_rel_CDS
  ) |>
  
  dplyr::inner_join(
    
    interactions |>
      dplyr::select(
        
        reference_locus,
        
        method,
        
        binding_left,
        
        binding_right
      ),
    
    by =
      "reference_locus"
  ) |>
  
  dplyr::rowwise() |>
  
  dplyr::mutate(
    
    motif_interaction_distance_nt =
      interval_distance(
        
        motif_start_rel_CDS,
        motif_end_rel_CDS,
        
        binding_left,
        binding_right
      )
  ) |>
  
  dplyr::ungroup()



motif_distance_summary <- motif_interaction_pairs |>
  dplyr::group_by(
    reference_locus
  ) |>
  
  dplyr::summarise(
    
    minimum_motif_interaction_distance_nt =
      min(
        motif_interaction_distance_nt,
        na.rm = TRUE
      ),
    
    motif_overlaps_predicted_interaction =
      any(
        motif_interaction_distance_nt ==
          0,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  )



# ============================================================
# STRING CONTEXT
# ============================================================

check_columns(
  
  string_nodes,
  
  c(
    "reference_locus",
    "degree",
    "functional_category"
  ),
  
  "STRING node table"
)


# Add betweenness if it is absent.
# This keeps the downstream table structure consistent.

if (!"betweenness" %in% names(string_nodes)) {
  
  string_nodes$betweenness <- NA_real_
}


string_summary <- string_nodes |>
  dplyr::transmute(
    
    reference_locus,
    
    STRING_mapped =
      TRUE,
    
    STRING_connected =
      degree > 0,
    
    STRING_degree =
      as.numeric(
        degree
      ),
    
    STRING_betweenness =
      as.numeric(
        betweenness
      ),
    
    STRING_functional_category =
      as.character(
        functional_category
      )
  )


# ============================================================
# AMINO-ACID/METABOLIC ANNOTATION FLAG
# ============================================================

amino_acid_pattern <- paste(
  
  c(
    "amino acid",
    "amino-acid",
    "aminotransferase",
    "transaminase",
    "branched-chain",
    "acetolactate",
    "isopropylmalate",
    "histidine",
    "histidinol",
    "tryptophan",
    "anthranilate",
    "serine",
    "phosphoserine",
    "threonine",
    "methionine",
    "cysteine",
    "lysine",
    "arginine",
    "proline",
    "leucine",
    "isoleucine",
    "valine",
    "diaminopimelate",
    "chorismate",
    "shikimate"
  ),
  
  collapse =
    "|"
)



# ============================================================
# BUILD DETAILED CANDIDATE TABLE
# ============================================================

detailed <- candidates |>

  dplyr::left_join(
    transcript_context,
    by = "reference_locus",
    suffix = c("", "_transcript")
  ) |>

  dplyr::left_join(
    interaction_summary,
    by = "reference_locus",
    suffix = c("", "_interaction")
  ) |>

  dplyr::left_join(
    appendix_interaction_context,
    by = "reference_locus"
  ) |>

  dplyr::left_join(
    motif_summary,
    by = "reference_locus"
  ) |>

  dplyr::left_join(
    motif_distance_summary,
    by = "reference_locus"
  ) |>

  dplyr::left_join(
    string_summary,
    by = "reference_locus"
  )


# ============================================================
# NORMALISE JOINED INTERACTION COLUMNS
# ============================================================
#
# Some Script 06 outputs may already contain summary fields
# with the same names as those generated above.
#
# left_join() therefore adds suffixes. Here we explicitly
# prefer the interaction-summary version generated in this
# script.

normalise_joined_column <- function(
  data,
  canonical_name,
  preferred_suffix = "_interaction"
) {

  preferred_name <- paste0(
    canonical_name,
    preferred_suffix
  )

  transcript_name <- paste0(
    canonical_name,
    "_transcript"
  )


  if (preferred_name %in% names(data)) {

    data[[canonical_name]] <-
      data[[preferred_name]]

    return(data)
  }


  if (canonical_name %in% names(data)) {

    return(data)
  }


  if (transcript_name %in% names(data)) {

    data[[canonical_name]] <-
      data[[transcript_name]]

    return(data)
  }


  stop(
    "Could not find expected column: ",
    canonical_name
  )
}



interaction_columns_to_normalise <- c(

  "prediction_methods",

  "n_method_level_interactions",

  "seed_regions",

  "any_R1_binding",

  "any_R2_binding",

  "any_conserved_seed_binding",

  "binding_contexts",

  "any_5UTR_binding",

  "any_candidate_internal_5UTR_binding",

  "any_start_codon_overlap",

  "any_CDS_binding",

  "minimum_binding_start_rel_CDS",

  "maximum_binding_end_rel_CDS"
)


for (
  column_name in
  interaction_columns_to_normalise
) {

  detailed <- normalise_joined_column(
    detailed,
    column_name
  )
}



# ============================================================
# ADD PRIORITISATION FEATURES
# ============================================================

detailed <- detailed |>
  dplyr::mutate(

    CA_motif_present_le1 =
      dplyr::coalesce(
        CA_motif_present_le1,
        FALSE
      ),


    exact_CA_motif_present =
      dplyr::coalesce(
        exact_CA_motif_present,
        FALSE
      ),


    one_mismatch_CA_motif_present =
      dplyr::coalesce(
        one_mismatch_CA_motif_present,
        FALSE
      ),


    motif_overlaps_predicted_interaction =
      dplyr::coalesce(
        motif_overlaps_predicted_interaction,
        FALSE
      ),


    STRING_mapped =
      dplyr::coalesce(
        STRING_mapped,
        FALSE
      ),


    STRING_connected =
      dplyr::coalesce(
        STRING_connected,
        FALSE
      ),


    amino_acid_related_annotation =

      stringr::str_detect(

        stringr::str_to_lower(
          dplyr::coalesce(
            annotation,
            ""
          )
        ),

        amino_acid_pattern

      ) |

      dplyr::coalesce(
        STRING_functional_category ==
          "Amino-acid transport and metabolism",
        FALSE
      )
  )



# ============================================================
# HUMAN-READABLE SELECTION RATIONALE
# ============================================================

detailed <- detailed |>
  dplyr::mutate(

    selection_rationale =
      purrr::pmap_chr(

        list(

          n_prediction_methods,

          any_conserved_seed_binding,

          any_5UTR_binding,

          any_start_codon_overlap,

          supported_gene_level_UTR,

          operon_favourable_for_gene_level_UTR,

          CA_motif_present_le1,

          motif_overlaps_predicted_interaction,

          amino_acid_related_annotation,

          STRING_connected
        ),


        function(
          n_methods,
          seed,
          utr_binding,
          start_overlap,
          supported_utr,
          operon_favourable,
          motif,
          motif_overlap,
          aa_related,
          string_connected
        ) {

          evidence <- character(0)


          if (
            !is.na(n_methods) &&
            n_methods >= 2
          ) {

            evidence <- c(
              evidence,
              paste0(
                n_methods,
                " prediction tools"
              )
            )
          }


          if (isTRUE(seed)) {

            evidence <- c(
              evidence,
              "R1/R2 seed interaction"
            )
          }


          if (isTRUE(utr_binding)) {

            evidence <- c(
              evidence,
              "predicted 5' UTR interaction"
            )
          }


          if (isTRUE(start_overlap)) {

            evidence <- c(
              evidence,
              "interaction overlaps translation start"
            )
          }


          if (isTRUE(supported_utr)) {

            evidence <- c(
              evidence,
              "supported gene-level 5' UTR"
            )
          }


          if (isTRUE(operon_favourable)) {

            evidence <- c(
              evidence,
              "single/first operon position"
            )
          }


          if (isTRUE(motif)) {

            evidence <- c(
              evidence,
              "C/A-rich motif <=1 mismatch"
            )
          }


          if (isTRUE(motif_overlap)) {

            evidence <- c(
              evidence,
              "motif overlaps predicted interaction"
            )
          }


          if (isTRUE(aa_related)) {

            evidence <- c(
              evidence,
              "amino-acid/metabolic annotation"
            )
          }


          if (isTRUE(string_connected)) {

            evidence <- c(
              evidence,
              "STRING-connected"
            )
          }


          if (length(evidence) == 0) {

            return(
              NA_character_
            )
          }


          paste(
            evidence,
            collapse = "; "
          )
        }
      )
  )

# ============================================================
# SELECT AND ORDER DETAILED GITHUB TABLE
# ============================================================

detailed_output <- detailed |>
  dplyr::transmute(
    
    No =
      candidate_order,
    
    gene,
    
    reference_locus,
    
    annotation,
    
    
    # --------------------------------------------------------
    # Prediction support
    # --------------------------------------------------------
    
    prediction_support,
    
    n_prediction_methods,
    
    CopraRNA,
    
    IntaRNA,
    
    TargetRNA2,
    
    
    # --------------------------------------------------------
    # GcvB seed-region context
    # --------------------------------------------------------
    
    seed_regions,
    
    any_R1_binding,
    
    any_R2_binding,
    
    any_conserved_seed_binding,
    
    
    # --------------------------------------------------------
    # mRNA interaction position
    # --------------------------------------------------------
    
    binding_contexts,
    
    interaction_context,
    
    minimum_binding_start_rel_CDS,
    
    maximum_binding_end_rel_CDS,
    
    any_5UTR_binding,
    
    any_candidate_internal_5UTR_binding,
    
    any_start_codon_overlap,
    
    any_CDS_binding,
    
    
    # --------------------------------------------------------
    # TSS / transcript architecture
    # --------------------------------------------------------
    
    TSS_support,
    
    representative_TSS,
    
    representative_TSS_rel_start,
    
    local_UTR5_length,
    
    operon_aware_UTR5_length,
    
    transcript_context_class,
    
    UTR_interpretation,
    
    appendix_TSS_evidence,
    
    
    # --------------------------------------------------------
    # Operon context
    # --------------------------------------------------------
    
    operon_id,
    
    operon_size,
    
    operon_position,
    
    operon_role,
    
    appendix_operon_position,
    
    operon_favourable_for_gene_level_UTR,
    
    supported_gene_level_UTR,
    
    
    # --------------------------------------------------------
    # C/A-rich motif
    # --------------------------------------------------------
    
    CA_motif_present_le1,
    
    exact_CA_motif_present,
    
    one_mismatch_CA_motif_present,
    
    best_CA_motif_mismatches,
    
    CA_motif_sequences,
    
    CA_motif_start_rel_CDS,
    
    CA_motif_end_rel_CDS,
    
    minimum_motif_interaction_distance_nt,
    
    motif_overlaps_predicted_interaction,
    
    
    # --------------------------------------------------------
    # Functional context
    # --------------------------------------------------------
    
    amino_acid_related_annotation,
    
    STRING_mapped,
    
    STRING_connected,
    
    STRING_degree,
    
    STRING_betweenness,
    
    STRING_functional_category,
    
    
    # --------------------------------------------------------
    # Transparent prioritisation summary
    # --------------------------------------------------------
    
    selection_rationale
  ) |>
  
  dplyr::arrange(
    No
  )



# ============================================================
# BUILD COMPACT APPENDIX TABLE
# ============================================================

appendix_table <- detailed |>
  dplyr::transmute(
    
    No =
      candidate_order,
    
    Gene =
      gene,
    
    `MIDG2331 locus tag` =
      reference_locus,
    
    Annotation =
      annotation,
    
    `Prediction support` =
      prediction_support,
    
    `Interaction context` =
      interaction_context,
    
    `TSS evidence` =
      appendix_TSS_evidence,
    
    `5' UTR length (nt)` =
      operon_aware_UTR5_length,
    
    `Operon position` =
      appendix_operon_position
  ) |>
  
  dplyr::arrange(
    No
  )



# ============================================================
# QC
# ============================================================

qc <- tibble::tibble(
  
  metric = c(
    
    "candidate_rows",
    
    "appendix_rows",
    
    "detailed_rows",
    
    "duplicate_candidate_loci",
    
    "candidates_with_2plus_prediction_methods",
    
    "candidates_with_3_prediction_methods",
    
    "candidates_with_R1_interaction",
    
    "candidates_with_R2_interaction",
    
    "candidates_with_any_R1_R2_interaction",
    
    "candidates_with_defined_gene_level_5UTR",
    
    "candidates_with_predicted_5UTR_binding",
    
    "candidates_single_or_first_in_operon",
    
    "candidates_with_exact_CA_motif",
    
    "candidates_with_CA_motif_le1_mismatch",
    
    "candidates_CA_motif_overlapping_predicted_interaction",
    
    "candidates_amino_acid_related",
    
    "candidates_STRING_mapped",
    
    "candidates_STRING_connected",
    
    "appendix_missing_gene",
    
    "appendix_missing_annotation",
    
    "appendix_missing_interaction_context"
  ),
  
  
  value = c(
    
    nrow(
      candidates
    ),
    
    nrow(
      appendix_table
    ),
    
    nrow(
      detailed_output
    ),
    
    sum(
      duplicated(
        candidates$reference_locus
      )
    ),
    
    sum(
      detailed$n_prediction_methods >=
        2,
      na.rm = TRUE
    ),
    
    sum(
      detailed$n_prediction_methods ==
        3,
      na.rm = TRUE
    ),
    
    sum(
      detailed$any_R1_binding,
      na.rm = TRUE
    ),
    
    sum(
      detailed$any_R2_binding,
      na.rm = TRUE
    ),
    
    sum(
      detailed$any_conserved_seed_binding,
      na.rm = TRUE
    ),
    
    sum(
      detailed$supported_gene_level_UTR,
      na.rm = TRUE
    ),
    
    sum(
      detailed$any_5UTR_binding,
      na.rm = TRUE
    ),
    
    sum(
      detailed$operon_favourable_for_gene_level_UTR,
      na.rm = TRUE
    ),
    
    sum(
      detailed$exact_CA_motif_present,
      na.rm = TRUE
    ),
    
    sum(
      detailed$CA_motif_present_le1,
      na.rm = TRUE
    ),
    
    sum(
      detailed$motif_overlaps_predicted_interaction,
      na.rm = TRUE
    ),
    
    sum(
      detailed$amino_acid_related_annotation,
      na.rm = TRUE
    ),
    
    sum(
      detailed$STRING_mapped,
      na.rm = TRUE
    ),
    
    sum(
      detailed$STRING_connected,
      na.rm = TRUE
    ),
    
    sum(
      is.na(
        appendix_table$Gene
      ) |
        appendix_table$Gene ==
        ""
    ),
    
    sum(
      is.na(
        appendix_table$Annotation
      ) |
        appendix_table$Annotation ==
        ""
    ),
    
    sum(
      is.na(
        appendix_table$
          `Interaction context`
      ) |
        appendix_table$
        `Interaction context` ==
        ""
    )
  )
)



# ============================================================
# WRITE OUTPUTS
# ============================================================

readr::write_csv(
  
  appendix_table,
  
  appendix_output_file,
  
  na = ""
)



readr::write_csv(
  
  detailed_output,
  
  detailed_output_file,
  
  na = ""
)



readr::write_csv(
  
  qc,
  
  qc_output_file,
  
  na = ""
)



# ============================================================
# CONSOLE REPORT
# ============================================================

cat(
  "\n",
  "============================================\n",
  "Candidate prioritisation tables complete\n",
  "============================================\n\n",
  sep = ""
)


print(
  qc,
  n = Inf,
  width = Inf
)



cat(
  "\nAppendix table preview:\n\n"
)


print(
  
  appendix_table |>
    dplyr::slice_head(
      n = 12
    ),
  
  n = 12,
  width = Inf
)



cat(
  "\nDetailed prioritisation preview:\n\n"
)


print(
  
  detailed_output |>
    dplyr::select(
      
      No,
      gene,
      prediction_support,
      seed_regions,
      interaction_context,
      appendix_TSS_evidence,
      operon_aware_UTR5_length,
      appendix_operon_position,
      CA_motif_present_le1,
      minimum_motif_interaction_distance_nt,
      amino_acid_related_annotation,
      STRING_degree,
      selection_rationale
    ) |>
    
    dplyr::slice_head(
      n = 12
    ),
  
  n = 12,
  width = Inf
)



cat(
  "\nOutputs:\n",
  appendix_output_file,
  "\n",
  detailed_output_file,
  "\n",
  qc_output_file,
  "\n",
  sep = ""
)
