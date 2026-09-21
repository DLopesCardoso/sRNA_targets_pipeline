# ============================================================
# 11_plot_string_network.R
#
# Aim
# ---
# Build a STRING functional-association network for the curated
# sRNA candidate target set.
#
#
# IMPORTANT
# ---------
#
# STRING edges represent functional associations and do not
# necessarily indicate direct physical protein-protein
# interactions.
#
# The complete candidate set is submitted to STRING without
# functional preselection.
#
# Functional categories are assigned only AFTER STRING mapping
# and network construction, for visualisation purposes.
#
#
# INPUT
# -----
#
# results/tables/
#   target_candidates_annotated.csv
#
#
# OUTPUTS
# -------
#
# results/figures/
#   STRING_functional_association_network.pdf
#   STRING_functional_association_network.png
#
# results/tables/
#   STRING_identifier_mapping.csv
#   STRING_network_nodes.csv
#   STRING_network_edges.csv
#   STRING_unmapped_targets.csv
#   STRING_network_qc.csv
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


candidate_file <-
  "results/tables/target_candidates_annotated.csv"


figure_dir <-
  "results/figures"


table_dir <-
  "results/tables"


string_cache_dir <-
  "data/cache/STRING"


# ------------------------------------------------------------
# STRING settings
# ------------------------------------------------------------

# STRING species identifier for
# Actinobacillus pleuropneumoniae.

string_species_id <-
  416269


# STRING version used for manuscript analysis.

string_version <-
  "12.0"


# Minimum combined association score.
#
# 400 = medium confidence
# 700 = high confidence

string_score_threshold <-
  400


# ------------------------------------------------------------
# Labelling settings
# ------------------------------------------------------------

# Single-method targets are labelled only if they are highly
# connected.

minimum_degree_for_single_method_label <-
  5L


# Two-method targets are labelled if they connect to at least
# this many other candidate proteins.

minimum_degree_for_two_method_label <-
  2L


# Targets that should always be labelled if present.

genes_to_always_label <- c(
  "ilvC",
  "ilvG"
)


# Reproducible force-directed layout.

network_seed <-
  123L



# ============================================================
# PACKAGES
# ============================================================

cran_packages <- c(
  "readr",
  "dplyr",
  "stringr",
  "tibble",
  "purrr",
  "igraph",
  "tidygraph",
  "ggraph",
  "ggplot2"
)


missing_cran <- cran_packages[
  !vapply(
    cran_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]


if (length(missing_cran) > 0) {
  
  install.packages(
    missing_cran
  )
}


if (
  !requireNamespace(
    "STRINGdb",
    quietly = TRUE
  )
) {
  
  if (
    !requireNamespace(
      "BiocManager",
      quietly = TRUE
    )
  ) {
    
    install.packages(
      "BiocManager"
    )
  }
  
  
  BiocManager::install(
    "STRINGdb",
    ask = FALSE,
    update = FALSE
  )
}


suppressPackageStartupMessages({
  
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(ggplot2)
  library(STRINGdb)
  
})



# ============================================================
# CREATE OUTPUT DIRECTORIES
# ============================================================

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  string_cache_dir,
  recursive = TRUE,
  showWarnings = FALSE
)



# ============================================================
# CHECK INPUT FILE
# ============================================================

if (!file.exists(candidate_file)) {
  
  stop(
    "Candidate table not found:\n",
    candidate_file,
    "\n\nRun Script 01 first."
  )
}



# ============================================================
# READ CANDIDATE TABLE
# ============================================================

targets <- readr::read_csv(
  candidate_file,
  show_col_types = FALSE
)



# ============================================================
# CHECK REQUIRED COLUMNS
# ============================================================

required_columns <- c(
  
  "reference_locus",
  "gene",
  "annotation",
  
  "CopraRNA",
  "TargetRNA2",
  "IntaRNA",
  
  "n_methods"
)


missing_columns <- setdiff(
  required_columns,
  names(targets)
)


if (length(missing_columns) > 0) {
  
  stop(
    "Candidate table is missing required column(s):\n",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}



# ============================================================
# STANDARDISE TARGET TABLE
# ============================================================

targets <- targets |>
  dplyr::mutate(
    
    reference_locus =
      stringr::str_trim(
        as.character(
          reference_locus
        )
      ),
    
    gene =
      stringr::str_trim(
        as.character(
          gene
        )
      ),
    
    annotation =
      stringr::str_trim(
        as.character(
          annotation
        )
      ),
    
    CopraRNA =
      as.logical(
        CopraRNA
      ),
    
    TargetRNA2 =
      as.logical(
        TargetRNA2
      ),
    
    IntaRNA =
      as.logical(
        IntaRNA
      ),
    
    n_methods =
      suppressWarnings(
        as.integer(
          n_methods
        )
      )
  ) |>
  
  dplyr::filter(
    !is.na(
      reference_locus
    ),
    reference_locus != ""
  ) |>
  
  dplyr::distinct(
    reference_locus,
    .keep_all = TRUE
  )



# ============================================================
# CHECK METHOD-SUPPORT CONSISTENCY
# ============================================================

targets <- targets |>
  dplyr::mutate(
    
    n_methods_recalculated =
      as.integer(
        CopraRNA
      ) +
      as.integer(
        TargetRNA2
      ) +
      as.integer(
        IntaRNA
      )
  )


support_mismatch <- targets |>
  dplyr::filter(
    n_methods !=
      n_methods_recalculated
  )


if (
  nrow(
    support_mismatch
  ) > 0
) {
  
  stop(
    nrow(
      support_mismatch
    ),
    " candidate(s) have inconsistent n_methods values."
  )
}



# ============================================================
# CREATE SUPPORT GROUPS
# ============================================================

targets <- targets |>
  dplyr::mutate(
    
    prediction_support =
      n_methods_recalculated,
    
    
    support_group =
      dplyr::case_when(
        
        prediction_support == 3L ~
          
          "Predicted by 3 tools",
        
        prediction_support == 2L ~
          
          "Predicted by 2 tools",
        
        prediction_support == 1L ~
          
          "Predicted by 1 tool",
        
        TRUE ~
          
          "Other"
      ),
    
    
    support_group =
      factor(
        
        support_group,
        
        levels = c(
          "Predicted by 1 tool",
          "Predicted by 2 tools",
          "Predicted by 3 tools"
        )
      ),
    
    
    display_name =
      dplyr::case_when(
        
        !is.na(gene) &
          gene != "" ~
          
          gene,
        
        TRUE ~
          
          reference_locus
      )
  )



# ============================================================
# QC: EXPECTED CANDIDATE SET
# ============================================================

cat(
  "\nCandidate loci imported: ",
  nrow(targets),
  "\n",
  sep = ""
)


cat(
  "\nPrediction-support distribution:\n"
)


print(
  
  targets |>
    dplyr::count(
      support_group,
      name =
        "n_targets"
    )
)



# ============================================================
# ASSIGN BROAD FUNCTIONAL CATEGORIES
# ============================================================
#
# These are used ONLY for network visualisation.
#
# Categories are assigned from existing annotation text after
# candidate selection.
#
# They do not influence STRING mapping, edge retrieval or
# network inclusion.

targets <- targets |>
  dplyr::mutate(
    
    annotation_lower =
      stringr::str_to_lower(
        dplyr::coalesce(
          annotation,
          ""
        )
      ),
    
    
    functional_category =
      dplyr::case_when(
        
        # ----------------------------------------------------
        # Amino-acid transport and metabolism
        # ----------------------------------------------------
        
        stringr::str_detect(
          
          annotation_lower,
          
          paste(
            
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
              "shikimate",
              "aminopeptidase",
              "peptide transport",
              "oligopeptide",
              "dipeptide"
            ),
            
            collapse = "|"
          )
        ) ~
          
          "Amino-acid transport and metabolism",
        
        
        # ----------------------------------------------------
        # Translation
        # ----------------------------------------------------
        
        stringr::str_detect(
          
          annotation_lower,
          
          paste(
            
            c(
              "ribosomal protein",
              "ribosome",
              "translation",
              "elongation factor",
              "initiation factor",
              "trna",
              "aminoacyl"
            ),
            
            collapse = "|"
          )
        ) ~
          
          "Translation",
        
        
        # ----------------------------------------------------
        # Regulation
        # ----------------------------------------------------
        
        stringr::str_detect(
          
          annotation_lower,
          
          paste(
            
            c(
              "transcriptional regulator",
              "response regulator",
              "two-component",
              "regulatory protein",
              "transcription factor",
              "repressor",
              "activator"
            ),
            
            collapse = "|"
          )
        ) ~
          
          "Regulation",
        
        
        # ----------------------------------------------------
        # Everything else
        # ----------------------------------------------------
        
        TRUE ~
          
          "Other"
      ),
    
    
    functional_category =
      factor(
        
        functional_category,
        
        levels = c(
          "Amino-acid transport and metabolism",
          "Translation",
          "Regulation",
          "Other"
        )
      )
  )



# ============================================================
# INITIALISE STRINGdb
# ============================================================

string_db <- STRINGdb$new(
  
  version =
    string_version,
  
  species =
    string_species_id,
  
  score_threshold =
    string_score_threshold,
  
  input_directory =
    string_cache_dir
)



# ============================================================
# PREPARE GENE-NAME MAPPING INPUT
# ============================================================

gene_mapping_input <- targets |>
  dplyr::filter(
    !is.na(gene),
    gene != ""
  ) |>
  
  dplyr::transmute(
    
    reference_locus =
      as.character(
        reference_locus
      ),
    
    identifier =
      as.character(
        gene
      ),
    
    mapping_source =
      "Gene name"
  ) |>
  
  dplyr::filter(
    !is.na(identifier),
    identifier != ""
  ) |>
  
  dplyr::distinct(
    identifier,
    .keep_all = TRUE
  ) |>
  
  as.data.frame(
    stringsAsFactors = FALSE
  )



# ============================================================
# MAP GENE NAMES TO STRING
# ============================================================

if (
  nrow(
    gene_mapping_input
  ) == 0
) {
  
  mapped_by_gene <-
    data.frame()
  
} else {
  
  mapped_by_gene <- string_db$map(
    
    gene_mapping_input,
    
    "identifier",
    
    takeFirst =
      TRUE,
    
    removeUnmappedRows =
      TRUE,
    
    quiet =
      FALSE
  )
}



# ============================================================
# PREPARE LOCUS-TAG MAPPING INPUT
# ============================================================

locus_mapping_input <- targets |>
  dplyr::transmute(
    
    reference_locus =
      as.character(
        reference_locus
      ),
    
    identifier =
      as.character(
        reference_locus
      ),
    
    mapping_source =
      "Locus tag"
  ) |>
  
  dplyr::filter(
    !is.na(identifier),
    identifier != ""
  ) |>
  
  dplyr::distinct(
    identifier,
    .keep_all = TRUE
  ) |>
  
  as.data.frame(
    stringsAsFactors = FALSE
  )



# ============================================================
# MAP LOCUS TAGS TO STRING
# ============================================================

if (
  nrow(
    locus_mapping_input
  ) == 0
) {
  
  mapped_by_locus <-
    data.frame()
  
} else {
  
  mapped_by_locus <- string_db$map(
    
    locus_mapping_input,
    
    "identifier",
    
    takeFirst =
      TRUE,
    
    removeUnmappedRows =
      TRUE,
    
    quiet =
      FALSE
  )
}



# ============================================================
# SELECT PREFERRED MAPPING PER LOCUS
# ============================================================
#
# Prefer gene-name mapping.
# Use locus-tag mapping as fallback.

mapping_combined <- dplyr::bind_rows(
  
  mapped_by_gene,
  mapped_by_locus
  
) |>
  
  dplyr::mutate(
    
    mapping_priority =
      dplyr::case_when(
        
        mapping_source ==
          "Gene name" ~
          
          1L,
        
        mapping_source ==
          "Locus tag" ~
          
          2L,
        
        TRUE ~
          
          3L
      )
  ) |>
  
  dplyr::arrange(
    reference_locus,
    mapping_priority
  ) |>
  
  dplyr::distinct(
    reference_locus,
    .keep_all = TRUE
  )



# ============================================================
# CHECK STRING MAPPING COLUMN
# ============================================================

if (
  !"STRING_id" %in%
  names(
    mapping_combined
  )
) {
  
  stop(
    "STRINGdb mapping did not return a STRING_id column."
  )
}



# ============================================================
# JOIN MAPPING TO TARGET METADATA
# ============================================================

preferred_name_column <-
  if (
    "preferred_name" %in%
    names(
      mapping_combined
    )
  ) {
    
    "preferred_name"
    
  } else {
    
    NULL
  }


mapping_columns <- c(
  
  "reference_locus",
  "identifier",
  "mapping_source",
  "STRING_id",
  
  preferred_name_column
)


mapped_targets <- mapping_combined |>
  dplyr::select(
    dplyr::all_of(
      mapping_columns
    )
  ) |>
  
  dplyr::left_join(
    targets,
    by =
      "reference_locus"
  )



# ============================================================
# IDENTIFY UNMAPPED TARGETS
# ============================================================

unmapped_targets <- targets |>
  dplyr::anti_join(
    
    mapped_targets |>
      dplyr::select(
        reference_locus
      ),
    
    by =
      "reference_locus"
  )



# ============================================================
# MAPPING SUMMARY
# ============================================================

mapping_percentage <-
  100 *
  nrow(
    mapped_targets
  ) /
  nrow(
    targets
  )


cat(
  "\nTargets mapped to STRING: ",
  nrow(mapped_targets),
  " / ",
  nrow(targets),
  " (",
  round(
    mapping_percentage,
    1
  ),
  "%)\n",
  sep = ""
)


cat(
  "Targets not mapped: ",
  nrow(
    unmapped_targets
  ),
  "\n",
  sep = ""
)



if (
  nrow(
    mapped_targets
  ) < 5
) {
  
  stop(
    "Fewer than five candidate targets mapped to STRING."
  )
}



# ============================================================
# RETRIEVE STRING ASSOCIATIONS
# ============================================================

string_interactions_raw <-
  string_db$get_interactions(
    mapped_targets$STRING_id
  )



# ============================================================
# IDENTIFY STRING SCORE COLUMN
# ============================================================

possible_score_columns <- c(
  "combined_score",
  "score"
)


score_column <-
  possible_score_columns[
    possible_score_columns %in%
      names(
        string_interactions_raw
      )
  ][1]


if (
  length(
    score_column
  ) == 0 ||
  is.na(
    score_column
  )
) {
  
  stop(
    "No recognised STRING association-score column was found."
  )
}



# ============================================================
# CLEAN STRING EDGE TABLE
# ============================================================

# Extract the STRING score explicitly before using dplyr.
# This avoids problems with dynamic column selection.

string_interactions_clean <- string_interactions_raw |>
  dplyr::mutate(
    combined_score =
      as.numeric(
        .data[[score_column]]
      )
  )


# Keep only associations where BOTH proteins are members of
# the mapped candidate-target set.

edge_table <- string_interactions_clean |>
  dplyr::filter(
    from %in% mapped_targets$STRING_id,
    to %in% mapped_targets$STRING_id,
    from != to
  ) |>
  dplyr::transmute(
    from = as.character(from),
    to = as.character(to),
    combined_score = combined_score
  ) |>
  dplyr::mutate(
    
    # STRING associations are undirected.
    # Standardise each pair so A--B and B--A are treated
    # as the same association.
    
    node_1 = pmin(
      from,
      to
    ),
    
    node_2 = pmax(
      from,
      to
    )
  ) |>
  dplyr::arrange(
    dplyr::desc(
      combined_score
    )
  ) |>
  dplyr::distinct(
    node_1,
    node_2,
    .keep_all = TRUE
  ) |>
  dplyr::select(
    from = node_1,
    to = node_2,
    combined_score
  )


cat(
  "\nNumber of STRING associations retained among candidate targets: ",
  nrow(edge_table),
  "\n",
  sep = ""
)


if (
  nrow(
    edge_table
  ) == 0
) {
  
  stop(
    "No STRING associations were found among the mapped ",
    "candidate targets at the selected threshold."
  )
}
  


# ============================================================
# CONSTRUCT NODE TABLE
# ============================================================

if (
  !"preferred_name" %in%
  names(
    mapped_targets
  )
) {
  
  mapped_targets$preferred_name <-
    NA_character_
}



node_table <- mapped_targets |>
  dplyr::mutate(
    
    display_name =
      dplyr::case_when(
        
        !is.na(
          preferred_name
        ) &
          preferred_name != "" ~
          
          preferred_name,
        
        !is.na(gene) &
          gene != "" ~
          
          gene,
        
        TRUE ~
          
          reference_locus
      )
  ) |>
  
  dplyr::transmute(
    
    name =
      STRING_id,
    
    STRING_id =
      STRING_id,
    
    reference_locus =
      reference_locus,
    
    gene =
      gene,
    
    display_name =
      display_name,
    
    annotation =
      annotation,
    
    CopraRNA =
      CopraRNA,
    
    TargetRNA2 =
      TargetRNA2,
    
    IntaRNA =
      IntaRNA,
    
    prediction_support =
      prediction_support,
    
    support_group =
      support_group,
    
    functional_category =
      functional_category,
    
    mapping_source =
      mapping_source
  ) |>
  
  dplyr::distinct(
    name,
    .keep_all =
      TRUE
  )



# ============================================================
# BUILD IGRAPH NETWORK
# ============================================================

network_igraph <- igraph::graph_from_data_frame(
  
  d =
    edge_table,
  
  directed =
    FALSE,
  
  vertices =
    node_table
)



network_igraph <- igraph::simplify(
  
  network_igraph,
  
  remove.multiple =
    TRUE,
  
  remove.loops =
    TRUE,
  
  edge.attr.comb =
    "max"
)



# ============================================================
# CONVERT TO TIDYGRAPH
# ============================================================

network_graph <- tidygraph::as_tbl_graph(
  network_igraph
) |>
  
  tidygraph::activate(
    nodes
  ) |>
  
  dplyr::mutate(
    
    degree =
      tidygraph::centrality_degree(),
    
    betweenness =
      tidygraph::centrality_betweenness(),
    
    connected =
      degree > 0,
    
    always_label =
      stringr::str_to_lower(
        dplyr::coalesce(
          gene,
          display_name
        )
      ) %in%
      stringr::str_to_lower(
        genes_to_always_label
      )
  )



# ============================================================
# RETAIN CONNECTED CANDIDATES FOR FIGURE
# ============================================================

network_graph_connected <- network_graph |>
  
  tidygraph::activate(
    nodes
  ) |>
  
  dplyr::filter(
    degree > 0
  ) |>
  
  dplyr::mutate(
    
    label_network =
      dplyr::case_when(
        
        # Always label explicitly prioritised genes.
        
        always_label ~
          
          display_name,
        
        
        # Label all three-tool predictions.
        
        prediction_support == 3L ~
          
          display_name,
        
        
        # Label reasonably connected two-tool predictions.
        
        prediction_support == 2L &
          degree >=
          minimum_degree_for_two_method_label ~
          
          display_name,
        
        
        # Label highly connected one-tool predictions.
        
        prediction_support == 1L &
          degree >=
          minimum_degree_for_single_method_label ~
          
          display_name,
        
        
        TRUE ~
          
          ""
      ),
    
    
    functional_category =
      droplevels(
        functional_category
      ),
    
    
    support_group =
      droplevels(
        support_group
      )
  )



# ============================================================
# EXPORT NODE TABLE WITH NETWORK STATISTICS
# ============================================================

network_node_output <- network_graph |>
  
  tidygraph::activate(
    nodes
  ) |>
  
  tibble::as_tibble() |>
  
  dplyr::arrange(
    
    dplyr::desc(
      prediction_support
    ),
    
    dplyr::desc(
      degree
    ),
    
    display_name
  )



# ============================================================
# EXPORT EDGE TABLE
# ============================================================

network_edge_output <- network_graph |>
  
  tidygraph::activate(
    edges
  ) |>
  
  tibble::as_tibble() |>
  
  dplyr::arrange(
    dplyr::desc(
      combined_score
    )
  )



# ============================================================
# NETWORK QC SUMMARY
# ============================================================

network_qc <- tibble::tibble(
  
  metric = c(
    
    "candidate_loci",
    
    "targets_mapped_to_STRING",
    
    "targets_not_mapped_to_STRING",
    
    "STRING_mapping_percent",
    
    "mapped_targets_with_at_least_one_association",
    
    "STRING_associations",
    
    "targets_predicted_by_at_least_two_tools",
    
    "targets_predicted_by_all_three_tools",
    
    "connected_amino_acid_nodes",
    
    "connected_translation_nodes",
    
    "connected_regulation_nodes",
    
    "connected_other_nodes"
  ),
  
  
  value = c(
    
    nrow(
      targets
    ),
    
    nrow(
      mapped_targets
    ),
    
    nrow(
      unmapped_targets
    ),
    
    mapping_percentage,
    
    sum(
      network_node_output$
        degree >
        0
    ),
    
    nrow(
      edge_table
    ),
    
    sum(
      targets$
        prediction_support >=
        2
    ),
    
    sum(
      targets$
        prediction_support ==
        3
    ),
    
    sum(
      network_node_output$
        degree >
        0 &
        network_node_output$
        functional_category ==
        "Amino-acid transport and metabolism",
      na.rm = TRUE
    ),
    
    sum(
      network_node_output$
        degree >
        0 &
        network_node_output$
        functional_category ==
        "Translation",
      na.rm = TRUE
    ),
    
    sum(
      network_node_output$
        degree >
        0 &
        network_node_output$
        functional_category ==
        "Regulation",
      na.rm = TRUE
    ),
    
    sum(
      network_node_output$
        degree >
        0 &
        network_node_output$
        functional_category ==
        "Other",
      na.rm = TRUE
    )
  )
)



# ============================================================
# OUTPUT PATHS
# ============================================================

mapping_file <- file.path(
  table_dir,
  "STRING_identifier_mapping.csv"
)


node_file <- file.path(
  table_dir,
  "STRING_network_nodes.csv"
)


edge_file <- file.path(
  table_dir,
  "STRING_network_edges.csv"
)


unmapped_file <- file.path(
  table_dir,
  "STRING_unmapped_targets.csv"
)


qc_file <- file.path(
  table_dir,
  "STRING_network_qc.csv"
)


pdf_file <- file.path(
  figure_dir,
  "STRING_functional_association_network.pdf"
)


png_file <- file.path(
  figure_dir,
  "STRING_functional_association_network.png"
)



# ============================================================
# WRITE TABLES
# ============================================================

readr::write_csv(
  mapped_targets,
  mapping_file,
  na = ""
)


readr::write_csv(
  network_node_output,
  node_file,
  na = ""
)


readr::write_csv(
  network_edge_output,
  edge_file,
  na = ""
)


readr::write_csv(
  unmapped_targets,
  unmapped_file,
  na = ""
)


readr::write_csv(
  network_qc,
  qc_file,
  na = ""
)



# ============================================================
# GENERATE PUBLICATION NETWORK
# ============================================================

set.seed(
  network_seed
)


string_network_plot <- ggraph::ggraph(
  network_graph_connected,
  layout =
    "stress"
) +
  
  
  # ----------------------------------------------------------
# STRING functional associations
# ----------------------------------------------------------

ggraph::geom_edge_link(
  
  ggplot2::aes(
    edge_width =
      combined_score
  ),
  
  edge_colour =
    "grey65",
  
  edge_alpha =
    0.70,
  
  show.legend =
    TRUE
) +
  
  
  ggraph::scale_edge_width(
    
    name =
      "STRING association score",
    
    range =
      c(
        0.30,
        1.8
      ),
    
    breaks =
      c(
        400,
        700,
        900
      ),
    
    labels =
      c(
        "400",
        "700",
        "900"
      )
  ) +
  
  
  # ----------------------------------------------------------
# Nodes
# ----------------------------------------------------------

ggraph::geom_node_point(
  
  ggplot2::aes(
    
    fill =
      functional_category,
    
    shape =
      support_group
  ),
  
  size =
    5.8,
  
  colour =
    "grey20",
  
  stroke =
    0.6
) +
  
  
  # ----------------------------------------------------------
# Prediction-support shapes
# ----------------------------------------------------------

ggplot2::scale_shape_manual(
  
  name =
    "Prediction support",
  
  values =
    c(
      
      "Predicted by 1 tool" =
        21,
      
      "Predicted by 2 tools" =
        22,
      
      "Predicted by 3 tools" =
        24
    ),
  
  labels =
    c(
      
      "Predicted by 1 tool" =
        "1 tool",
      
      "Predicted by 2 tools" =
        "2 tools",
      
      "Predicted by 3 tools" =
        "3 tools"
    ),
  
  drop =
    TRUE
) +
  
  
  # ----------------------------------------------------------
# Functional categories
# ----------------------------------------------------------

ggplot2::scale_fill_brewer(
  
  name =
    "Functional category",
  
  palette =
    "Set2",
  
  drop =
    TRUE,
  
  na.value =
    "grey80"
) +
  
  
  # ----------------------------------------------------------
# Gene labels
# ----------------------------------------------------------

ggraph::geom_node_text(
  
  ggplot2::aes(
    label =
      label_network
  ),
  
  repel =
    TRUE,
  
  size =
    4.0,
  
  family =
    "sans",
  
  fontface =
    "italic",
  
  colour =
    "grey10",
  
  box.padding =
    0.50,
  
  point.padding =
    0.40,
  
  max.overlaps =
    Inf,
  
  segment.colour =
    "grey65",
  
  segment.linewidth =
    0.30
) +
  
  
# ----------------------------------------------------------
# Title
# ----------------------------------------------------------

ggplot2::labs(
  
  title =
    "Functional association network of predicted GcvB targets",
  
  subtitle =
    paste0(
      "STRING associations among connected candidate targets ",
      "(minimum combined score = ",
      string_score_threshold,
      ")"
    ),
  
  caption =
    NULL
) +
  
  
  # ----------------------------------------------------------
# Theme
# ----------------------------------------------------------

ggraph::theme_graph(
  
  base_family =
    "sans",
  
  base_size =
    13
) +
  
  
  ggplot2::theme(
    
    plot.title =
      ggplot2::element_text(
        
        face =
          "bold",
        
        size =
          16,
        
        hjust =
          0
      ),
    
    
    plot.subtitle =
      ggplot2::element_text(
        
        size =
          12,
        
        hjust =
          0,
        
        margin =
          ggplot2::margin(
            b = 10
          )
      ),
    
    
    legend.position =
      "bottom",
    
    legend.direction =
      "vertical",
    
    legend.box =
      "vertical",
    
    legend.box.just =
      "left",
    
    legend.justification =
      "left",
    
    
    legend.spacing.x =
      grid::unit(
        0.35,
        "cm"
      ),
    
    
    legend.spacing.y =
      grid::unit(
        0.10,
        "cm"
      ),
    
    
    legend.title =
      ggplot2::element_text(
        
        face =
          "bold",
        
        size =
          10.5
      ),
    
    
    legend.text =
      ggplot2::element_text(
        size = 10
      ),
    
    
    legend.key.width =
      grid::unit(
        0.65,
        "cm"
      ),
    
    
    legend.key.height =
      grid::unit(
        0.45,
        "cm"
      ),
    
    
    legend.box.background =
      ggplot2::element_rect(
        
        colour =
          "grey65",
        
        fill =
          "white",
        
        linewidth =
          0.5
      ),
    
    
    legend.box.margin =
      ggplot2::margin(
        t = 7,
        r = 10,
        b = 7,
        l = 10
      ),
    
    
    plot.margin =
      ggplot2::margin(
        t = 12,
        r = 20,
        b = 25,
        l = 20
      )
  ) +
  
  
  # ----------------------------------------------------------
# Legend order
# ----------------------------------------------------------

ggplot2::guides(
  
  fill =
    ggplot2::guide_legend(
      
      order =
        1,
      
      nrow =
        1,
      
      byrow =
        TRUE,
      
      override.aes =
        list(
          
          size =
            5.5,
          
          shape =
            21,
          
          colour =
            "grey20"
        )
    ),
  
  
  shape =
    ggplot2::guide_legend(
      
      order =
        2,
      
      nrow =
        1,
      
      override.aes =
        list(
          
          size =
            5.5,
          
          fill =
            "white",
          
          colour =
            "grey20"
        )
    ),
  
  
  edge_width =
    ggplot2::guide_legend(
      
      order =
        3,
      
      nrow =
        1,
      
      title.position =
        "top"
    )
)



# ============================================================
# DISPLAY FIGURE
# ============================================================

print(
  string_network_plot
)



# ============================================================
# SAVE FIGURE
# ============================================================

ggplot2::ggsave(
  
  filename =
    pdf_file,
  
  plot =
    string_network_plot,
  
  width =
    10,
  
  height =
    8,
  
  units =
    "in",
  
  device =
    cairo_pdf,
  
  limitsize =
    FALSE
)



ggplot2::ggsave(
  
  filename =
    png_file,
  
  plot =
    string_network_plot,
  
  width =
    10,
  
  height =
    8,
  
  units =
    "in",
  
  dpi =
    600,
  
  limitsize =
    FALSE
)



# ============================================================
# CONSOLE REPORT
# ============================================================

cat(
  "\n",
  "============================================\n",
  "STRING network analysis complete\n",
  "============================================\n\n",
  sep = ""
)


print(
  network_qc,
  n = Inf,
  width = Inf
)


cat(
  "\nUnmapped STRING targets:\n"
)


if (
  nrow(
    unmapped_targets
  ) == 0
) {
  
  cat(
    "None\n"
  )
  
} else {
  
  print(
    
    unmapped_targets |>
      dplyr::select(
        reference_locus,
        gene,
        annotation,
        prediction_support
      ),
    
    n = Inf,
    width = Inf
  )
}


cat(
  "\nOutputs:\n",
  pdf_file,
  "\n",
  png_file,
  "\n",
  mapping_file,
  "\n",
  node_file,
  "\n",
  edge_file,
  "\n",
  unmapped_file,
  "\n",
  qc_file,
  "\n",
  sep = ""
)

