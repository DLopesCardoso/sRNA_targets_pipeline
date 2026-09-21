# ============================================================
# 10_plot_functional_enrichment.R
#
# Aim
# ---
# Generate the manuscript functional-enrichment figure from
# ShinyGO output for the curated sRNA target set.
#
#
# IMPORTANT
# ---------
#
# Enrichment statistics are NOT recalculated in this script.
#
# Fold enrichment, FDR and gene counts are taken directly from
# the supplied ShinyGO output.
#
# Representative terms are selected only to reduce visual
# redundancy among overlapping metabolic categories.
#
#
# INPUT
# -----
#
# data/raw/functional_enrichment/
#   enrichment_all.csv
#
#
# OUTPUTS
# -------
#
# results/figures/
#   functional_enrichment.pdf
#   functional_enrichment.png
#
# results/tables/
#   functional_enrichment_selected_terms.csv
#   functional_enrichment_qc.csv
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


input_file <-
  "data/raw/functional_enrichment/enrichment_all.csv"


figure_dir <-
  "results/figures"


table_dir <-
  "results/tables"


pdf_file <-
  "results/figures/functional_enrichment.pdf"


png_file <-
  "results/figures/functional_enrichment.png"


selected_output_file <-
  "results/tables/functional_enrichment_selected_terms.csv"


qc_output_file <-
  "results/tables/functional_enrichment_qc.csv"



# ============================================================
# REPRESENTATIVE TERMS FOR MANUSCRIPT FIGURE
# ============================================================
#
# These terms are selected from the complete ShinyGO output
# to reduce visual redundancy.
#
# Their enrichment statistics are preserved exactly as
# reported by ShinyGO.

selected_terms <- c(
  
  "Biosynthesis of amino acids",
  
  "Branched-chain amino acid biosynthesis",
  
  "Histidine biosynthetic process",
  
  "L-serine biosynthetic process",
  
  "Aromatic amino acid biosynthesis",
  
  "Diaminopimelate biosynthetic process",
  
  "Cysteine metabolic process",
  
  "Dicarboxylic acid biosynthetic process",
  
  "Arginine metabolic process, and Lysine biosynthesis"
)



# ============================================================
# PACKAGES
# ============================================================

required_pkgs <- c(
  "readr",
  "dplyr",
  "ggplot2",
  "forcats",
  "viridis",
  "scales",
  "cowplot",
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
  library(ggplot2)
  library(forcats)
  library(viridis)
  library(scales)
  library(cowplot)
  library(tibble)
  
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



# ============================================================
# CHECK INPUT FILE
# ============================================================

if (!file.exists(input_file)) {
  
  stop(
    "ShinyGO input file not found:\n",
    input_file
  )
}



# ============================================================
# READ SHINYGO OUTPUT
# ============================================================

enrichment <- readr::read_csv(
  input_file,
  show_col_types = FALSE
)



# ============================================================
# CHECK REQUIRED COLUMNS
# ============================================================

required_cols <- c(
  
  "Enrichment FDR",
  
  "nGenes",
  
  "Pathway Genes",
  
  "Fold Enrichment",
  
  "Pathway",
  
  "Genes"
)


missing_cols <- setdiff(
  required_cols,
  names(enrichment)
)


if (length(missing_cols) > 0) {
  
  stop(
    "ShinyGO input is missing required column(s):\n",
    paste(
      missing_cols,
      collapse = ", "
    )
  )
}



# ============================================================
# STANDARDISE NUMERIC COLUMNS
# ============================================================

enrichment <- enrichment |>
  dplyr::mutate(
    
    `Enrichment FDR` =
      suppressWarnings(
        as.numeric(
          `Enrichment FDR`
        )
      ),
    
    nGenes =
      suppressWarnings(
        as.integer(
          nGenes
        )
      ),
    
    `Pathway Genes` =
      suppressWarnings(
        as.integer(
          `Pathway Genes`
        )
      ),
    
    `Fold Enrichment` =
      suppressWarnings(
        as.numeric(
          `Fold Enrichment`
        )
      ),
    
    Pathway =
      as.character(
        Pathway
      ),
    
    Genes =
      as.character(
        Genes
      )
  )



# ============================================================
# BASIC INPUT QC
# ============================================================

invalid_rows <- enrichment |>
  dplyr::filter(
    
    is.na(
      Pathway
    ) |
      
      is.na(
        `Enrichment FDR`
      ) |
      
      is.na(
        `Fold Enrichment`
      ) |
      
      is.na(
        nGenes
      ) |
      
      `Enrichment FDR` <= 0 |
      
      `Enrichment FDR` > 1
  )


if (
  nrow(
    invalid_rows
  ) > 0
) {
  
  warning(
    nrow(
      invalid_rows
    ),
    " ShinyGO row(s) contain missing or invalid plotting values."
  )
}



# ============================================================
# CHECK MANUSCRIPT TERMS
# ============================================================

missing_terms <- setdiff(
  selected_terms,
  enrichment$Pathway
)


if (length(missing_terms) > 0) {
  
  stop(
    "Selected manuscript term(s) absent from ShinyGO output:\n",
    paste(
      missing_terms,
      collapse = "\n"
    )
  )
}



# ============================================================
# CHECK DUPLICATE SELECTED TERMS
# ============================================================

duplicate_selected_terms <- enrichment |>
  dplyr::filter(
    Pathway %in%
      selected_terms
  ) |>
  dplyr::count(
    Pathway,
    name = "n_rows"
  ) |>
  dplyr::filter(
    n_rows > 1
  )


if (
  nrow(
    duplicate_selected_terms
  ) > 0
) {
  
  stop(
    "One or more selected pathways occur more than once ",
    "in the ShinyGO output."
  )
}



# ============================================================
# PREPARE PLOTTING TABLE
# ============================================================

plot_data <- enrichment |>
  dplyr::filter(
    Pathway %in%
      selected_terms
  ) |>
  dplyr::mutate(
    
    neg_log10_FDR =
      -log10(
        `Enrichment FDR`
      ),
    
    Pathway =
      forcats::fct_reorder(
        Pathway,
        `Fold Enrichment`
      )
  )



# ============================================================
# EXPORT SELECTED TERMS
# ============================================================

selected_output <- plot_data |>
  dplyr::transmute(
    
    pathway =
      as.character(
        Pathway
      ),
    
    fold_enrichment =
      `Fold Enrichment`,
    
    FDR =
      `Enrichment FDR`,
    
    neg_log10_FDR =
      neg_log10_FDR,
    
    candidate_genes =
      nGenes,
    
    pathway_genes =
      `Pathway Genes`,
    
    genes =
      Genes
  ) |>
  
  dplyr::arrange(
    dplyr::desc(
      fold_enrichment
    )
  )


readr::write_csv(
  selected_output,
  selected_output_file,
  na = ""
)



# ============================================================
# MAIN PLOT
# ============================================================
#
# The size guide is deliberately suppressed here.
#
# A custom bubble-size key is overlaid inside the plotting
# panel later so that the external legend does not compress
# the figure.

main_plot <- ggplot2::ggplot(
  
  plot_data,
  
  ggplot2::aes(
    
    x =
      `Fold Enrichment`,
    
    y =
      Pathway,
    
    size =
      nGenes,
    
    colour =
      neg_log10_FDR
  )
) +
  
  ggplot2::geom_point(
    alpha = 0.92
  ) +
  
  ggplot2::scale_size_continuous(
    
    range =
      c(
        4.2,
        10.5
      ),
    
    guide =
      "none"
  ) +
  
  viridis::scale_colour_viridis_c(
    
    name =
      expression(
        -log[10](FDR)
      ),
    
    option =
      "D",
    
    direction =
      1
  ) +
  
  ggplot2::scale_x_continuous(
    
    expand =
      ggplot2::expansion(
        mult =
          c(
            0.04,
            0.08
          )
      )
  ) +
  
  ggplot2::labs(
    
    x =
      "Fold enrichment",
    
    y =
      NULL
  ) +
  
  ggplot2::guides(
    
    colour =
      ggplot2::guide_colourbar(
        
        title.position =
          "right",
        
        title.hjust =
          0.5,
        
        barheight =
          grid::unit(
            68,
            "mm"
          ),
        
        barwidth =
          grid::unit(
            6,
            "mm"
          ),
        
        ticks =
          TRUE,
        
        frame.colour =
          "black"
      )
  ) +
  
  ggplot2::theme_bw(
    base_size = 12
  ) +
  
  ggplot2::theme(
    
    # --------------------------------------------------------
    # Boxed plotting panel
    # --------------------------------------------------------
    
    panel.border =
      ggplot2::element_rect(
        
        colour =
          "black",
        
        fill =
          NA,
        
        linewidth =
          0.8
      ),
    
    
    # --------------------------------------------------------
    # Vertical grid only
    # --------------------------------------------------------
    
    panel.grid.major.x =
      ggplot2::element_line(
        
        colour =
          "grey88",
        
        linewidth =
          0.4
      ),
    
    panel.grid.major.y =
      ggplot2::element_blank(),
    
    panel.grid.minor =
      ggplot2::element_blank(),
    
    
    # --------------------------------------------------------
    # Axes
    # --------------------------------------------------------
    
    axis.text.y =
      ggplot2::element_text(
        
        size =
          11,
        
        colour =
          "black",
        
        margin =
          ggplot2::margin(
            r = 8
          )
      ),
    
    axis.text.x =
      ggplot2::element_text(
        
        size =
          10.5,
        
        colour =
          "black"
      ),
    
    axis.title.x =
      ggplot2::element_text(
        
        size =
          11.5,
        
        margin =
          ggplot2::margin(
            t = 8
          )
      ),
    
    
    # --------------------------------------------------------
    # External colour legend
    # --------------------------------------------------------
    
    legend.position =
      "right",
    
    legend.title =
      ggplot2::element_text(
        size = 10.5
      ),
    
    legend.text =
      ggplot2::element_text(
        size = 9.5
      ),
    
    legend.margin =
      ggplot2::margin(
        0,
        0,
        0,
        2
      ),
    
    
    plot.margin =
      ggplot2::margin(
        8,
        4,
        8,
        8
      )
  )



# ============================================================
# CANDIDATE-GENE BUBBLE KEY
# ============================================================
#
# This is overlaid inside the main plot.
#
# It therefore does not reduce the width of the plotting panel.

size_key_data <- data.frame(
  
  n =
    c(
      5,
      15,
      30
    ),
  
  y =
    c(
      3,
      2,
      1
    ),
  
  label =
    c(
      "5 genes",
      "15 genes",
      "30 genes"
    )
)



size_key <- ggplot2::ggplot(
  
  size_key_data,
  
  ggplot2::aes(
    
    x =
      1,
    
    y =
      y,
    
    size =
      n
  )
) +
  
  ggplot2::geom_point(
    
    shape =
      21,
    
    fill =
      "grey15",
    
    colour =
      "black",
    
    stroke =
      0.4
  ) +
  
  ggplot2::geom_text(
    
    ggplot2::aes(
      
      x =
        1.38,
      
      label =
        label
    ),
    
    hjust =
      0,
    
    size =
      3.7
  ) +
  
  ggplot2::annotate(
    
    "text",
    
    x =
      0.72,
    
    y =
      3.75,
    
    label =
      "Candidate genes",
    
    hjust =
      0,
    
    size =
      4.1
  ) +
  
  ggplot2::scale_size_continuous(
    
    range =
      c(
        4.2,
        10.5
      )
  ) +
  
  ggplot2::xlim(
    0.65,
    2.25
  ) +
  
  ggplot2::ylim(
    0.45,
    4.05
  ) +
  
  ggplot2::theme_void() +
  
  ggplot2::theme(
    
    legend.position =
      "none",
    
    plot.background =
      ggplot2::element_rect(
        
        fill =
          NA,
        
        colour =
          NA
      ),
    
    panel.background =
      ggplot2::element_rect(
        
        fill =
          NA,
        
        colour =
          NA
      )
  )



# ============================================================
# ASSEMBLE FINAL FIGURE
# ============================================================

final_plot <- cowplot::ggdraw(
  main_plot
) +
  
  cowplot::draw_plot(
    
    size_key,
    
    x =
      0.665,
    
    y =
      0.085,
    
    width =
      0.205,
    
    height =
      0.255
  )



# ============================================================
# SAVE FIGURE
# ============================================================

ggplot2::ggsave(
  
  filename =
    pdf_file,
  
  plot =
    final_plot,
  
  width =
    255,
  
  height =
    160,
  
  units =
    "mm",
  
  device =
    cairo_pdf
)



ggplot2::ggsave(
  
  filename =
    png_file,
  
  plot =
    final_plot,
  
  width =
    255,
  
  height =
    160,
  
  units =
    "mm",
  
  dpi =
    600
)



# ============================================================
# QC SUMMARY
# ============================================================

qc_summary <- tibble::tibble(
  
  metric = c(
    
    "ShinyGO_rows_imported",
    
    "selected_terms_requested",
    
    "selected_terms_found",
    
    "selected_terms_missing",
    
    "duplicate_selected_terms",
    
    "invalid_ShinyGO_rows",
    
    "minimum_selected_FDR",
    
    "maximum_selected_FDR",
    
    "minimum_selected_fold_enrichment",
    
    "maximum_selected_fold_enrichment",
    
    "minimum_selected_gene_count",
    
    "maximum_selected_gene_count"
  ),
  
  
  value = c(
    
    nrow(
      enrichment
    ),
    
    length(
      selected_terms
    ),
    
    nrow(
      plot_data
    ),
    
    length(
      missing_terms
    ),
    
    nrow(
      duplicate_selected_terms
    ),
    
    nrow(
      invalid_rows
    ),
    
    min(
      plot_data$
        `Enrichment FDR`,
      na.rm = TRUE
    ),
    
    max(
      plot_data$
        `Enrichment FDR`,
      na.rm = TRUE
    ),
    
    min(
      plot_data$
        `Fold Enrichment`,
      na.rm = TRUE
    ),
    
    max(
      plot_data$
        `Fold Enrichment`,
      na.rm = TRUE
    ),
    
    min(
      plot_data$
        nGenes,
      na.rm = TRUE
    ),
    
    max(
      plot_data$
        nGenes,
      na.rm = TRUE
    )
  )
)



readr::write_csv(
  qc_summary,
  qc_output_file,
  na = ""
)



# ============================================================
# CONSOLE REPORT
# ============================================================

cat(
  "\n",
  "============================================\n",
  "Functional enrichment plot complete\n",
  "============================================\n\n",
  sep = ""
)


cat(
  "ShinyGO rows imported: ",
  nrow(
    enrichment
  ),
  "\n",
  sep = ""
)


cat(
  "Representative manuscript terms: ",
  nrow(
    plot_data
  ),
  "\n\n",
  sep = ""
)



cat(
  "Selected enrichment terms:\n"
)


print(
  
  selected_output |>
    dplyr::select(
      
      pathway,
      fold_enrichment,
      FDR,
      candidate_genes
    ),
  
  n = Inf,
  width = Inf
)



cat(
  "\nOutputs:\n",
  pdf_file,
  "\n",
  png_file,
  "\n",
  selected_output_file,
  "\n",
  qc_output_file,
  "\n",
  sep = ""
)



# ============================================================
# DISPLAY FIGURE IN RSTUDIO
# ============================================================

print(
  final_plot
)

