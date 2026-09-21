# ============================================================
# 00_extract_gbk_annotation.R
#
# Aim
# ---
# Extract CDS annotations and genomic coordinates from a GenBank
# (.gbk/.gbff) reference genome annotation.
#
# Two reference tables are generated:
#
#   1) reference_annotation.csv
#
#      A compact annotation table used during target-prediction
#      harmonisation.
#
#      Columns:
#
#        reference_locus
#        gene
#        annotation
#
#
#   2) reference_gene_coordinates.csv
#
#      A coordinate-aware reference table used for analyses that
#      require gene position and strand information, including
#      transcription start site (TSS) analysis.
#
#      Columns:
#
#        reference_locus
#        old_locus_tag
#        gene
#        annotation
#        strand
#        gene_start
#        gene_end
#
#
# GenBank CDS coordinates are interpreted in genomic orientation:
#
#   gene_start = lower genomic coordinate
#   gene_end   = higher genomic coordinate
#
# Strand is recorded separately as "+" or "-".
#
# For example:
#
#   CDS  123..456
#
# becomes:
#
#   strand      +
#   gene_start  123
#   gene_end    456
#
#
#   CDS  complement(123..456)
#
# becomes:
#
#   strand      -
#   gene_start  123
#   gene_end    456
#
#
# The current RefSeq locus tag is stored as reference_locus.
#
# Where available, the previous locus identifier is retained as
# old_locus_tag for traceability between annotation versions.
#
# No Bioconductor packages are required.
#
#
# INPUT
# -----
#
# data/reference/reference_genome.gbff
#
#
# OUTPUTS
# -------
#
# data/reference/reference_annotation.csv
#
# data/reference/reference_gene_coordinates.csv
#
#
# RUNNING THE SCRIPT
# ------------------
#
# 1. Place the GenBank/GBFF reference file in data/reference/.
# 2. Rename it reference_genome.gbff, or edit the USER SETTINGS.
# 3. Run this script from the repository root.
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


gbk_file <-
  "data/reference/reference_genome.gbff"


annotation_output_file <-
  "data/reference/reference_annotation.csv"


coordinate_output_file <-
  "data/reference/reference_gene_coordinates.csv"



# ============================================================
# PACKAGES
# ============================================================

required_pkgs <- c(
  "dplyr",
  "readr",
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
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})



# ============================================================
# CHECK INPUT
# ============================================================

if (!file.exists(gbk_file)) {
  
  stop(
    "GenBank file not found:\n",
    gbk_file
  )
}



# ============================================================
# READ GENBANK FILE
# ============================================================

lines <- readLines(
  gbk_file,
  warn = FALSE
)



# ============================================================
# HELPER: EXTRACT FEATURE QUALIFIER
# ============================================================

extract_qualifier <- function(
    feature_lines,
    qualifier
) {
  
  # Locate qualifiers such as:
  #
  # /locus_tag="MIDG2331_RS10200"
  # /old_locus_tag="MIDG2331_02040"
  # /gene="ilvC"
  # /product="ketol-acid reductoisomerase"
  #
  # Product descriptions can span multiple lines. Continuation
  # lines are joined until the closing quotation mark is found.
  
  
  pattern <- paste0(
    '^\\s*/',
    qualifier,
    '="'
  )
  
  
  start <- which(
    str_detect(
      feature_lines,
      pattern
    )
  )
  
  
  if (length(start) == 0) {
    return(
      NA_character_
    )
  }
  
  
  start <- start[1]
  
  
  value <- str_replace(
    feature_lines[start],
    pattern,
    ""
  )
  
  
  # Qualifier closes on the same line.
  
  if (
    str_detect(
      value,
      '"\\s*$'
    )
  ) {
    
    value <- str_replace(
      value,
      '"\\s*$',
      ""
    )
    
    
    return(
      str_squish(
        value
      )
    )
  }
  
  
  # Otherwise collect continuation lines until the closing quote.
  
  if (
    start <
    length(feature_lines)
  ) {
    
    for (
      i in seq(
        from = start + 1,
        to = length(feature_lines)
      )
    ) {
      
      continuation <- str_trim(
        feature_lines[i]
      )
      
      
      value <- paste(
        value,
        continuation
      )
      
      
      if (
        str_detect(
          continuation,
          '"\\s*$'
        )
      ) {
        break
      }
    }
  }
  
  
  value <- str_replace(
    value,
    '"\\s*$',
    ""
  )
  
  
  str_squish(
    value
  )
}



# ============================================================
# HELPER: EXTRACT CDS LOCATION
# ============================================================

extract_feature_location <- function(
    feature_line
) {
  
  # Remove the feature label itself:
  #
  #      CDS             123..456
  #
  # becomes:
  #
  #      123..456
  
  location <- str_replace(
    feature_line,
    "^\\s{5}CDS\\s+",
    ""
  )
  
  
  str_trim(
    location
  )
}



# ============================================================
# HELPER: PARSE CDS STRAND
# ============================================================

parse_strand <- function(
    location
) {
  
  if (
    is.na(location) |
    location == ""
  ) {
    
    return(
      NA_character_
    )
  }
  
  
  if (
    str_detect(
      location,
      "^complement\\("
    )
  ) {
    
    return(
      "-"
    )
  }
  
  
  "+"
}



# ============================================================
# HELPER: PARSE CDS COORDINATES
# ============================================================

parse_coordinates <- function(
    location
) {
  
  # Extract every integer present in the GenBank location.
  #
  # This supports ordinary locations such as:
  #
  #   123..456
  #
  # as well as forms such as:
  #
  #   complement(123..456)
  #   <123..456
  #   123..>456
  #   join(123..200,300..456)
  #
  # The outermost genomic coordinates are retained.
  
  numbers <- str_extract_all(
    as.character(location),
    "\\d+"
  )[[1]]
  
  
  if (
    length(numbers) == 0
  ) {
    
    return(
      tibble(
        gene_start = NA_integer_,
        gene_end = NA_integer_
      )
    )
  }
  
  
  numbers <- suppressWarnings(
    as.integer(
      numbers
    )
  )
  
  
  numbers <- numbers[
    !is.na(numbers)
  ]
  
  
  if (
    length(numbers) == 0
  ) {
    
    return(
      tibble(
        gene_start = NA_integer_,
        gene_end = NA_integer_
      )
    )
  }
  
  
  tibble(
    
    gene_start =
      min(numbers),
    
    gene_end =
      max(numbers)
  )
}



# ============================================================
# IDENTIFY CDS FEATURES
# ============================================================

# GenBank feature lines normally begin with five spaces followed
# by the feature type:
#
#      CDS             123..456

cds_starts <- which(
  str_detect(
    lines,
    "^\\s{5}CDS\\s+"
  )
)


if (
  length(cds_starts) == 0
) {
  
  stop(
    "No CDS features were found in the GenBank file."
  )
}



# ============================================================
# EXTRACT CDS RECORDS
# ============================================================

records <- vector(
  "list",
  length(
    cds_starts
  )
)



for (
  i in seq_along(
    cds_starts
  )
) {
  
  feature_start <- cds_starts[i]
  
  
  # ----------------------------------------------------------
  # Find the next GenBank feature
  # ----------------------------------------------------------
  
  subsequent_features <- which(
    
    seq_along(lines) >
      feature_start &
      
      str_detect(
        lines,
        "^\\s{5}[A-Za-z][A-Za-z0-9_'-]*\\s+"
      )
  )
  
  
  if (
    length(
      subsequent_features
    ) > 0
  ) {
    
    feature_end <-
      subsequent_features[1] - 1
    
  } else {
    
    feature_end <-
      length(lines)
  }
  
  
  feature_lines <- lines[
    feature_start:
      feature_end
  ]
  
  
  # ----------------------------------------------------------
  # Extract CDS location
  # ----------------------------------------------------------
  
  location <- extract_feature_location(
    lines[
      feature_start
    ]
  )
  
  
  coordinates <- parse_coordinates(
    location
  )
  
  
  strand <- parse_strand(
    location
  )
  
  
  # ----------------------------------------------------------
  # Extract qualifiers
  # ----------------------------------------------------------
  
  records[[i]] <- tibble(
    
    reference_locus =
      extract_qualifier(
        feature_lines,
        "locus_tag"
      ),
    
    old_locus_tag =
      extract_qualifier(
        feature_lines,
        "old_locus_tag"
      ),
    
    gene =
      extract_qualifier(
        feature_lines,
        "gene"
      ),
    
    annotation =
      extract_qualifier(
        feature_lines,
        "product"
      ),
    
    strand =
      strand,
    
    gene_start =
      coordinates$gene_start,
    
    gene_end =
      coordinates$gene_end
  )
}



# ============================================================
# COMBINE AND CLEAN
# ============================================================

reference_coordinates <- bind_rows(
  records
) %>%
  
  filter(
    !is.na(
      reference_locus
    ),
    reference_locus != ""
  ) %>%
  
  mutate(
    
    reference_locus =
      str_to_upper(
        str_trim(
          reference_locus
        )
      ),
    
    old_locus_tag =
      str_trim(
        old_locus_tag
      ),
    
    gene =
      str_trim(
        gene
      ),
    
    annotation =
      str_squish(
        annotation
      ),
    
    old_locus_tag =
      na_if(
        old_locus_tag,
        ""
      ),
    
    gene =
      na_if(
        gene,
        ""
      ),
    
    annotation =
      na_if(
        annotation,
        ""
      )
  ) %>%
  
  distinct(
    reference_locus,
    .keep_all = TRUE
  ) %>%
  
  arrange(
    reference_locus
  )



# ============================================================
# CREATE COMPACT ANNOTATION TABLE
# ============================================================

reference_annotation <- reference_coordinates %>%
  select(
    reference_locus,
    gene,
    annotation
  )



# ============================================================
# CREATE COORDINATE TABLE
# ============================================================

reference_gene_coordinates <- reference_coordinates %>%
  select(
    reference_locus,
    old_locus_tag,
    gene,
    annotation,
    strand,
    gene_start,
    gene_end
  )



# ============================================================
# BASIC QC
# ============================================================

missing_coordinates <- reference_gene_coordinates %>%
  filter(
    is.na(strand) |
      is.na(gene_start) |
      is.na(gene_end)
  )


unexpected_strand <- reference_gene_coordinates %>%
  filter(
    !is.na(strand),
    !strand %in% c(
      "+",
      "-"
    )
  )


missing_old_locus <- reference_gene_coordinates %>%
  filter(
    is.na(old_locus_tag)
  )



if (
  nrow(
    unexpected_strand
  ) > 0
) {
  
  warning(
    nrow(
      unexpected_strand
    ),
    " CDS record(s) contain unexpected strand values."
  )
}



# ============================================================
# EXPORT
# ============================================================

dir.create(
  dirname(
    annotation_output_file
  ),
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  dirname(
    coordinate_output_file
  ),
  recursive = TRUE,
  showWarnings = FALSE
)



write_csv(
  reference_annotation,
  annotation_output_file,
  na = ""
)


write_csv(
  reference_gene_coordinates,
  coordinate_output_file,
  na = ""
)



# ============================================================
# RUN SUMMARY
# ============================================================

cat(
  "\n============================================\n",
  "GenBank reference extraction complete\n",
  "============================================\n\n",
  sep = ""
)



cat(
  "CDS features detected: ",
  length(
    cds_starts
  ),
  "\n",
  sep = ""
)



cat(
  "Reference loci exported: ",
  nrow(
    reference_gene_coordinates
  ),
  "\n",
  sep = ""
)



cat(
  "Records without complete coordinates: ",
  nrow(
    missing_coordinates
  ),
  "\n",
  sep = ""
)



cat(
  "Records without old locus tag: ",
  nrow(
    missing_old_locus
  ),
  "\n\n",
  sep = ""
)



cat(
  "Annotation output:\n",
  annotation_output_file,
  "\n\n",
  sep = ""
)



cat(
  "Coordinate output:\n",
  coordinate_output_file,
  "\n",
  sep = ""
)

