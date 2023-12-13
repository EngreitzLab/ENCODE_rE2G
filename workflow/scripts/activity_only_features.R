## Create Activity-only feature table for an ABC sample

# required packages and functions
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  source(file.path(snakemake@scriptdir, "get_fill_values.R"))
})

# load feature config file
config <- fread(snakemake@input$feature_table_file)
# load ABC table
abc <- fread(snakemake@input$abc)

# load individual feature tables
NumCandidateEnhGene <- fread(snakemake@input$NumCandidateEnhGene,
  header = TRUE,
  col.names = c("name", "TargetGene", "numCandidateEnhGene")
)
NumTSSEnhGene <- fread(snakemake@input$NumTSSEnhGene,
  header = TRUE,
  col.names = c("name", "TargetGene", "numTSSEnhGene")
)
NumEnhancersEG5kb <- fread(snakemake@input$NumEnhancersEG5kb,
  header = FALSE, sep = "\t",
  col.names = c("name", "numNearbyEnhancers")
)
SumEnhancersEG5kb <- fread(snakemake@input$SumEnhancersEG5kb,
  header = FALSE, sep = "\t",
  col.names = c("name", "sumNearbyEnhancers")
)
ubiqExprGenes <- fread(snakemake@input$ubiqExprGenes)

# process abc input and compute squared contact and activity measurements
to_square <- dplyr::filter(config, square_source_col != "NA")
for (i in seq_len(nrow(to_square))) {
  abc[[to_square$input_col[i]]] <- abc[[to_square$square_source_col[i]]]^2
}

# extract Activity-only feature columns and output names as vector
activity_feat_cols <- config %>%
  select(feature, input_col) %>%
  deframe()

# core columns from ABC for output
core_cols <- c(
  "chr", "start", "end", "name", "class", "TargetGene", "TargetGeneTSS",
  "TargetGeneIsExpressed", "TargetGeneEnsembl_ID", "isSelfPromoter",
  "CellType", "distance"
)
# filter to existing core columns for backwards compatability
core_cols <- core_cols[core_cols %in% colnames(abc)]

# select all core columns from ABC table for output
output <- select(abc, all_of(core_cols))

# assemble feature table ---------------------------------------------------------------------------

# add features from ABC table
output <- abc %>%
  select(name, TargetGene, any_of(activity_feat_cols)) %>%
  left_join(output, ., by = c("name", "TargetGene"))

# add number of candidate enhancers between enhancers and genes
output <- NumCandidateEnhGene %>%
  select(name, TargetGene, any_of(activity_feat_cols)) %>%
  left_join(output, ., by = c("name", "TargetGene"))

# add number of protein-coding TSSs between enhancer and target
output <- NumTSSEnhGene %>%
  select(name, TargetGene, any_of(activity_feat_cols)) %>%
  left_join(output, ., by = c("name", "TargetGene"))

# add number of nearby enhancers
output <- NumEnhancersEG5kb %>%
  select(name, any_of(activity_feat_cols)) %>%
  left_join(output, ., by = "name")

# add sum of activities of nearby enhancers
output <- SumEnhancersEG5kb %>%
  select(name, any_of(activity_feat_cols)) %>%
  left_join(output, ., by = "name")

# add ubiquitously expressed gene info
output <- ubiqExprGenes %>%
  select(TargetGene = GeneSymbol, any_of(activity_feat_cols)) %>%
  left_join(output, ., by = "TargetGene")

# fill in NAs and write to output

# get fill values for each feature
config <- filter(config, feature %in% colnames(output))
fill_values <- get_fill_values(output, config = config)

# replace NAs with fill values
output <- replace_na(output, replace = fill_values)

# reorder columns for output
output <- select(output, all_of(core_cols), all_of(names(activity_feat_cols)))

# save output to file
fwrite(output, file = snakemake@output[[1]], sep = "\t", na = "NA", quote = FALSE)
