## Create Activity-only feature table for an ABC sample

# required packages and functions
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

combine_crispr_features <- function(space_sep_string) {
  crispr_feature_files <- space_sep_list %>% strsplit(" ") %>% trimws()
  crispr_features <- lapply(crispr_feature_files, fread)

  # get common columns (e.g. features in all)
  col_names_list <- lapply(crispr_features, colnames)
  common_cols <- Reduce(intersect, col_names_list)

  # filter crispr features to selected columns then combine
  crispr_features <- lapply(crispr_features, function(df) select(df, all_of(common_cols))) %>% 
    rbindlist() %>% as.data.frame()

  return(crispr_features)
}



# load inputs
crispr_features <- combine_crispr_features(snakemake@input$crispr_features)
crispr_missing <- combine_crispr_features(snakemake@input$crispr_missing)


genes <- fread(snakemake@params$genes)
colnames(genes) <- c('chr', 'start', 'end', 'gene', 'score', 'strand', 'Ensembl_ID', 'gene_type')

# rename some column names for later utility
crispr_features <- crispr_features %>% 
  rename('chr'='chrom', 'start'='chromStart', 'end'='chromEnd', 'TargetGene'='measuredGeneSymbol')
crispr_missing <- crispr_missing %>% 
  rename('chr'='chrom', 'start'='chromStart', 'end'='chromEnd', 'TargetGene'='measuredGeneSymbol')

# drop elements where Regulated=NA
crispr_features = filter(crispr_features, Regulated != "NA" & !is.na(Regulated))
crispr_missing = filter(crispr_missing, Regulated != "NA" & !is.na(Regulated))

# filter to elements with target gene in gene universe (should be redundant with applying filter_genes=tss_universe in feature overlap)
crispr_features = filter(crispr_features, TargetGene %in% genes$gene)
crispr_missing = filter(crispr_missing, TargetGene %in% genes$gene)

# save output to file
fwrite(crispr_features, file = snakemake@output$processed, sep = "\t", na = "NA", quote = FALSE)
fwrite(crispr_missing, file = snakemake@output$missing, sep = "\t", na = "NA", quote = FALSE)
