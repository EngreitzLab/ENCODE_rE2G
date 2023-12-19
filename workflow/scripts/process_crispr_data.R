## Create Activity-only feature table for an ABC sample

# required packages and functions
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

# load inputs
df <- fread(snakemake@input$crispr_features)
filter_features <- (snakemake@params$filter_features) %>% strsplit(" ") %>% unlist()
genes <- fread(snakemake@params$genes)
colnames(genes) <- c('chr', 'start', 'end', 'gene', 'score', 'strand')

# rename some column names for later utility
df = dplyr::rename(df, 'chr'='chrom', 'start'='chromStart', 'end'='chromEnd', 'TargetGene'='measuredGeneSymbol')

# drop elements where Regulated=NA
df = dplyr::filter(df, Regulated!="NA")

# # if filter features are present, keep entries where they are both nonzero
# ff_present = filter_features[filter_features %in% colnames(df)]
# df= dplyr::filter(df, if_all(ff_present, ~ . !=0))

# filter to elements with target gene in gene universe
df = dplyr::filter(df, TargetGene %in% genes$gene)

# save output to file
fwrite(df, file = snakemake@output$processed, sep = "\t", na = "NA", quote = FALSE)
