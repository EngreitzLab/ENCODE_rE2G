# required packages
library(data.table)
library(dplyr)

# inputs from snakemake
model_config = fread(snakemake@input$model_config)
ds = snakemake@wildcards$dataset

# merge feature tables
models_this = dplyr::filter(model_config, dataset==ds)
for (i in 1:nrow(models_this)){
    ft = fread(models_this$feature_table[i])
    if (i==1){df = ft} else {df = rbind(df, ft)}
}

# for sc-E2G pipeline
if (("ARC.E2G.Score" %in% df$feature) | ("Kendall" %in% df$feature)){
	gene_expr_row = data.frame("meanLogNormRNA", "mean_log_normalized_rna", NA, "mean", 0, "Mean log normalized RNA expression"); names(gene_expr_row) = colnames(df)
	Kendall_row = data.frame("Kendall", "Kendall", NA, "max", 0, "Kendall correlation"); names(Kendall_row) = colnames(df)
	ARC_row = data.frame("ARC.E2G.Score", "ARC.E2G.Score", NA, "mean", 0, "ARC-E2G score"); names(ARC_row) = colnames(df)
	df = bind_rows(df, gene_expr_row, Kendall_row, ARC_row)
}

df = dplyr::distinct(df)
fwrite(df, file = snakemake@output$dataset_features, sep = "\t")
