# required packages
library(data.table)
library(dplyr)

# inputs from snakemake
model_config = fread(snakemake@input$model_config)
ds = snakemake@wildcards$dataset

# merge feature tables for models with this dataset
ft_files = c()
for (i in 1:nrow(model_config)){
	model_datasets = model_config$dataset[i] %>% strsplit(",") %>% trimws()
	if (ds %in% model_datasets) {
		ft_files = c(ft_files, model_config$feature_table[i])
	}
}

df <- lapply(ft_files, fread) %>% rbindlist() %>% as.data.frame()

# for sc-E2G pipeline
if (("ARC.E2G.Score" %in% df$feature) | ("Kendall" %in% df$feature)){
	ARC_rows = data.frame(c("RNA_meanLogNorm", "RNA_pseudobulkTPM", "RNA_percentCellsDetected", "Kendall", "ARC.E2G.Score", "ABC.Score", "normalizedATAC_enh"),
		c("mean_log_normalized_rna", "RnaPseudobulkTPM", "RnaDetectedPercent", "Kendall", "ARC.E2G.Score", "ABC.Score", "normalized_atac_enh"),
		c(NA, NA, NA, NA, NA, NA, NA),
		c("mean", "mean", "mean", "max", "mean", "sum", "mean"),
		c(0, 0, 0, 0, 0, 0, 0),
		c("Mean log normalized RNA expression", "RNA pseudobulk TPM", "RNA percent cells detected", "Kendall correlation", "ARC-E2G score", "ABC score", "ATAC signal at E"));

	ARC_rows = ARC_rows %>% setNames(colnames(df)) %>%
		dplyr::filter(!(feature %in% df$feature)) # get rid of duplicate rows (e.g. ABC.Score if already there)
	
	df = rbind(df, ARC_rows)
}

df = dplyr::distinct(df)
fwrite(df, file = snakemake@output$dataset_features, sep = "\t")
