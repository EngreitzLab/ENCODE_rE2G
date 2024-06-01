# required packages
library(data.table)
library(dplyr)

# inputs from snakemake
model_dirs = snakemake@params$model_dirs %>% strsplit(" ") %>% unlist()

for (i in 1:length(model_dirs)){
	ft = fread(file.path(model_dirs[i],"feature_table.tsv"))
    if (i==1){df = ft} else {df = rbind(df, ft)}
}

# for sc-E2G pipeline
if ("ARC.E2G.Score" %in% df$feature){
	Kendall_row = data.frame("Kendall", "Kendall", NA, "max", "overlap", "Kendall correlation"); names(Kendall_row) = colnames(df)
	df = rbind(df, Kendall_row)
}

df = dplyr::distinct(df)
fwrite(df, file = snakemake@output$biosample_features, sep = "\t")
