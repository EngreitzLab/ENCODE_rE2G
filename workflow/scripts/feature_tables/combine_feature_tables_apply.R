# required packages
library(data.table)
library(dplyr)

# inputs from snakemake
model_dirs = snakemake@params$model_dirs %>% strsplit(" ") %>% unlist()

for (i in 1:length(model_dirs)){
	ft = fread(file.path(model_dirs[i],"feature_table.tsv"))
    if (i==1){df = ft} else {df = rbind(df, ft)}
}
df = dplyr::distinct(df)
fwrite(df, file = snakemake@output$biosample_features, sep = "\t")
