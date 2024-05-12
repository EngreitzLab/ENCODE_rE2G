# required packages
library(data.table)
library(dplyr)

# inputs from snakemake
biosample_config = fread(snakemake@input$biosample_config)
b = snakemake@wildcards$biosample

# merge feature tables
models_this = dplyr::filter(biosample_config, biosample==b)
models_this$feature_table = file.path(models_this$model_dir, "feature_table.tsv")
for (i in 1:nrow(models_this)){
    ft = fread(models_this$feature_table[i])
    if (i==1){df = ft} else {df = rbind(df, ft)}
}
df = dplyr::distinct(df)
fwrite(df, file = snakemake@output$biosampe_features, sep = "\t")
