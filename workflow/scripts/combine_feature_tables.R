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
df = dplyr::distinct(df)
fwrite(df, file = snakemake@output$dataset_features, sep = "\t")
