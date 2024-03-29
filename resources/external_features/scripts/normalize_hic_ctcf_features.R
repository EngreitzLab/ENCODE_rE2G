## Normalize CTCF PET and Hi-C loop count features

# required packages
library(data.table)
library(dplyr)

compute_norm_features <- function(input_file, output_file) {
    df = fread(input_file, header=TRUE)

    df <- df %>%
        mutate(pet_outside_norm = (pet_outside + 1) / (prom_PET + 1),
        pet_cross_norm   = (pet_cross + 1) / (prom_PET + 1),
        hicloop_outside_norm = (hicloop_outside + 1) / (prom_hicloop + 1),
        hicloop_cross_norm   = (hicloop_cross + 1) / (prom_hicloop + 1))
    
    fwrite(df, output_file, sep="\t")

}

input_file = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/resources/external_features/K562/EnhancerPredictionsAllPutative_CTCF_HiC_CIA_hg38.txt.gz"
output_file = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/resources/external_features/K562/EnhancerPredictionsAllPutative_CTCF_HiC_CIA_hg38.with_norm.txt.gz"
compute_norm_features(input_file, output_file)
