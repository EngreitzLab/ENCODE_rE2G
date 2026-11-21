import pickle

import click
import numpy as np
import pandas as pd

MODEL = "ENCODE-rE2G"


def make_e2g_predictions(df_enhancers, feature_list, model_dir, epsilon):
    # transform the features
    X = df_enhancers.loc[:, feature_list]
    X = np.log(np.abs(X) + epsilon)
    chr_list = np.unique(df_enhancers["chr"])

    for chr in chr_list:
        idx_test = df_enhancers[df_enhancers["chr"] == chr].index.values

        if len(idx_test) > 0:
            X_test = X.loc[idx_test, :]

            with open(f"{model_dir}/model_test_{chr}.pkl", "rb") as f:
                model = pickle.load(f)

            probs = model.predict_proba(X_test)
            df_enhancers.loc[idx_test, MODEL + ".Score"] = probs[:, 1]

    return df_enhancers


def main():
    gw_features_file = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/results/2024_0329_final_models_two_inputs/DNase_both_megamap/genomewide_features.tsv.gz"  # two bam files (gw features)
    feature_table_file = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/resources/feature_tables/final_feature_set_DNase_hic.tsv"  # dnase/hic final 8
    model_dir = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/results/2024_0329_final_models_two_inputs/DNase_ENCFF860XAE_megamap/DNase_ENCFF860XAE_megamap_8features/model"  # crossval model trained on one bam file, 'model_test_chrN.pkl'
    output_file = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/results/2024_0329_final_models_two_inputs/encode_re2g_K562.genome_wide_predictions.crossvalidation_models.tsv.gz"

    # output_file = '/oak/stanford/groups/engreitz/Projects/Benchmarking/Revisions/Predictors/ENCODE-rE2G/dhs_only/encode_re2g_K562.genome_wide_predictions.crossvalidation_models.tsv.gz'

    feature_table = pd.read_csv(feature_table_file, sep="\t")
    feature_list = feature_table["feature"]

    gw_features = pd.read_csv(gw_features_file, sep="\t")
    gw_features = gw_features.replace([np.inf, -np.inf], np.nan)
    gw_features = gw_features.fillna(0)

    predictions = make_e2g_predictions(
        gw_features, feature_list, model_dir, epsilon=0.01
    )
    predictions.to_csv(output_file, compression="gzip", sep="\t", index=False)


if __name__ == "__main__":
    main()
