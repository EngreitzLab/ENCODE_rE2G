import pickle
import click
import numpy as np
import pandas as pd
import scipy
from training_functions import (
    statistic_aupr,
    statistic_precision_at_threshold,
    train_and_predict_once,
    threshold_70_pct_recall
)


def compare_feature_sets(df_dataset, feature_table, epsilon, params, n_boot):
    feature_list = feature_table["feature"]
    X = df_dataset.loc[:, feature_list]
    X = np.log(np.abs(X) + epsilon)
    Y_true = df_dataset["Regulated"].values.astype(np.int64)

    # specify all feature sets
    n_features = len(feature_list)
    n_sets_iter = 2**n_features
    df = pd.DataFrame(columns=feature_list)

    for i in range(1, n_sets_iter):
        binary_n = str(bin(i).replace("0b", ""))  # convert to binary
        fill_zeros = list(
            map(int, [*binary_n.zfill(n_features)])
        )  # convert to list of ints of length n_features
        df.loc[len(df)] = fill_zeros

    df["features"] = df.apply(lambda row: list(df.columns[row == 1]), axis=1)
    df["n_features"] = 0.0
    df["AUPRC"] = 0.0
    df["AUPRC_95CI_low"] = 0.0
    df["AUPRC_95CI_high"] = 0.0
    df["precision_70pct_recall"] = 0.0
    df["precision_95CI_low"] = 0.0
    df["precision_95CI_high"] = 0.0

    for i in range(len(df)):  # iterate through sets
        model_name = "row_" + str(i)
        if i % 100 == 0:
            print(model_name)
        features = df.loc[i, "features"]
        df.loc[i, "n_features"] = len(features)
        df_dataset = train_and_predict_once(
            df_dataset, X, Y_true, features, model_name, params
        )
        Y_pred = df_dataset[model_name + ".Score"]

        res_aupr = scipy.stats.bootstrap(
            (Y_true, Y_pred),
            statistic_aupr,
            n_resamples=n_boot,
            paired=True,
            confidence_level=0.95,
            method="BCa",
        )
        df.loc[i, "AUPRC"] = np.mean(res_aupr.bootstrap_distribution)
        df.loc[i, "AUPRC_95CI_low"] = res_aupr.confidence_interval[0]
        df.loc[i, "AUPRC_95CI_high"] = res_aupr.confidence_interval[1]

        thresh = threshold_70_pct_recall(Y_true, Y_pred)

        res_prec = scipy.stats.bootstrap(
            (Y_true, Y_pred),
            lambda Y_true, Y_pred: statistic_precision_at_threshold(Y_true, Y_pred, thresh),
            n_resamples=n_boot,
            paired=True,
            confidence_level=0.95,
            method="BCa",
        )
        df.loc[i, "precision_70pct_recall"] = np.mean(res_prec.bootstrap_distribution)
        df.loc[i, "precision_95CI_low"] = res_prec.confidence_interval[0]
        df.loc[i, "precision_95CI_high"] = res_prec.confidence_interval[1]

    # sort table by AUPRC
    df = df.sort_values(by="AUPRC", ascending=False)

    return df


@click.command()
@click.option("--crispr_features_file", required=True)
@click.option("--feature_table_file", required=True)
@click.option("--out_dir", required=True)
@click.option("--epsilon", type=float, default=0.01)
@click.option("--params_file", required=True)
def main(crispr_features_file, feature_table_file, out_dir, epsilon, params_file):
    df_dataset = pd.read_csv(crispr_features_file, sep="\t")
    feature_table = pd.read_csv(feature_table_file, sep="\t")

    with open(params_file, "rb") as handle:
        params = pickle.load(handle)

    res = compare_feature_sets(df_dataset, feature_table, epsilon, params, n_boot=1000)

    res.to_csv(out_dir + "/all_feature_sets.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
