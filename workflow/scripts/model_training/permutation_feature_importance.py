import pickle
import click
import numpy as np
import pandas as pd
import scipy
from sklearn.metrics import precision_recall_curve, auc, log_loss, roc_auc_score
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LogisticRegression
from training_functions import (
    statistic_aupr,
    statistic_precision,
    train_and_predict_once,
    bootstrap_pvalue,
    statistic_delta_aupr,
    statistic_delta_precision_at_threshold,
    threshold_70_pct_recall,
)


def permutation_feature_importance(
    df_dataset,
    feature_table,
    model_name,
    epsilon,
    params,
    n_repeats=20,
    polynomial=False,
):
    feature_list_core = feature_table["feature"]
    model_name_core = model_name

    # add polynomial features
    if polynomial:
        X_core = df_dataset.loc[:, feature_list_core]
        poly = PolynomialFeatures(degree=2)
        X_2 = poly.fit_transform(X_core)
        polynomial_names = poly.get_feature_names_out(X_core.columns)
        X_2 = pd.DataFrame(X_2, columns=polynomial_names)
        X = X_2
        feature_list = pd.Series(polynomial_names)
    else:
        X = df_dataset.loc[:, feature_list_core]
        feature_list = feature_list_core

    # # tranforms features
    X = np.log(np.abs(X) + epsilon)
    Y_true = df_dataset["Regulated"].values.astype(np.int64)

    # train full model
    model_name = model_name_core + "_full"
    df_dataset = train_and_predict_once(
        df_dataset, X, Y_true, feature_list, model_name, params
    )
    Y_full = df_dataset[model_name + ".Score"]

    df = pd.DataFrame(columns=["feature_permuted", "delta_aupr", "delta_precision"])

    for i in range(len(feature_list)):  # iterate through features
        print(feature_list[i])
        delta_aupr_feature = []
        delta_precision_feature = []
        original_values = X[feature_list[i]]

        for j in range(n_repeats):
            X[feature_list[i]] = np.random.permutation(X[feature_list[i]])
            model_name = model_name_core + "_shuff_" + feature_list[i]
            df_dataset = train_and_predict_once(
                df_dataset, X, Y_true, feature_list, model_name, params
            )
            Y_shuffle = df_dataset[model_name + ".Score"]

            delta_aupr = statistic_delta_aupr(Y_true, Y_full, Y_shuffle)
            thresh_full = threshold_70_pct_recall(Y_true, Y_full)
            thresh_shuffle = threshold_70_pct_recall(Y_true, Y_shuffle)
            delta_precision = statistic_delta_precision_at_threshold(
                Y_true, Y_full, Y_shuffle, thresh_full, thresh_shuffle
            )

            delta_aupr_feature = delta_aupr_feature + [delta_aupr]
            delta_precision_feature = delta_precision_feature + [delta_precision]

        df_temp = pd.DataFrame(
            {
                "delta_aupr": delta_aupr_feature,
                "delta_precision": delta_precision_feature,
            }
        )
        df_temp["feature_permuted"] = feature_list[i]
        df = pd.concat([df, df_temp])

        X[feature_list[i]] = original_values  # reset feature values

    return df


@click.command()
@click.option("--crispr_features_file", required=True)
@click.option("--feature_table_file", required=True)
@click.option("--out_dir", required=True)
@click.option("--polynomial", type=bool, default=False)
@click.option("--epsilon", type=float, default=0.01)
@click.option("--n_repeats", type=float, default=20)
@click.option("--params_file", required=True)
def main(
    crispr_features_file,
    feature_table_file,
    out_dir,
    polynomial,
    epsilon,
    n_repeats,
    params_file,
):
    model_name = "ENCODE-rE2G"
    n_repeats = int(n_repeats)
    df_dataset = pd.read_csv(crispr_features_file, sep="\t")
    feature_table = pd.read_csv(feature_table_file, sep="\t")

    with open(params_file, "rb") as handle:
        params = pickle.load(handle)

    res = permutation_feature_importance(
        df_dataset, feature_table, model_name, epsilon, params, n_repeats, polynomial
    )

    res.to_csv(out_dir + "/permutation_feature_importance.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
