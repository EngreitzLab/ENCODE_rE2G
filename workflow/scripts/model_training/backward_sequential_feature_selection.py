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
    threshold_70_pct_recall,
    statistic_precision_at_threshold,
    statistic_delta_precision_at_threshold,
)


def SBFS(df_dataset, feature_table, model_name, epsilon, params, polynomial=False):
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

    n_features = len(feature_list)

    feature_removed = []
    all_auprs = []
    all_precisions = []

    # train all feature model
    model_name = model_name_core + "_full"
    df_dataset = train_and_predict_once(
        df_dataset, X, Y_true, feature_list, model_name, params
    )
    Y_pred = df_dataset[model_name + ".Score"]
    aupr = statistic_aupr(Y_true, Y_pred)
    thresh = threshold_70_pct_recall(Y_true, Y_pred)
    precision_at_70_pct_recall = statistic_precision_at_threshold(
        Y_true, Y_pred, thresh
    )

    feature_removed = feature_removed + ["None"]
    all_auprs = all_auprs + [aupr]
    all_precisions = all_precisions + [precision_at_70_pct_recall]

    for i in range(n_features - 1):  # iterate (# of feature) times
        best_k = -1
        best_aupr = 0
        best_precision = 0
        for k in range(
            len(feature_list)
        ):  # iterate through each feature, remove one with least effect
            new_feature = feature_list[k]  # feature to remove
            print(new_feature)
            features = feature_list.loc[feature_list != new_feature].reset_index(
                drop=True
            )
            print(features)
            model_name = model_name_core + "_" + str(i + 1) + "_" + str(k + 1)
            print(model_name)

            df_dataset = train_and_predict_once(
                df_dataset, X, Y_true, features, model_name, params
            )

            # evaluate model once trained
            Y_pred = df_dataset[model_name + ".Score"]
            aupr = statistic_aupr(Y_true, Y_pred)
            thresh = threshold_70_pct_recall(Y_true, Y_pred)
            precision_at_70_pct_recall = statistic_precision_at_threshold(
                Y_true, Y_pred, thresh
            )

            if aupr >= best_aupr:
                best_aupr = aupr
                best_k = k
                best_precision = precision_at_70_pct_recall

        # update best feature
        feature_removed = feature_removed + [feature_list[best_k]]
        all_auprs = all_auprs + [best_aupr]
        all_precisions = all_precisions + [best_precision]
        feature_list = feature_list.loc[
            feature_list != feature_list[best_k]
        ].reset_index(drop=True)
    # add remaining feature
    feature_removed = feature_removed + [feature_list[0]]
    print(feature_removed)
    print(all_auprs)
    print(all_precisions)

    return feature_removed


def SBFS_significance(
    df_dataset,
    feature_table,
    model_name,
    epsilon,
    params,
    polynomial=False,
    n_boot=1000,
):
    feature_list = SBFS(
        df_dataset, feature_table, model_name, epsilon, params, polynomial
    )  # get order of features
    feature_list_core = feature_table["feature"]
    model_name_core = model_name

    # calculate polynomial features
    if polynomial:
        X_core = df_dataset.loc[:, feature_list_core]
        poly = PolynomialFeatures(degree=2)
        X_2 = poly.fit_transform(X_core)
        polynomial_names = poly.get_feature_names_out(X_core.columns)
        X_2 = pd.DataFrame(X_2, columns=polynomial_names)
        X = X_2
    else:
        X = df_dataset.loc[:, feature_list_core]

    X = np.log(np.abs(X) + epsilon)
    Y_true = df_dataset["Regulated"].values.astype(np.int64)

    df = pd.DataFrame(
        columns=[
            "feature_removed",
            "aupr",
            "delta_aupr",
            "delta_aupr_low",
            "delta_aupr_high",
            "pval_aupr",
            "precision",
            "delta_precision",
            "delta_precision_low",
            "delta_precision_high",
            "pval_precision",
        ]
    )
    # train first model (all features)
    model_name = model_name_core + "_full"
    feature_list.remove("None")
    df_dataset = train_and_predict_once(
        df_dataset, X, Y_true, feature_list, model_name, params
    )
    Y_new = df_dataset[model_name + ".Score"]

    data = (Y_true, Y_new)
    res_aupr = scipy.stats.bootstrap(
        data,
        statistic_aupr,
        n_resamples=n_boot,
        paired=True,
        confidence_level=0.95,
        method="BCa",
    )
    thresh = threshold_70_pct_recall(Y_true, Y_new)
    res_precision = scipy.stats.bootstrap(
        data,
        lambda Y_true, Y_new: statistic_precision_at_threshold(Y_true, Y_new, thresh),
        n_resamples=n_boot,
        paired=True,
        confidence_level=0.95,
        method="BCa",
    )
    aupr = np.mean(res_aupr.bootstrap_distribution)
    precision = np.mean(res_precision.bootstrap_distribution)

    df_temp = pd.DataFrame(
        {
            "feature_removed": "None",
            "aupr": aupr,
            "delta_aupr": 0,
            "delta_aupr_low": 0,
            "delta_aupr_high": 0,
            "pval_aupr": 1,
            "precision": precision,
            "delta_precision": 0,
            "delta_precision_low": 0,
            "delta_precision_high": 0,
            "pval_precision": 1,
        },
        index=[0],
    )
    df = pd.concat([df, df_temp])

    Y_last = Y_new

    # loop feature list from SBFS through to train next model and compare
    for i in range(len(feature_list)):
        to_remove = feature_list[0]
        print(to_remove)
        if len(feature_list) == 1:
            df_dataset["all_true"] = 1
            Y_new = df_dataset["all_true"]
        else:
            feature_list.remove(to_remove)
            model_name = model_name_core + "_" + str(i + 1)
            df_dataset = train_and_predict_once(
                df_dataset, X, Y_true, feature_list, model_name, params
            )
            Y_new = df_dataset[model_name + ".Score"]

        data_compare = (Y_true, Y_last, Y_new)
        thresh_last = threshold_70_pct_recall(Y_true, Y_last)
        thresh_new = threshold_70_pct_recall(Y_true, Y_new)
        res_delta_aupr = scipy.stats.bootstrap(
            data_compare,
            statistic_delta_aupr,
            n_resamples=n_boot,
            paired=True,
            confidence_level=0.95,
            method="BCa",
        )
        res_delta_precision = scipy.stats.bootstrap(
            data_compare,
            lambda Y_true, Y_last, Y_new,: statistic_delta_precision_at_threshold(
                Y_true, Y_last, Y_new, thresh_last, thresh_new
            ),
            n_resamples=n_boot,
            paired=True,
            confidence_level=0.95,
            method="BCa",
        )
        delta_aupr = np.mean(
            res_delta_aupr.bootstrap_distribution
        )  # aupr_new - aupr_last
        delta_precision = np.mean(res_delta_precision.bootstrap_distribution)
        pval_aupr = bootstrap_pvalue(delta_aupr, res_delta_aupr)
        pval_precision = bootstrap_pvalue(delta_precision, res_delta_precision)

        data = (Y_true, Y_new)
        res_aupr = scipy.stats.bootstrap(
            data,
            statistic_aupr,
            n_resamples=n_boot,
            paired=True,
            confidence_level=0.95,
            method="BCa",
        )
        thresh = threshold_70_pct_recall(Y_true, Y_new)
        res_precision = scipy.stats.bootstrap(
            data,
            lambda Y_true, Y_new: statistic_precision_at_threshold(
                Y_true, Y_new, thresh
            ),
            n_resamples=n_boot,
            paired=True,
            confidence_level=0.95,
            method="BCa",
        )
        aupr = np.mean(res_aupr.bootstrap_distribution)
        precision = np.mean(res_precision.bootstrap_distribution)

        df_temp = pd.DataFrame(
            {
                "feature_removed": to_remove,
                "aupr": aupr,
                "delta_aupr": delta_aupr,
                "delta_aupr_low": res_delta_aupr.confidence_interval[0],
                "delta_aupr_high": res_delta_aupr.confidence_interval[1],
                "pval_aupr": pval_aupr,
                "precision": precision,
                "delta_precision": delta_precision,
                "delta_precision_low": res_delta_precision.confidence_interval[0],
                "delta_precision_high": res_delta_precision.confidence_interval[1],
                "pval_precision": pval_precision,
            },
            index=[0],
        )
        df = pd.concat([df, df_temp])
        Y_last = Y_new

    return df


@click.command()
@click.option("--crispr_features_file", required=True)
@click.option("--feature_table_file", required=True)
@click.option("--out_dir", required=True)
@click.option("--polynomial", type=bool, default=False)
@click.option("--epsilon", type=float, default=0.01)
@click.option("--params_file", required=True)
def main(
    crispr_features_file, feature_table_file, out_dir, polynomial, params_file, epsilon
):
    model_name = "ENCODE-rE2G"
    df_dataset = pd.read_csv(crispr_features_file, sep="\t")
    feature_table = pd.read_csv(feature_table_file, sep="\t")
    with open(params_file, "rb") as handle:
        params = pickle.load(handle)

    res = SBFS_significance(
        df_dataset, feature_table, model_name, epsilon, params, polynomial, n_boot=1000
    )

    res.to_csv(out_dir + "/backward_feature_selection.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
