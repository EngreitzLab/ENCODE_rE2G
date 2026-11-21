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


def calculate_ablation(
    df_dataset, feature_table, missing_df, params, out_dir, epsilon, polynomial=False
):
    feature_list_core = feature_table["feature"]
    X = df_dataset.loc[:, feature_list_core]

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
    model_name_full = "rE2G_full"
    df_dataset = train_and_predict_once(
        df_dataset, X, Y_true, feature_list, model_name_full, params
    )
    Y_full = df_dataset[model_name_full + ".Score"]

    # process "missing data" and calculate threshold for full model
    missing_df[model_name_full + ".Score"] = 0
    Y_true_missing = missing_df["Regulated"].values.astype(np.int64)
    Y_pred_missing = missing_df[model_name_full + ".Score"]
    Y_true_all = np.concatenate((Y_true, Y_true_missing))
    Y_full_all = np.concatenate((Y_full, Y_pred_missing))
    thresh_full = threshold_70_pct_recall(Y_true_all, Y_full_all)

    # get Y_true and Y_full for distance bins
    df_bin1 = df_dataset.loc[df_dataset["distanceToTSS"] < 10e3]
    Y_true_bin1 = np.concatenate(
        (df_bin1["Regulated"].values.astype(np.int64), Y_true_missing)
    )
    Y_full_bin1 = np.concatenate([df_bin1[model_name_full + ".Score"], Y_pred_missing])
    df_bin2 = df_dataset.loc[
        (df_dataset["distanceToTSS"] >= 10e3) & (df_dataset["distanceToTSS"] < 100e3)
    ]
    Y_true_bin2 = np.concatenate(
        (df_bin2["Regulated"].values.astype(np.int64), Y_true_missing)
    )
    Y_full_bin2 = np.concatenate([df_bin2[model_name_full + ".Score"], Y_pred_missing])
    df_bin3 = df_dataset.loc[df_dataset["distanceToTSS"] >= 100e3]
    Y_true_bin3 = np.concatenate(
        (df_bin3["Regulated"].values.astype(np.int64), Y_true_missing)
    )
    Y_full_bin3 = np.concatenate([df_bin3[model_name_full + ".Score"], Y_pred_missing])

    # initialize results df
    df = pd.DataFrame(
        columns=[
            "feature_ablated",
            "distance_bin",
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

    # iterate through features
    for i in range(len(feature_list)):
        print(feature_list[i])
        to_ablate = feature_list[i]

        # train model with without feature "to_ablate"
        feature_list_ablated = feature_list.loc[feature_list != to_ablate].reset_index(
            drop=True
        )
        model_name = "rE2G_without_" + to_ablate
        df_dataset = train_and_predict_once(
            df_dataset, X, Y_true, feature_list_ablated, model_name, params
        )
        Y_ablated = df_dataset[model_name + ".Score"]
        Y_ablated_all = np.concatenate([Y_ablated, Y_pred_missing])

        # calculate comparisons
        row_all_dist = get_performance_metrics(
            to_ablate,
            "all",
            Y_true_all,
            Y_full_all,
            Y_ablated_all,
            thresh_full,
            n_boot=1000,
        )

        df_bin1 = df_dataset.loc[df_dataset["distanceToTSS"] < 10e3]
        Y_ablated_bin1 = np.concatenate(
            [df_bin1[model_name + ".Score"], Y_pred_missing]
        )
        df_bin2 = df_dataset.loc[
            (df_dataset["distanceToTSS"] >= 10e3)
            & (df_dataset["distanceToTSS"] < 100e3)
        ]
        Y_ablated_bin2 = np.concatenate(
            [df_bin2[model_name + ".Score"], Y_pred_missing]
        )
        df_bin3 = df_dataset.loc[df_dataset["distanceToTSS"] >= 100e3]
        Y_ablated_bin3 = np.concatenate(
            [df_bin3[model_name + ".Score"], Y_pred_missing]
        )

        row_all_dist = get_performance_metrics(
            to_ablate,
            "all",
            Y_true_all,
            Y_full_all,
            Y_ablated_all,
            thresh_full,
            n_boot=1000,
        )
        row_bin1 = get_performance_metrics(
            to_ablate,
            "0-10 Kb",
            Y_true_bin1,
            Y_full_bin1,
            Y_ablated_bin1,
            thresh_full,
            n_boot=1000,
        )
        row_bin2 = get_performance_metrics(
            to_ablate,
            "10-100 Kb",
            Y_true_bin2,
            Y_full_bin2,
            Y_ablated_bin2,
            thresh_full,
            n_boot=1000,
        )
        row_bin3 = get_performance_metrics(
            to_ablate,
            ">100 Kb",
            Y_true_bin3,
            Y_full_bin3,
            Y_ablated_bin3,
            thresh_full,
            n_boot=1000,
        )
        df = pd.concat([df, row_all_dist, row_bin1, row_bin2, row_bin3])

    df_dataset.to_csv(
        out_dir + "/crispr_data_with_ablated_scores.tsv.gz", sep="\t", index=False
    )
    df.sort_values(by="distance_bin").to_csv(
        out_dir + "/ablation_by_feature_by_distance.tsv", sep="\t", index=False
    )
    return df


def get_performance_metrics(
    to_ablate, distance_bin, Y_true, Y_full, Y_ablated, thresh_full, n_boot=1000
):
    data_compare = (
        Y_true,
        Y_full,
        Y_ablated,
    )  #  will calculate performance(ablated - full)
    thresh_ablated = threshold_70_pct_recall(Y_true, Y_ablated)

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
        lambda Y_true, Y_full, Y_ablated: statistic_delta_precision_at_threshold(
            Y_true, Y_full, Y_ablated, thresh_full, thresh_ablated
        ),
        n_resamples=n_boot,
        paired=True,
        confidence_level=0.95,
        method="BCa",
    )
    delta_aupr = np.mean(
        res_delta_aupr.bootstrap_distribution
    )  # aupr_ablated - aupr_full
    delta_precision = np.mean(res_delta_precision.bootstrap_distribution)
    pval_aupr = bootstrap_pvalue(delta_aupr, res_delta_aupr)
    pval_precision = bootstrap_pvalue(delta_precision, res_delta_precision)

    this_row = pd.DataFrame(
        {
            "feature_ablated": to_ablate,
            "distance_bin": distance_bin,
            "delta_aupr": delta_aupr,
            "delta_aupr_low": res_delta_aupr.confidence_interval[0],
            "delta_aupr_high": res_delta_aupr.confidence_interval[1],
            "pval_aupr": pval_aupr,
            "delta_precision": delta_precision,
            "delta_precision_low": res_delta_precision.confidence_interval[0],
            "delta_precision_high": res_delta_precision.confidence_interval[1],
            "pval_precision": pval_precision,
        },
        index=[0],
    )
    return this_row


# @click.command()
# @click.option("--crispr_features_file", required=True)
# @click.option("--feature_table_file", required=True)
# @click.option("--out_dir", required=True)
# @click.option("--polynomial", type=bool, default=False)
# @click.option("--epsilon", type=float, default=0.01)
# @click.option("--n_repeats", type=float, default=20)
# @click.option("--params_file", required=True)


# def main(crispr_features_file, feature_table_file, out_dir, polynomial, epsilon, n_repeats, params_file):
def main():
    crispr_features_file = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/results/2024_0329_final_models_two_inputs/DNase_ENCFF860XAE_megamap/for_training.EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_NAfilled.tsv.gz"
    missing_file = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/results/2024_0329_final_models_two_inputs/DNase_ENCFF860XAE_megamap/missing.EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_NAfilled.tsv.gz"
    params_file = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/results/2024_0329_final_models_two_inputs/DNase_ENCFF860XAE_megamap/DNase_ENCFF860XAE_megamap_8features/model/training_params.pkl"
    out_dir = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/results/2024_0329_final_models_two_inputs/DNase_ENCFF860XAE_megamap/DNase_ENCFF860XAE_megamap_8features/ablation_analysis_by_feature"
    feature_table_file = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/results/2024_0329_final_models_two_inputs/DNase_ENCFF860XAE_megamap/feature_table.tsv"

    df_dataset = pd.read_csv(crispr_features_file, sep="\t")
    feature_table = pd.read_csv(feature_table_file, sep="\t")
    missing_df = pd.read_csv(missing_file, sep="\t")
    with open(params_file, "rb") as handle:
        params = pickle.load(handle)

    calculate_ablation(
        df_dataset, feature_table, missing_df, params, out_dir, epsilon=0.01
    )


if __name__ == "__main__":
    main()
