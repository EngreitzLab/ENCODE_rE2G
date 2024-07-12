import pickle
import os
import click
import numpy as np
import pandas as pd
import scipy
from training_functions import (
    statistic_aupr,
    statistic_precision_at_threshold,
    threshold_70_pct_recall,
)


def performance_summary(
    model_id, dataset, model_name, out_dir, crispr_data="", n_boot=1000
):
    # read in predicitons
    if model_id == "distance":
        crispr_data = pd.read_csv(crispr_data, sep="\t")
        crispr_data["distance"] = np.abs(
            (crispr_data["chromStart"] + crispr_data["chromEnd"]) / 2
            - (crispr_data["startTSS"] + crispr_data["endTSS"]) / 2
        )
        crispr_data = crispr_data.dropna(subset=["Regulated", "distance"])
        Y_true_all = crispr_data["Regulated"].values.astype(np.int64)
        Y_pred_all = crispr_data["distance"] * -1
        pct_missing = 0
    else:  # normal models
        pred_file = os.path.join(
            out_dir, dataset, model_id, "model", "training_predictions.tsv"
        )
        missing_file = os.path.join(
            out_dir,
            dataset,
            "missing.EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_NAfilled.tsv.gz",
        )
        pred_df = pd.read_csv(pred_file, sep="\t")
        missing_df = pd.read_csv(missing_file, sep="\t")

        # extract relevant data from predictions
        Y_true = pred_df["Regulated"].values.astype(np.int64)
        Y_pred = pred_df[model_name + ".Score"]

        # add rows for crispr pairs not in predictions/training data
        n_missing = len(missing_df)
        if n_missing > 0:
            missing_df[model_name + ".Score"] = 0
            Y_true_missing = missing_df["Regulated"].values.astype(np.int64)
            Y_pred_missing = missing_df[model_name + ".Score"]
            Y_true_all = np.concatenate((Y_true, Y_true_missing))
            Y_pred_all = np.concatenate((Y_pred, Y_pred_missing))
        else:
            Y_true_all = Y_true
            Y_pred_all = Y_pred
        pct_missing = n_missing / len(Y_true_all)

    # evaluate
    res_aupr = scipy.stats.bootstrap(
        (Y_true_all, Y_pred_all),
        statistic_aupr,
        n_resamples=n_boot,
        paired=True,
        confidence_level=0.95,
        method="BCa",
    )
    thresh = threshold_70_pct_recall(
        Y_true_all, Y_pred_all
    )  # will return None if max recall < 70%
    if thresh is not None:
        res_prec = scipy.stats.bootstrap(
            (Y_true_all, Y_pred_all),
            lambda Y_true, Y_pred: statistic_precision_at_threshold(
                Y_true, Y_pred, thresh
            ),
            n_resamples=n_boot,
            paired=True,
            confidence_level=0.95,
            method="BCa",
        )
    prec_mean = 0 if thresh is None else np.mean(res_prec.bootstrap_distribution)
    prec_low = 0 if thresh is None else res_prec.confidence_interval[0]
    prec_high = 0 if thresh is None else res_prec.confidence_interval[1]

    res_row = pd.DataFrame(
        {
            "model": model_id,
            "dataset": dataset,
            "AUPRC": np.mean(res_aupr.bootstrap_distribution),
            "AUPRC_95CI_low": res_aupr.confidence_interval[0],
            "AUPRC_95CI_high": res_aupr.confidence_interval[1],
            "precision": prec_mean,
            "precision_95CI_low": prec_low,
            "precision_95CI_high": prec_high,
            "threshold_70_pct_recall": thresh,
            "pct_missing_elements": pct_missing,
        },
        index=[0],
    )
    return res_row


@click.command()
@click.option("--model_config_file", required=True)
@click.option("--output_file", required=True)
@click.option("--crispr_data", required=True)
@click.option("--out_dir", required=True)
def main(model_config_file, output_file, crispr_data, out_dir):
    model_name = "ENCODE-rE2G"
    model_config = (
        pd.read_table(model_config_file, na_values="")
        .fillna("None")
        .set_index("model", drop=False)
    )
    # initiate final df
    df = pd.DataFrame(
        columns=[
            "model",
            "dataset",
            "AUPRC",
            "AUPRC_95CI_low",
            "AUPRC_95CI_high",
            "precision",
            "precision_95CI_low",
            "precision_95CI_high",
            "threshold_70_pct_recall",
            "pct_missing_elements",
        ]
    )

    # iterate through rows of model config and add results to final df
    for row in model_config.itertuples(index=False):
        res_row = performance_summary(row.model, row.dataset, model_name, out_dir)
        df = pd.concat([df, res_row])
    # add row for distance
    res_row = performance_summary(
        "distance", "baseline", "", out_dir, crispr_data=crispr_data
    )
    df = pd.concat([df, res_row])

    # sort table by AUPRC
    df = df.sort_values(by="AUPRC", ascending=False)
    df.to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":
    main()
