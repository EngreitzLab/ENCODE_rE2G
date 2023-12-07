import warnings

import numpy as np
import pandas as pd
from sklearn.kernel_approximation import RBFSampler
from sklearn.linear_model import LinearRegression
from sklearn.metrics import auc, precision_recall_curve
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import PolynomialFeatures, StandardScaler

warnings.filterwarnings("ignore")


CRISPR_DATA_FILE = "/oak/stanford/groups/engreitz/Users/atan5133/encode_e2g_features/results/K562_dnase_no_qnorm/EPCrisprBenchmark_ensemble_data_GRCh38.K562_ActivityOnly_features_NAfilled.tsv.gz"
TPM_FILE = (
    "/oak/stanford/groups/engreitz/Users/atan5133/encode_e2g_features/rna_tpm.tsv"
)
FEATURES = ["DHS.RPKM", "3DContactAvgHicTrack2"]
E2G_FEATURES = [
    "numTSSEnhGene",
    "distanceToTSS",
    "normalizedDNase_enh",
    "numNearbyEnhancers",
    "sumNearbyEnhancers",
    "ubiquitousExpressedGene",
    "numCandidateEnhGene",
    "3DContactAvgHicTrack2",
    "ABC.Score",
    # "TPM",
]
FEATURES_TO_LOG = [
    # "numTSSEnhGene",
    "distanceToTSS",
    # "normalizedDNase_enh",
    # "numNearbyEnhancers",
    # "sumNearbyEnhancers",
    # "ubiquitousExpressedGene",
    # "numCandidateEnhGene",
    "3DContactAvgHicTrack2",
    # "ABC.Score"
]
SCORE_COL = "Prediction.Score"
TRAINING_SCORE_COL = "Prediction.Training.Score"
SINGLE_CHROM_SCORE_COL = "Prediction.Single.Score"


def linear_regression_model(training_data, training_labels):
    pipeline = Pipeline(
        [
            ("scaler", StandardScaler()),
            ("poly", PolynomialFeatures(degree=1)),
            ("linear_regression", LinearRegression()),
        ]
    )
    pipeline.fit(training_data, training_labels)
    return pipeline


def get_data_and_labels(data_df):
    X = data_df[E2G_FEATURES]
    X.loc[:, FEATURES_TO_LOG] = np.log10(X[FEATURES_TO_LOG] + 0.01)
    # Y = np.log(1.000001 - data_df["pValueAdjusted"])
    Y = 1 - data_df["pValueAdjusted"]
    return X, Y


def training_set_predictions(data_df):
    print("Computing predictions on training set")
    chr_list = np.unique(data_df["chrom"])
    X, Y = get_data_and_labels(data_df)
    model = linear_regression_model(X, Y)
    predictions = model.predict(X)
    data_df[TRAINING_SCORE_COL] = predictions


def held_out_chromosome_prediction(data_df):
    print("Computing held out chromosome predictions")
    chr_list = np.unique(data_df["chrom"])
    X, Y = get_data_and_labels(data_df)
    indices = np.arange(len(Y))
    for chr in chr_list:
        test_indices = data_df[data_df["chrom"] == chr].index.values
        if len(test_indices) > 0:
            train_indices = np.delete(indices, test_indices)
            X_test = X.loc[test_indices, :]
            X_train = X.loc[train_indices, :]
            Y_train = Y[train_indices]

            model = linear_regression_model(X_train, Y_train)
            predictions = model.predict(X_test)
            data_df.loc[test_indices, SCORE_COL] = predictions


def single_trained_chromosome_prediction(data_df):
    print("Computing predictions from training on 1 chromosome")
    chr_list = np.unique(data_df["chrom"])
    X, Y = get_data_and_labels(data_df)
    train_indices = data_df[data_df["chrom"] == "chr1"].index.values
    X_train = X.loc[train_indices, :]
    Y_train = Y[train_indices]

    model = linear_regression_model(X_train, Y_train)
    predictions = model.predict(X)
    data_df[SINGLE_CHROM_SCORE_COL] = predictions[:, 1]

    # Y_indices = np.delete(np.arange(len(Y)), train_indices)
    # Y_true = data_df["Significant"].iloc[Y_indices].values.astype(np.int64)
    # Y_pred = data_df[SINGLE_CHROM_SCORE_COL].iloc[Y_indices]
    # precision, recall, thresholds = precision_recall_curve(Y_true, Y_pred)
    # auprc = auc(recall, precision)
    # print(f"{SINGLE_CHROM_SCORE_COL}:\nPrecision: {precision}\nRecall:{recall}\nAUPRC: {auprc}\n")


def compute_AUPRC(data_df, score_column):
    Y_true = data_df["Significant"].values.astype(np.int64)
    Y_pred = data_df[score_column]
    precision, recall, thresholds = precision_recall_curve(Y_true, Y_pred)
    auprc = auc(recall, precision)
    print(f"{score_column}:\nPrecision: {precision}\nRecall:{recall}\nAUPRC: {auprc}\n")


def main():
    data_df = pd.read_csv(CRISPR_DATA_FILE, sep="\t")
    # Only consider negative effect size samples
    data_df = data_df[
        (~(((data_df["EffectSize"] > 0) & (data_df["Significant"]))))
        & (~data_df["pValueAdjusted"].isna())
    ].reset_index(drop=True)
    print(f"Working with {len(data_df)} crispr samples")

    tpm_df = pd.read_csv(TPM_FILE, sep="\t", index_col=0)
    # merge DF based on measuredGeneSymbol and index
    data_df = data_df.merge(
        tpm_df, how="left", left_on="measuredGeneSymbol", right_index=True
    )
    data_df["TPM"].fillna(0, inplace=True)
    held_out_chromosome_prediction(data_df)
    compute_AUPRC(data_df, SCORE_COL)

    training_set_predictions(data_df)
    compute_AUPRC(data_df, TRAINING_SCORE_COL)

    # single_trained_chromosome_prediction(data_df)
    # data_df.to_csv("predictions.tsv", sep="\t", index=False)
    # compute_AUPRC(data_df, SINGLE_CHROM_SCORE_COL)


if __name__ == "__main__":
    main()
