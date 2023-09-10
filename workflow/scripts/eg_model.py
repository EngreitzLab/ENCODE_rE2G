import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from sklearn.metrics import auc, precision_recall_curve
from sklearn.kernel_approximation import RBFSampler


CRISPR_DATA_FILE = "/oak/stanford/groups/engreitz/Users/atan5133/encode_e2g_features/results/K562_dnase_no_qnorm/EPCrisprBenchmark_ensemble_data_GRCh38.K562_ActivityOnly_features_NAfilled.tsv.gz"
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
    "ABC.Score"
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

def lr_model(training_data, training_labels):
    pipeline = Pipeline(
        [
            ("scaler", StandardScaler()),
            ("poly", PolynomialFeatures(degree=1)),
            ("clf", LogisticRegression(max_iter=100000)),
        ]
    )
    pipeline.fit(training_data, training_labels)
    return pipeline

def get_data_and_labels(data_df):
    X = data_df[E2G_FEATURES]
    X.loc[:,FEATURES_TO_LOG] = np.log10(X[FEATURES_TO_LOG] + .01)
    # X = data_df[FEATURES]

    Y = data_df["Significant"].values.astype(np.int64)
    return X, Y

def training_set_predictions(data_df):
    print("Computing predictions on training set")
    chr_list = np.unique(data_df["chrom"])
    X,Y = get_data_and_labels(data_df)
    model = lr_model(X, Y)
    predictions = model.predict_proba(X)
    data_df[TRAINING_SCORE_COL] = predictions[:, 1]

def held_out_chromosome_prediction(data_df):
    print("Computing held out chromosome predictions")
    chr_list = np.unique(data_df["chrom"])
    X,Y = get_data_and_labels(data_df)
    indices = np.arange(len(Y))
    for chr in chr_list:
        test_indices = data_df[data_df["chrom"] == chr].index.values
        if len(test_indices) > 0:
            train_indices = np.delete(indices, test_indices)
            X_test = X.loc[test_indices, :]
            X_train = X.loc[train_indices, :]
            Y_train = Y[train_indices]

            model = lr_model(X_train, Y_train)
            predictions = model.predict_proba(X_test)
            data_df.loc[test_indices, SCORE_COL] = predictions[:, 1]

def single_trained_chromosome_prediction(data_df):
    print("Computing predictions from training on 1 chromosome")
    chr_list = np.unique(data_df["chrom"])
    X,Y = get_data_and_labels(data_df)
    train_indices = data_df[data_df["chrom"] == "chr1"].index.values
    X_train = X.loc[train_indices, :]
    Y_train = Y[train_indices]

    model = lr_model(X_train, Y_train)
    predictions = model.predict_proba(X)
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
    held_out_chromosome_prediction(data_df)
    compute_AUPRC(data_df, SCORE_COL)

    training_set_predictions(data_df)
    compute_AUPRC(data_df, TRAINING_SCORE_COL)
    
    # single_trained_chromosome_prediction(data_df)
    # data_df.to_csv("predictions.tsv", sep="\t", index=False)
    # compute_AUPRC(data_df, SINGLE_CHROM_SCORE_COL)

if __name__ == "__main__":
    main()