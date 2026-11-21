import pickle
import click
import numpy as np
import pandas as pd
import shap
from sklearn.metrics import precision_recall_curve, auc, log_loss, roc_auc_score
from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from sklearn.linear_model import LogisticRegression
from training_functions import statistic_aupr
import os


def train_and_predict(
    df_dataset, feature_table, model_name, out_dir, epsilon, params, polynomial=False
):
    # specify feature list
    feature_list_core = feature_table["feature"]
    if polynomial:
        X_core = df_dataset.loc[:, feature_list_core]
        poly = PolynomialFeatures(degree=2)
        X_2 = poly.fit_transform(X_core)
        polynomial_names = poly.get_feature_names_out(X_core.columns)
        X_2 = pd.DataFrame(X_2, columns=polynomial_names)
        X = X_2
    else:
        X = df_dataset.loc[:, feature_list_core]

    # transform features
    scaler = StandardScaler()
    X = pd.DataFrame(scaler.fit_transform(X))
    # X = np.log(np.abs(X) + epsilon)
    Y = df_dataset["Regulated"].values.astype(np.int64)

    # initialize df for feature weights & metrics
    df_shap = pd.DataFrame()
    X_test_all = pd.DataFrame()

    # train aggregate model across all chromosomes
    model_full = LogisticRegression(**params).fit(X, Y)
    probs_full = model_full.predict_proba(X)
    df_dataset[model_name + ".Score_full"] = probs_full[:, 1]

    # logistic regression predictions on chromosome-wise cross validation
    idx = np.arange(len(Y))
    chr_list = np.unique(df_dataset["chr"])
    if len(chr_list) > 1:
        for chr in chr_list:
            idx_test = df_dataset[df_dataset["chr"] == chr].index.values

            if len(idx_test) > 0:
                idx_train = np.delete(idx, idx_test)
                X_test = X.loc[idx_test, :]
                X_train = X.loc[idx_train, :]
                Y_test = Y[idx_test]
                Y_train = Y[idx_train]

                model = LogisticRegression(**params).fit(X_train, Y_train)

                probs = model.predict_proba(X_test)
                df_dataset.loc[idx_test, model_name + ".Score"] = probs[:, 1]

                # shap scores
                background = shap.kmeans(
                    X_train, 10
                )  # use 10 kmeans samples from train samples as background
                explainer = shap.KernelExplainer(model.predict_proba, background)
                shap_values = explainer.shap_values(
                    X_test
                )  # d1 = eg pairs; d2 = features; d3 = neg/pos (want pos)
                print(shap_values)
                shap_values_this = pd.DataFrame(shap_values[:, :, 1])

                df_shap = pd.concat([df_shap, shap_values_this])
                X_test_all = pd.concat([X_test_all, X_test])

    # save dfs
    # df_dataset.to_csv(out_dir + '/training_predictions.tsv', sep='\t', index=False)
    # df_coef.to_csv(out_dir + '/model_coefficients.tsv', sep="\t", index=False)
    # df_metrics.to_csv(out_dir + '/performance_metrics.tsv', sep="\t", index=False)

    df_shap.columns = X_train.columns
    df_shap.to_csv(out_dir + "/shap_scores_std.tsv", sep="\t", index=False)
    X_test_all.to_csv(
        out_dir + "/training_data_in_order_std.tsv", sep="\t", index=False
    )


def main():
    dataset_dir = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/results/2024_0210_DNase_megamap_regularization/DNase_megamap/"
    model_file_name = "DNase_megamap_noReg"

    crispr_features_file = os.path.join(
        dataset_dir,
        "for_training.EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_NAfilled.tsv.gz",
    )
    feature_table_file = os.path.join(dataset_dir, "feature_table.tsv")
    model_name = "ENCODE-rE2G"
    epsilon = 0.01
    out_dir = os.path.join(dataset_dir, model_file_name, "model")
    params_file = os.path.join(
        dataset_dir, model_file_name, "model", "training_params.pkl"
    )

    df_dataset = pd.read_csv(crispr_features_file, sep="\t")
    feature_table = pd.read_csv(feature_table_file, sep="\t")
    with open(params_file, "rb") as handle:
        params = pickle.load(handle)

    train_and_predict(df_dataset, feature_table, model_name, out_dir, epsilon, params)


# @click.command()
# @click.option("--crispr_features_file", required=True)
# @click.option("--feature_table_file", required=True)
# @click.option("--out_dir", required=True)
# @click.option("--polynomial", type=bool, default=False)
# @click.option("--epsilon", type=float, default=0.01)
# @click.option("--params_file", required=True)

# def main(crispr_features_file, feature_table_file, out_dir, polynomial, params_file, epsilon):
#     model_name = 'ENCODE-rE2G'
#     df_dataset = pd.read_csv(crispr_features_file, sep="\t")
#     feature_table = pd.read_csv(feature_table_file, sep="\t")
#     with open(params_file, 'rb') as handle:
#         params = pickle.load(handle)

#     train_and_predict(df_dataset, feature_table, model_name, out_dir, epsilon, params, polynomial)

if __name__ == "__main__":
    main()
