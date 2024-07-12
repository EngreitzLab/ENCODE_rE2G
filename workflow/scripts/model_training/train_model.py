import pickle
import click
import numpy as np
import pandas as pd
import shap
from sklearn.metrics import precision_recall_curve, auc, log_loss, roc_auc_score
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LogisticRegression
from training_functions import statistic_aupr


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
    X = np.log(np.abs(X) + epsilon)
    Y = df_dataset["Regulated"].values.astype(np.int64)

    # initialize df for feature weights & metrics
    df_coef = pd.DataFrame(columns=["feature", "coefficient", "test_chr"])
    df_metrics = pd.DataFrame(
        columns=[
            "test_chr",
            "log_loss_test_full",
            "log_loss_train",
            "log_loss_test",
            "AUROC_test_full",
            "AUROC_train",
            "AUROC_test",
            "AUPRC_test_full",
            "AUPRC_train",
            "AUPRC_test",
            "n_test_pos",
            "n_test_neg",
            "n_train_neg",
            "n_train_pos",
        ]
    )
    # df_shap = pd.DataFrame()
    # X_test_all = pd.DataFrame()

    # train aggregate model across all chromosomes, calc weights and performance metrics, save full model
    model_full = LogisticRegression(**params).fit(X, Y)
    probs_full = model_full.predict_proba(X)
    df_dataset[model_name + ".Score_full"] = probs_full[:, 1]
    coefficients = model_full.coef_[0]
    df_temp = pd.DataFrame(
        {"feature": X.columns, "coefficient": coefficients, "test_chr": "none"}
    )
    df_coef = pd.concat([df_coef, df_temp])
    with open(out_dir + f"/model_full.pkl", "wb") as f:
        pickle.dump(model_full, f)

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

                with open(out_dir + f"/model_test_{chr}.pkl", "wb") as f:
                    pickle.dump(model, f)

                probs = model.predict_proba(X_test)
                df_dataset.loc[idx_test, model_name + ".Score"] = probs[:, 1]

                # performance metrics
                n_train_pos = np.sum(Y_train)
                n_train_neg = len(Y_train) - n_train_pos
                n_test_pos = np.sum(Y_test)
                n_test_neg = len(Y_test) - n_test_pos

                ll_train = log_loss(Y_train, model.predict_proba(X_train)[:, 1])
                ll_test_full = (
                    log_loss(Y_test, probs_full[idx_test, 1])
                    if n_test_pos > 0
                    else np.NaN
                )
                ll_test = log_loss(Y_test, probs[:, 1]) if n_test_pos > 0 else np.NaN
                auroc_train = roc_auc_score(Y_train, model.predict_proba(X_train)[:, 1])
                auroc_test_full = (
                    roc_auc_score(Y_test, probs_full[idx_test, 1])
                    if n_test_pos > 0
                    else np.NaN
                )
                auroc_test = (
                    roc_auc_score(Y_test, probs[:, 1]) if n_test_pos > 0 else np.NaN
                )
                auprc_train = statistic_aupr(
                    Y_train, model.predict_proba(X_train)[:, 1]
                )
                auprc_test_full = (
                    statistic_aupr(Y_test, probs_full[idx_test, 1])
                    if n_test_pos > 0
                    else np.NaN
                )
                auprc_test = (
                    statistic_aupr(Y_test, probs[:, 1]) if n_test_pos > 0 else np.NaN
                )

                df_temp = pd.DataFrame(
                    {
                        "test_chr": [chr],
                        "log_loss_test_full": [ll_test_full],
                        "log_loss_train": [ll_train],
                        "log_loss_test": [ll_test],
                        "AUROC_test_full": [auroc_test_full],
                        "AUROC_train": [auroc_train],
                        "AUROC_test": [auroc_test],
                        "AUPRC_test_full": [auprc_test_full],
                        "AUPRC_train": [auprc_train],
                        "AUPRC_test": [auprc_test],
                        "n_test_pos": [n_test_pos],
                        "n_test_neg": [n_test_neg],
                        "n_train_pos": [n_train_pos],
                        "n_train_neg": [n_train_neg],
                    }
                )
                df_metrics = pd.concat([df_metrics, df_temp])

                # save model weights
                coefficients = model.coef_[0]
                df_temp = pd.DataFrame(
                    {"feature": X.columns, "coefficient": coefficients, "test_chr": chr}
                )
                df_coef = pd.concat([df_coef, df_temp])

                # shap scores
                # background = shap.kmeans(X_train, 10)  # use 10 kmeans samples from train samples as background
                # explainer = shap.KernelExplainer(model.predict_proba, background)
                # shap_values = explainer.shap_values(X_test) # d1 = eg pairs; d2 = features; d3 = neg/pos (want pos)
                # shap_values_this =pd.DataFrame(shap_values[:,:,1])

                # df_shap = pd.concat([df_shap, shap_values_this])
                # X_test_all = pd.concat([X_test_all, X_test])

        # calc performance metrics across chromosomes
        total_ll_test = log_loss(Y, df_dataset[model_name + ".Score"])
        total_ll_test_full = log_loss(Y, probs_full[:, 1])
        total_auroc_test = roc_auc_score(Y, df_dataset[model_name + ".Score"])
        total_auroc_test_full = log_loss(Y, probs_full[:, 1])
        total_auprc_test = statistic_aupr(Y, df_dataset[model_name + ".Score"])
        total_auprc_test_full = statistic_aupr(Y, probs_full[:, 1])
        n_pos = np.sum(Y)
        n_neg = len(Y) - n_pos

        df_temp = pd.DataFrame(
            {
                "test_chr": ["all"],
                "log_loss_test_full": [total_ll_test_full],
                "log_loss_train": [np.NaN],
                "log_loss_test": [total_ll_test],
                "AUROC_test_full": [total_auroc_test_full],
                "AUROC_train": [np.NaN],
                "AUROC_test": [total_auroc_test],
                "AUPRC_test_full": [total_auprc_test_full],
                "AUPRC_train": [np.NaN],
                "AUPRC_test": [total_auprc_test],
                "n_test_pos": [n_pos],
                "n_test_neg": [n_neg],
                "n_train_pos": [np.NaN],
                "n_train_neg": [np.NaN],
            }
        )
        df_metrics = pd.concat([df_metrics, df_temp])

    # save dfs
    df_dataset.to_csv(out_dir + "/training_predictions.tsv", sep="\t", index=False)
    df_coef.to_csv(out_dir + "/model_coefficients.tsv", sep="\t", index=False)
    df_metrics.to_csv(out_dir + "/performance_metrics.tsv", sep="\t", index=False)

    # df_shap.columns = X_train.columns
    # df_shap.to_csv(out_dir +'/shap_scores.tsv', sep="\t", index=False)
    # X_test_all.to_csv(out_dir + '/training_data_in_order.tsv',  sep="\t", index=False)


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

    train_and_predict(
        df_dataset, feature_table, model_name, out_dir, epsilon, params, polynomial
    )


if __name__ == "__main__":
    main()
