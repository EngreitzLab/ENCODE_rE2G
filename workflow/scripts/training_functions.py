import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_curve, auc
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LogisticRegression


## statistic functions for delta auPR/precision to be used for scipy.stats.bootstrap

def statistic_delta_aupr(y_true, y_pred_full, y_pred_ablated): # return aupr_ablated-aupr_full
    precision_full, recall_full, thresholds_full = precision_recall_curve(y_true, y_pred_full)
    aupr_full = auc(recall_full, precision_full)

    precision_ablated, recall_ablated, thresholds_ablated = precision_recall_curve(y_true, y_pred_ablated)
    aupr_ablated = auc(recall_ablated, precision_ablated)
    delta_aupr = aupr_ablated - aupr_full

    return delta_aupr

def statistic_aupr(y_true, y_pred_full):
    precision_full, recall_full, thresholds_full = precision_recall_curve(y_true, y_pred_full)
    aupr_full = auc(recall_full, precision_full)

    return aupr_full

def statistic_delta_precision(y_true, y_pred_full, y_pred_ablated): # return precision_ablated-precision_full
    precision_full, recall_full, thresholds_full = precision_recall_curve(y_true, y_pred_full)
    idx_recall_full_70_pct = np.argsort(np.abs(recall_full - 0.7))[0]
    precision_full_at_70_pct_recall = precision_full[idx_recall_full_70_pct]

    precision_ablated, recall_ablated, thresholds_ablated = precision_recall_curve(y_true, y_pred_ablated)
    idx_recall_ablated_70_pct = np.argsort(np.abs(recall_ablated - 0.7))[0]
    precision_ablated_at_70_pct_recall = precision_ablated[idx_recall_ablated_70_pct]
    delta_precision = precision_ablated_at_70_pct_recall - precision_full_at_70_pct_recall

    return delta_precision

def statistic_precision(y_true, y_pred_full):
    precision_full, recall_full, thresholds_full = precision_recall_curve(y_true, y_pred_full)
    idx_recall_full_70_pct = np.argsort(np.abs(recall_full - 0.7))[0]
    precision_full_at_70_pct_recall = precision_full[idx_recall_full_70_pct]
    return precision_full_at_70_pct_recall

# bootstrap p-values for delta (aupr/precision) 
def bootstrap_pvalue(delta, res_delta):  
    # Generate boostrap distribution of delta under null hypothesis (important centering step to get sampling distribution under the null)
    delta_boot_distribution = res_delta.bootstrap_distribution - res_delta.bootstrap_distribution.mean()

    # Calculate proportion of bootstrap samples with at least as strong evidence against null    
    pval = np.mean(np.abs(delta_boot_distribution) >= np.abs(delta))

    return pval


# assumes necessary features are present in X and Y, features are already transformed
def train_and_predict_once(df_dataset, X, Y, feature_list, model_name):
    X = X.loc[:, feature_list]
    idx = np.arange(len(Y)) # number of elements 

    chr_list = np.unique(df_dataset['chr'])
    for chr in chr_list:
        idx_test = df_dataset[df_dataset['chr']==chr].index.values
        if len(idx_test) > 0:
            idx_train = np.delete(idx, idx_test)
            X_test = X.loc[idx_test, :]
            X_train = X.loc[idx_train, :]
            Y_train = Y[idx_train]
            model = LogisticRegression(random_state=0, class_weight=None, solver='lbfgs', max_iter=100000).fit(X_train, Y_train)
            probs = model.predict_proba(X_test) # calculate scores
            df_dataset.loc[idx_test, model_name+'.Score'] = probs[:,1]

    return df_dataset