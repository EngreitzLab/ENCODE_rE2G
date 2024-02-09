import pickle
import os
import click
import numpy as np
import pandas as pd
import scipy
from training_functions import statistic_aupr, statistic_precision



def performance_summary(model_id, dataset, model_name, out_dir, n_boot=1000):
    # read in predicitons
    pred_file = os.path.join(out_dir, dataset, model_id, "model", "training_predictions.tsv")
    pred_df = pd.read_csv(pred_file, sep="\t")

    # extract relevant data
    Y_true = pred_df['Regulated'].values.astype(np.int64)
    Y_pred = pred_df[model_name+'.Score']
    data = (Y_true, Y_pred)

    # evaluate
    res_aupr =  scipy.stats.bootstrap(data, statistic_aupr, n_resamples=n_boot, paired=True, confidence_level=0.95, method='percentile')
    res_prec = scipy.stats.bootstrap(data, statistic_precision, n_resamples=n_boot, paired=True, confidence_level=0.95, method='percentile')

    res_row = pd.DataFrame({'model': model_id, 'dataset': dataset,  'AUPRC': np.mean(res_aupr.bootstrap_distribution), 'AUPRC_95CI_low': res_aupr.confidence_interval[0], 'AUPRC_95CI_high': res_aupr.confidence_interval[1],
                                 'precision': np.mean(res_prec.bootstrap_distribution), 'precision_95CI_low': res_prec.confidence_interval[0], 'precision_95CI_high': res_prec.confidence_interval[1] }, index=[0])
    return res_row

@click.command()
@click.option("--model_config_file", required=True)
@click.option("--output_file", required=True)
@click.option("--out_dir", required=True)

def main(model_config_file, output_file, out_dir):
    model_name  = 'ENCODE-rE2G'
    model_config = pd.read_table(model_config_file, na_values="").fillna("None").set_index("model", drop=False)

    # initiate final df
    df = pd.DataFrame(columns = ['model', 'dataset', 'AUPRC', 'AUPRC_95CI_low', 'AUPRC_95CI_high',
                                 'precision', 'precision_95CI_low', 'precision_95CI_high' ])
    
    # iterate through rows of model config and add results to final df
    for row in model_config.itertuples(index=False):
         res_row = performance_summary(row.model, row.dataset, model_name, out_dir)
         df = pd.concat([df, res_row])

    # sort table by AUPRC
    df = df.sort_values(by='AUPRC', ascending=False)
    df.to_csv(output_file, sep = '\t', index=False)

if __name__ == "__main__":
    main()
