import pickle

import click
import numpy as np
import pandas as pd

MODEL = "ENCODE-E2G"


def threshold_predictions(all_putative, threshold, keep_self_promoters):
    # filter by threshold and promoter category
    filtered_predictions = all_putative[all_putative['ENCODE-E2G.Score'] >= threshold]

    # process promoters & self-promoters
    if keep_self_promoters:
        return filtered_predictions[(filtered_predictions["class"] != "promoter") | filtered_predictions["isSelfPromoter"]]
    else:
        return filtered_predictions[filtered_predictions["class"] != "promoter"]
            
    return filtered_predictions


@click.command()
@click.option("--all_predictions_file", required=True)
@click.option("--threshold", required=True)
@click.option("--keep_self_promoters", type=bool, default=True)
@click.option("--output_file", required=True)

def main(all_predictions_file, threshold, keep_self_promoters, output_file):
    all_predictions = pd.read_csv(all_predictions_file, sep="\t")

    filtered_predictions = threshold_predictions(all_predictions, threshold, keep_self_promoters)

    filtered_predictions.to_csv(output_file, compression="gzip", sep="\t", index=False)

if __name__ == "__main__":
    main()
