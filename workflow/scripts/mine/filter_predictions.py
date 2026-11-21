import click
import pandas as pd


def threshold_predictions(
    pred_file, out_file, threshold, score_col, include_self_promoter
):
    # filter by threshold and promoter category
    all_putative = pd.read_csv(pred_file, sep="\t")
    filtered_predictions = all_putative[all_putative[score_col] >= float(threshold)]

    # process promoters & self-promoters
    if include_self_promoter:
        filtered_predictions = filtered_predictions[
            (filtered_predictions["class"] != "promoter")
            | filtered_predictions["isSelfPromoter"]
        ]
    else:
        filtered_predictions = filtered_predictions[
            filtered_predictions["class"] != "promoter"
        ]

    filtered_predictions.to_csv(out_file, compression="gzip", sep="\t", index=False)


@click.command()
@click.option("--pred_file", required=True)
@click.option("--out_file", required=True)
@click.option("--threshold", required=True)
@click.option("--score_col", required=True)
def main(pred_file, out_file, threshold, score_col):
    threshold_predictions(
        pred_file, out_file, threshold, score_col, include_self_promoter=True
    )


if __name__ == "__main__":
    main()
