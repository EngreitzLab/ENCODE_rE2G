import os
import click
import pandas as pd


def determine_num_candidate_enh_gene(pred_df, out_file):
    df = pred_df.copy()
    df["midpoint"] = ((df["start"] + df["end"]) / 2).astype("int")
    df = df.sort_values(by=["TargetGene", "midpoint"], ascending=[True, True])

    is_downstream = df["midpoint"] > df["TargetGeneTSS"]
    is_upstream = df["midpoint"] < df["TargetGeneTSS"]

    # Rank upstream enhancers counting away from TSS (descending position order)
    df.loc[is_upstream, "NumCandidateEnhGene"] = (
        df[is_upstream]
        .sort_values("midpoint", ascending=False)
        .groupby("TargetGene").cumcount()
        + 1
    )

    # Rank downstream enhancers counting away from TSS (ascending order is already correct)
    df.loc[is_downstream, "NumCandidateEnhGene"] = (
        df[is_downstream]
        .groupby("TargetGene").cumcount()
        + 1
    )

    # Enhancers exactly at the TSS get rank 0
    df = df.fillna(value=0)
    df["NumCandidateEnhGene"] = df["NumCandidateEnhGene"].astype("int")

    df = df.sort_values(by=["chr", "midpoint"], ascending=True).reset_index(drop=True)
    df[["name", "TargetGene", "NumCandidateEnhGene"]].to_csv(
        out_file,
        sep="\t",
        index=False,
    )


@click.command()
@click.option("--abc_predictions")
@click.option("--out_file")
def main(abc_predictions, out_file):
    pred_df = pd.read_csv(abc_predictions, sep="\t", compression="gzip")
    if len(pred_df) == 0:
        raise Exception("Did not find any enhancers in the Predictions file")

    determine_num_candidate_enh_gene(pred_df, out_file)


if __name__ == "__main__":
    main()
