import os
import click
import pandas as pd


def determine_num_candidate_enh_gene(pred_df, out_file):
    # Concept: Start by making a working copy of our data so we don't accidentally
    # modify the original DataFrame that was passed into the function.
    df = pred_df.copy()

    # Concept: To determine if an enhancer is upstream or downstream, we need a single
    # point to represent it. The midpoint is a simple and effective choice.
    df["midpoint"] = ((df["start"] + df["end"]) / 2).astype("int")

    # Concept: For our counting logic to work correctly, all enhancers must be
    # sorted by their genomic position. This creates a predictable 'map' of the
    # chromosome that we can operate on. We sort by gene first, then position.
    df = df.sort_values(by=["TargetGene", "midpoint"], ascending=[True, True])

    # Concept: We need to handle enhancers before the TSS differently from those after it.
    # Let's create two 'masks'—like stencils—to quickly select all the upstream
    # enhancers or all the downstream enhancers in one go.
    is_downstream = df["midpoint"] > df["TargetGeneTSS"]
    is_upstream = df["midpoint"] < df["TargetGeneTSS"]

    # --- Upstream Calculation ---
    # Concept: To rank upstream enhancers, we need to count *away* from the TSS.
    # This means starting with the enhancer closest to the TSS (the one with the
    # highest 'midpoint' coordinate) and counting backwards.
    # We achieve this by creating a temporary, sorted copy of just the upstream enhancers.
    df.loc[is_upstream, "NumCandidateEnhGene"] = (
        # 1. Select only the upstream enhancers. This creates a temporary table.
        df[is_upstream]
        # 2. Sort this temporary table in DESCENDING order of position. Now, the
        #    enhancer closest to the TSS is at the top of each gene's group.
        .sort_values("midpoint", ascending=False)
        # 3. Group by gene and then perform a cumulative count. `cumcount` automatically
        #    creates the ranking (0, 1, 2, ...) for us within each gene's group.
        .groupby("TargetGene").cumcount()
        # 4. `cumcount` starts at 0, but we want our ranks to start at 1.
        + 1
    )

    # --- Downstream Calculation ---
    # Concept: For downstream enhancers, the default sort order (ascending) is already
    # correct for counting away from the TSS. The enhancer with the lowest
    # 'midpoint' coordinate is the one closest to the TSS.
    df.loc[is_downstream, "NumCandidateEnhGene"] = (
        # 1. Select only the downstream enhancers.
        df[is_downstream]
        # 2. Group by gene and perform the cumulative count. No special sorting needed.
        .groupby("TargetGene").cumcount()
        # 3. Add 1 to start the rank from 1 instead of 0.
        + 1
    )

    # Concept: Any enhancers located exactly at the TSS (or any other unhandled cases)
    # will have a missing ('NaN') rank. We'll fill these with 0 to be safe.
    df = df.fillna(value=0)
    df["NumCandidateEnhGene"] = df["NumCandidateEnhGene"].astype("int")

    # Concept: Before saving, restore the desired chromosomal sort order
    df = df.sort_values(by=["chr", "midpoint"], ascending=True).reset_index(drop=True)

    # Concept: Finally, select only the columns we need for the final report and
    # save the result to a file.
    df[["name", "TargetGene", "NumCandidateEnhGene"]].to_csv(
        out_file,
        sep="\t",
        index=False,
    )
    print("Saved num candidate enhancers")


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
