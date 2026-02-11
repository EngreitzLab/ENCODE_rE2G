import os
import click
import pandas as pd


def determine_num_tss_enh_gene(
    pred_df, ref_gene_tss, extended_enhancers, out_file  # enhancer_tss_int
):
    #  make the end be midpoint of enhancer + distance (This gives you the end coordinate of distance range)
    pred_df["midpoint"] = ((pred_df["start"] + pred_df["end"]) / 2).astype("int")
    pred_df["new_end"] = (pred_df["midpoint"] + pred_df["distance"]).astype("int")

    # if gene is located upstream of enhancer, modify the start to be the beginning of the TargetGeneTSS and the end be the midpoint of the enhancer
    downstream_enh = pred_df[pred_df["TargetGeneTSS"] < pred_df["midpoint"]]
    pred_df.loc[downstream_enh.index, "new_end"] = pred_df.loc[
        downstream_enh.index, "end"
    ]
    pred_df.loc[downstream_enh.index, "start"] = pred_df.loc[
        downstream_enh.index, "TargetGeneTSS"
    ]

    #  File to intersect with TSS annotations to count #  protein-coding TSSs between enhancer and promoter  (0 = closest TSS)
    pred_df[["chr", "start", "new_end", "name", "TargetGene"]].to_csv(
        extended_enhancers, sep="\t", index=False
    )

    # This pipeline performs the intersection, counting, and saves the final small result.
    # It is the command-line equivalent of groupby().size() and is very memory-efficient.
    header = "name\tgene\tcount\n"
    cmd = f"""
        printf "{header}" > {out_file};
        bedtools intersect -a {extended_enhancers} -b {ref_gene_tss} -wa -c \
        | cut -f4,5,6 \
        >> {out_file}
    """
    os.system(cmd)

    print("Saved num TSS between enh and gene")


@click.command()
@click.option("--abc_predictions")
@click.option("--ref_gene_tss")
@click.option("--extended_enhancers")
@click.option("--out_file")
def main(abc_predictions, ref_gene_tss, extended_enhancers, out_file):
    pred_df = pd.read_csv(abc_predictions, sep="\t", compression="gzip")
    if len(pred_df) == 0:
        raise Exception("Did not find any enhancers in the Predictions file")

    determine_num_tss_enh_gene(pred_df, ref_gene_tss, extended_enhancers, out_file)


if __name__ == "__main__":
    main()
