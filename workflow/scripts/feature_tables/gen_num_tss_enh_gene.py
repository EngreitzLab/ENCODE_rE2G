import subprocess
import click
import pandas as pd


def determine_num_tss_enh_gene(
    pred_df, ref_gene_tss, extended_enhancers, enhancer_tss_int, out_file
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

    # bedtools intersects extended enhancers with reference TSS and counts overlaps
    header = "name\tgene\tcount\n"
    # AMANDA EDIT: extended_enhancers still has the pandas-written header row (see to_csv above).
    # `bedtools intersect -a` was reading that row as data, which errors out ("unable to
    # determine types for file") since "start"/"new_end" aren't integers. Piping through
    # `sed '1d'` strips it first, same as the pre-34de2d2 implementation did. Also switched
    # from os.system to subprocess.run(check=True) because os.system swallowed that error
    # silently -- bedtools failed on every call, so out_file always ended up with just the
    # header line and 0 data rows, and downstream activity_only_features.R crashed on the
    # empty table while Snakemake reported it as a generic (memory-related-looking) failure.
    cmd = f"""
        printf "{header}" > {out_file};
        sed '1d' {extended_enhancers} \
        | bedtools intersect -a stdin -b {ref_gene_tss} -wa -c \
        | cut -f4,5,6 \
        >> {out_file}
    """
    subprocess.run(cmd, shell=True, check=True)

    print("Saved num TSS between enh and gene")


@click.command()
@click.option("--abc_predictions")
@click.option("--ref_gene_tss")
@click.option("--extended_enhancers")
@click.option("--enhancer_tss_int")
@click.option("--out_file")
def main(abc_predictions, ref_gene_tss, extended_enhancers, enhancer_tss_int, out_file):
    pred_df = pd.read_csv(abc_predictions, sep="\t")
    if len(pred_df) == 0:
        raise Exception("Did not find any enhancers in the Predictions file")

    determine_num_tss_enh_gene(
        pred_df, ref_gene_tss, extended_enhancers, enhancer_tss_int, out_file
    )


if __name__ == "__main__":
    main()
