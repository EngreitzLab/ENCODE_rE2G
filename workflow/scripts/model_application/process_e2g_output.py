import click
import pandas as pd
import os

def test_variant_overlap(all_putative, variant_overlap_output, score_column, threshold, chrom_sizes):
    variant_overlap_file = os.path.join(
	variant_overlap_output    
	)
    # generate predictions for variant overlap
    score_t = all_putative[score_column] > float(threshold)
    not_promoter = all_putative["class"] != "promoter"
    is_promoter = all_putative["class"] == "promoter"
    score_one = all_putative[score_column] > 0.1
    all_putative[(score_t & not_promoter) | (is_promoter & score_one)]
    variant_overlap = all_putative[(score_t & not_promoter) | (is_promoter & score_one)]
    # remove nan predictions
    variant_overlap_pred = variant_overlap.dropna(subset=[score_column])
    variant_overlap = variant_overlap_pred.loc[
        variant_overlap_pred["distance"] <= 2000000
    ]
    variant_overlap.to_csv(
        variant_overlap_file + ".tmp",
        sep="\t",
        index=False,
        header=True,
        compression="gzip",
        float_format="%.6f",
    )
    # shrink regions
    os.system(
        "zcat {}.tmp 2>/dev/null | head -1 | gzip > {}".format(
            variant_overlap_file, variant_overlap_file
        )
    )
    os.system(
        "zcat {}.tmp | sed 1d | bedtools slop -b -150 -g {} | gzip >> {}".format(
            variant_overlap_file, chrom_sizes, variant_overlap_file
        )
    )

    print("Done.")

def write_connections_bedpe_format(pred, outfile, score_column):
    # Output a 2d annotation file with EP connections in bedpe format for loading into IGV
    pred = pred.drop_duplicates()

    towrite = pd.DataFrame()

    towrite["chr1"] = pred["chr"]
    towrite["x1"] = pred["start"]
    towrite["x2"] = pred["end"]
    towrite["chr2"] = pred["chr"]
    towrite["y1"] = pred["TargetGeneTSS"]
    towrite["y2"] = pred["TargetGeneTSS"]
    towrite["name"] = pred["TargetGene"] + "_" + pred["name"]
    towrite["score"] = pred[score_column]
    towrite["strand1"] = "."
    towrite["strand2"] = "."

    towrite.to_csv(outfile, header=False, index=False, sep="\t")

def getE2GPairs(pred, outfile):
    e2g_df=pred[["chr","start","end", "TargetGene"]].drop_duplicates()
    e2g_df.to_csv(outfile, header=True, compression="gzip", index=False, sep="\t")

@click.command()
@click.option("--predictions_file", required=True)
@click.option("--score_column", required=True)
@click.option("--chrom_sizes", required=True)
@click.option("--threshold", required=True)
@click.option("--variant_overlap_output", required=True)
@click.option("--bedpe_output", required=True)
#@click.option("--e2g_output", required=True)

def main(predictions_file, score_column, variant_overlap_output, threshold,chrom_sizes, bedpe_output):
    all_putative = pd.read_csv(predictions_file, sep="\t")
    all_putative = all_putative.loc[all_putative.TargetGeneIsExpressed]
    test_variant_overlap(all_putative, variant_overlap_output, score_column, threshold, chrom_sizes)
    write_connections_bedpe_format(all_putative, bedpe_output, score_column)
    #getE2GPairs(all_putative, e2g_output)

if __name__ == "__main__":
    main()
