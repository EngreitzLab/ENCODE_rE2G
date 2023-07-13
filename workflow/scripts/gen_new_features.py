import os

import click
import numpy as np
import pandas as pd


def create_intermediate_dir(results_dir):
    intermediate_dir = os.path.join(results_dir, "intermediate_files")
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    return intermediate_dir


@click.command()
@click.option("--enhancer_list")
@click.option("--ref_gene_tss")
@click.option("--abc_predictions")
@click.option("--results_dir")
def main(enhancer_list, ref_gene_tss, abc_predictions, results_dir):
    intermediate_dir = create_intermediate_dir(results_dir)
    # Read in input prediction file
    subset = pd.read_csv(
        enhancer_list, sep="\t", usecols=["chr", "start", "end", "name"]
    )
    ##### columns in input prediction file: chr	start	end	name	class	activity_base	TargetGene	TargetGeneTSS	TargetGeneExpression	TargetGenePromoterActivityQuantile	TargetGeneIsExpressed	distance	isSelfPromoter	powerlaw_contact	powerlaw_contact_reference	hic_contact	hic_contact_pl_scaled	hic_pseudocount	hic_contact_pl_scaled_adj	ABC.Score.Numerator	ABC.Score	powerlaw.Score.Numerator	powerlaw.Score	CellType
    #### subsetting data to just columns required
    # subset = data[['chr', 'start', 'end', 'name', 'class', 'activity_base', 'TargetGene', 'TargetGeneTSS', 'TargetGeneExpression', 'TargetGenePromoterActivityQuantile', 'TargetGeneIsExpressed', 'distance', 'hic_contact', 'powerlaw_contact', 'hic_contact_pl_scaled', 'ABC.Score.Numerator', 'ABC.Score', 'CellType']]
    ##### remove entries that are promoters
    subset = subset.loc[subset["class"] != "promoter"]
    if len(subset) > 1:
        ##### Get midpoint of enhancer region
        subset["midpoint"] = subset["start"] + 0.5 * (subset["end"] - subset["start"])
        subset["midpoint"] = subset["midpoint"].astype("int")

        ##### Make the end be midpoint of enhancer + distance (This gives you the end coordinate of distance range)
        subset["end1"] = subset["midpoint"] + subset["distance"]

        ## If gene is located upstream of enhancer, modify the start to be the beginning of the TargetGeneTSS and the end be the midpoint of the enhancer
        entries = subset.loc[subset["TargetGeneTSS"] < subset["midpoint"]].index.astype(
            "int"
        )
        subset.loc[entries, "end1"] = subset.loc[entries, "end"]
        subset.loc[entries, "start"] = subset.loc[entries, "TargetGeneTSS"]
        subset["start"] = subset["start"].astype("int")
        subset["end1"] = subset["end1"].astype("int")
        # This file will be used to intersect gene annotations to count How many protein-coding TSSs away is the enhancer from the promoter?  (i.e., how many protein-coding gene TSSs are located between the enhancer and promoter?  0 = closest TSS)
        # Intersect this file with gene TSS
        extended_enhregions = os.path.join(
            intermediate_dir, "EnhancerRegions_extended.txt"
        )
        subset[["chr", "start", "end1", "name", "TargetGene"]].to_csv(
            extended_enhregions, sep="\t", index=False
        )

        ## Intersect midpoint of enhancer to target gene regions with GeneTSS
        ## This will include overlaps with the TargetGene
        filename = os.path.join(
            intermediate_dir, "EnhancerRegions_extended_RefSeq.intersected.txt.gz"
        )
        os.system(
            "sed '1d' {} | bedtools intersect -a stdin -b {} -wa -wb | cut -f4,5 | gzip > {}".format(
                extended_enhregions, ref_gene_tss, filename
            )
        )
        ##
        candidate_enhancers = pd.read_csv(enhancer_list, sep="\t")
        candidate_enhancers = candidate_enhancers.loc[
            candidate_enhancers["class"] != "promoter"
        ]
        candidate_enhancersfile = os.path.join(
            intermediate_dir, "EnhancerList_noPromoter.txt"
        )
        candidate_enhancers[["chr", "start", "end"]].to_csv(
            candidate_enhancersfile, sep="\t", index=False
        )
        #### How many candidate enhancer regions are located between the enhancer and the promoter?
        filename_enh = os.path.join(
            intermediate_dir, "EnhancerRegions_extended_CandidateReg.txt.gz"
        )
        os.system(
            "sed '1d' {} | bedtools intersect -a stdin -b {} -wa -wb | cut -f4,5 | gzip > {}".format(
                extended_enhregions, candidate_enhancersfile, filename_enh
            )
        )
        print("Reading in {}".format(filename_enh))
        predictions_can = pd.read_csv(filename_enh, sep="\t", names=["class", "gene"])
        num_candidate_enhancers = (
            predictions_can.groupby(["class", "gene"]).size().reset_index()
        )

        num_candidate_enhancers.to_csv(
            os.path.join(results_dir, "NumCandidateEnhGene.txt"),
            sep="\t",
            index=False,
        )
        print("Saved num candidate enhancers")

        print("Reading in {}".format(filename))
        ## Read in intersected midpoint of enhancer to target gene regions with GeneTSS file
        predictions = pd.read_csv(filename, sep="\t", names=["class", "gene"])
        # Calculate the number of TSS regions that fall within the enhancer to target gene regions.
        num_tss_between_enh_and_gene = (
            predictions.groupby(["class", "gene"]).size().reset_index()
        )
        num_tss_between_enh_and_gene.to_csv(
            os.path.join(results_dir, "NumTSSEnhGene.txt"),
            sep="\t",
            index=False,
        )
        print("Saved num TSS between enh and gene")

        ############ Generate Num/Sum Enhancers within 5kb/10kb ############
        ##### Make the end be midpoint of enhancer + distance (This gives you the end coordinate of distance range)
        enh_list_midpoint_file = os.path.join(
            intermediate_dir, "EnhancerList_midpoint.bed"
        )
        subset[["chr", "midpoint", "midpoint", "name"]].to_csv(
            enh_list_midpoint_file,
            sep="\t",
            index=False,
            header=False,
        )
        slop_5kb_file = os.path.join(intermediate_dir, "EnhancerRegions.slop_5kb.bed")
        os.system(
            "bedtools slop -b 5000 -i {} -g /oak/stanford/groups/akundaje/kmualim/ABC-Enhancer-Gene-Prediction/sherlock_scripts/hg38.chrom.sizes > {}".format(
                enh_list_midpoint_file, slop_5kb_file
            )
        )
        slop_10kb_file = os.path.join(intermediate_dir, "EnhancerRegions.slop_10kb.bed")
        os.system(
            "bedtools slop -b 10000 -i {} -g /oak/stanford/groups/akundaje/kmualim/ABC-Enhancer-Gene-Prediction/sherlock_scripts/hg38.chrom.sizes > {}".format(
                enh_list_midpoint_file, slop_10kb_file
            )
        )
        prediction_slim_file = os.path.join(
            intermediate_dir, "EnhancerPredictionsAllPutative.slim.bed"
        )
        os.system(
            "zcat {} | cut -f1,2,3,4,6 | sed '1d' > {}".format(
                abc_predictions, prediction_slim_file
            )
        )
        os.system(
            "bedtools intersect -a {} -b {} -wa -wb > {}".format(
                slop_5kb_file,
                prediction_slim_file,
                os.path.join(intermediate_dir, "NumEnhancers5kb.txt"),
            )
        )
        os.system(
            "bedtools intersect -a {} -b {} -wa -wb > {}".format(
                slop_10kb_file,
                prediction_slim_file,
                os.path.join(intermediate_dir, "NumEnhancers10kb.txt"),
            )
        )
        data = pd.read_csv(
            os.path.join(intermediate_dir, "NumEnhancers5kb.txt"),
            sep="\t",
            header=None,
        )
        data1 = data.loc[data[3] != data[7]]
        data2 = data1.groupby([3]).size().reset_index(name="count")
        data3 = data1.groupby([3])[8].sum().reset_index(name="sum")
        data2.to_csv(
            os.path.join(results_dir, "NumEnhancersEG5kb.txt"),
            sep="\t",
            header=False,
            index=False,
        )
        data3.to_csv(
            os.path.join(results_dir, "SumEnhancersEG5kb.txt"),
            sep="\t",
            header=False,
            index=False,
        )
        data = pd.read_csv(
            os.path.join(intermediate_dir, "NumEnhancers10kb.txt"),
            sep="\t",
            header=None,
        )
        data1 = data.loc[data[3] != data[7]]
        data2 = data1.groupby([3]).size().reset_index(name="count")
        data3 = data1.groupby([3])[8].sum().reset_index(name="sum")
        data2.to_csv(
            os.path.join(results_dir, "NumEnhancersEG10kb.txt"),
            sep="\t",
            header=False,
            index=False,
        )
        data3.to_csv(
            os.path.join(results_dir, "SumEnhancersEG10kb.txt"),
            sep="\t",
            header=False,
            index=False,
        )


if __name__ == "__main__":
    main()
