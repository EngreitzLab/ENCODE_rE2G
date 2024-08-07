import numpy as np
import subprocess
from io import StringIO
import click
import pandas as pd

NORMAL_CHROMOSOMES = set(["chr" + str(x) for x in range(1, 23)] + ["chrX"] + ["chrY"])


def count_bam_total(bam_file: str) -> int:
    cmd = ["samtools", "idxstat", bam_file]
    result = subprocess.check_output(cmd).decode("utf-8")
    tsv_io = StringIO(result)
    df = pd.read_csv(
        tsv_io, sep="\t", names=["chr", "size", "mapped_reads", "unmapped_reads"]
    )
    no_alt_chrom_df = df[df["chr"].isin(NORMAL_CHROMOSOMES)]
    return no_alt_chrom_df["mapped_reads"].sum()


def get_num_reads(accessibility_files):
    total_counts = 0
    for access_in in accessibility_files:
        if not access_in.endswith(".bam"):
            print("Only support num reads for bam files")
            return 0
        total_counts += count_bam_total(access_in)
    return total_counts


def get_num_enh(df):
    return len(df[["chr", "start", "end"]].drop_duplicates())


def get_num_genes(df):
    return len(df["TargetGene"].drop_duplicates())


def get_num_enh_gene_links(df):
    return len(df)


def get_num_genes_with_1_enh_min(df):
    genes = df.groupby(["TargetGene"]).size()
    return len(genes[genes > 0])


def get_mean_num_genes_per_enh(df):
    enh = df[["chr", "start", "end"]].groupby(["chr", "start", "end"]).size()
    return enh.mean()


def get_mean_num_enh_per_gene(df):
    genes = df.groupby(["TargetGene"]).size()
    return genes.mean()


def get_mean_num_enh_per_gene_no_prom(df):
    df = df[df["class"] != "promoter"]
    genes = df.groupby(["TargetGene"]).size()
    return genes.mean()


def get_mean_log_dist_to_tss(df):
    log_dist = df["distanceToTSS"].apply(np.log10)
    log_dist = log_dist.replace(-np.inf, 0)
    return log_dist.mean()


def get_mean_enh_region_size(df):
    enh = df.groupby(["chr", "start", "end"])
    sizes = []
    for (_, start, end), _ in enh:
        sizes.append(end - start)
    sizes = pd.Series(sizes)
    return sizes.mean()


@click.command()
@click.option("--predictions", type=str, required=True)
@click.option("--accessibility", type=str, required=True)
@click.option("--output_file", type=str, default="stats.tsv")
def main(predictions, accessibility, output_file):
    df = pd.read_csv(predictions, sep="\t")
    accessibility_files = [f.strip() for f in accessibility.split(" ")]
    stats = [
        ("num_sequencing_reads", get_num_reads(accessibility_files)),
        ("num_enh", get_num_enh(df)),
        ("num_genes", get_num_genes(df)),
        ("num_enh_gene_links", get_num_enh_gene_links(df)),
        ("num_genes_with_1_enh_min", get_num_genes_with_1_enh_min(df)),
        ("mean_num_genes_per_enh", get_mean_num_genes_per_enh(df)),
        ("mean_num_enh_per_gene", get_mean_num_enh_per_gene(df)),
        ("mean_num_enh_per_gene_no_prom", get_mean_num_enh_per_gene_no_prom(df)),
        ("mean_log10_dist_to_tss", get_mean_log_dist_to_tss(df)),
        ("mean_enh_region_size", get_mean_enh_region_size(df)),
    ]
    metric = [stat[0] for stat in stats]
    values = [stat[1] for stat in stats]
    df = pd.DataFrame({"Metric": metric, "Value": values})
    df.to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":
    main()
