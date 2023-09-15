import click
import scanpy
import pandas as pd

COLS = ["chr", "start", "end", "name"]

def get_gene_lengths(gene_ref_bed):
    # lengths in KB
    df = pd.read_table(
        gene_ref_bed, names=COLS, header=None, comment="#", usecols=COLS
    )
    gene_lengths = {}
    for _, row in df.iterrows():
        gene_lengths[row["name"]] = (row["end"] - row["start"]) / 1e3
    return gene_lengths

def get_gene_rna_counts(rna_10x):
    h5 = scanpy.read_10x_h5(rna_10x)
    df = h5.to_df().transpose()
    df["Counts"] = df.sum(axis=1)
    return df[["Counts"]]

def filter_unknown_genes(rna_count_df, gene_lengths):
    filtered_df = rna_count_df[rna_count_df.index.isin(gene_lengths)]
    return filtered_df

def compute_tpm(rna_count_df, gene_lengths):
    for _, row in rna_count_df.iterrows():
        row["Counts"] /= gene_lengths[row.name]
    total_reads = rna_count_df["Counts"].sum()
    rna_count_df.loc[:,"TPM"] = (rna_count_df["Counts"] / total_reads) * 1e6

@click.command()
@click.option("--rna_10x", type=str, required=True)
@click.option("--gene_ref_bed", type=str, required=True)
@click.option("--output_file", type=str, default="rna_tpm.tsv")
def main(rna_10x, gene_ref_bed, output_file):
    gene_lengths = get_gene_lengths(gene_ref_bed)
    rna_count_df = get_gene_rna_counts(rna_10x)
    rna_count_df = filter_unknown_genes(rna_count_df, gene_lengths)

    # Calculate TPM
    compute_tpm(rna_count_df, gene_lengths)

    rna_count_df.to_csv(output_file, sep="\t")

if __name__ == "__main__":
    main()
