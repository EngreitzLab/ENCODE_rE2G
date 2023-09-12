import os
import shutil

import click
import pandas as pd


def get_num_enh_gene_links(df):
    return

def get_genes_with_1_enh_min(df):
    return

def get_num_genes_per_enh(df):
    return

def get_num_enh_per_gene(df):
    return

def get_mean_dist_to_tss(df):
    return

def get_mean_enh_region_size(df):
    return

@click.command()
@click.option("--predictions")
def main(predictions):
    df = pd.read_csv(predictions, sep="\t")
    # column is metric, value

if __name__ == "__main__":
    main()
