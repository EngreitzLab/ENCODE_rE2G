import os
from pathlib import Path
from typing import List

import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

STATS_SUFFIX = "_stats.tsv"
VALUE_KEY = "Value"
METRICS = [
    "num_sequencing_reads",
    "num_enh",
    "num_genes",
    "num_enh_gene_links",
    "num_genes_with_1_enh_min",
    "mean_num_genes_per_enh",
    "mean_num_enh_per_gene",
    "mean_num_enh_per_gene_no_prom",
    "mean_log10_dist_to_tss",
    "mean_enh_region_size",
]


def load_stat_files(stat_files) -> List[pd.DataFrame]:
    results = {}
    for stat_file in stat_files:
        cell_cluster = Path(stat_file).parts[-3]
        results[cell_cluster] = pd.read_csv(stat_file, sep="\t").set_index("Metric")
    return results


def snake_to_title(snake_case):
    words = snake_case.split("_")
    # Capitalize the first letter of each word except the first one
    title = " ".join(word.capitalize() for word in words)
    return title


def plot_violin(stat_dfs, metric):
    plt.clf()
    points = [df.loc[metric, VALUE_KEY] for df in stat_dfs.values()]
    mean, median = np.mean(points), np.median(points)
    ax = sns.violinplot(y=points)
    ax.set_title(snake_to_title(metric))
    plt.scatter([], [], label=f"n={len(points)}\nMean={mean:.2f}\nMedian={median:.2f}")
    plt.legend()
    return ax.get_figure()


def plot_swarm_box(stat_dfs, metric):
    plt.clf()
    points = [df.loc[metric, VALUE_KEY] for df in stat_dfs.values()]
    mean, median = np.mean(points), np.median(points)
    ax = sns.swarmplot(y=points)
    ax = sns.boxplot(y=points)
    ax.set_title(snake_to_title(metric))
    plt.scatter([], [], label=f"n={len(points)}\nMean={mean:.2f}\nMedian={median:.2f}")
    plt.legend()
    return ax.get_figure()


def plot_scatter(stat_dfs, metadata_df, metric):
    plt.clf()
    X, Y = [], []
    K562_X, K562_Y = [], []
    for cell_cluster in stat_dfs:
        row = metadata_df[metadata_df["CellClusterID"] == cell_cluster].iloc[0]
        num_fragments = row["nCells"] * row["MeanATACFragmentsPerCell"]
        metric_val = stat_dfs[cell_cluster].loc[metric, VALUE_KEY]
        if "K562" in row["ManualAnnotationLabel"]:
            K562_X.append(int(num_fragments))
            K562_Y.append(metric_val)
        else:
            X.append(int(num_fragments))
            Y.append(metric_val)

    ax = sns.scatterplot(x=X, y=Y, label=f"n={len(X)}")
    ax = sns.scatterplot(x=K562_X, y=K562_Y, color="red", label=f"K562 n={len(K562_X)}")
    plt.xscale("log")
    plt.axvline(x=2e6, color="red", linestyle="-", label=f"2 million fragments")
    ax.set_xlabel("Num Fragments")
    ax.set_ylabel(snake_to_title(metric))
    ax.set_title(snake_to_title(metric))
    plt.legend()

    return ax.get_figure()


def add_intro_page(pdf_writer):
    plt.clf()
    plt.axis("off")
    plt.text(
        0.5,
        1,
        "QC Plots for ENCODE_rE2G\n",
        ha="center",
        va="center",
        fontsize="large",
        wrap=True,
    )

    plt.text(
        0.5,
        0.5,
        "N represents number of datasets\n\n"
        "The metrics represent E-G pairs after applying a threshold\n\n"
        "Metric definitions are found here: https://github.com/EngreitzLab/ENCODE_rE2G/blob/main/workflow/scripts/get_stats.py",
        ha="center",
        va="center",
        fontsize="medium",
        wrap=True,
    )
    pdf_writer.savefig()


def save_outlier_stats(stat_dfs, metric, pdf_writer):
    cell_cluster_values = {}
    for cell_cluster, stat_df in stat_dfs.items():
        cell_cluster_values[cell_cluster] = stat_df.loc[metric, VALUE_KEY]

    sorted_items = sorted(cell_cluster_values.items(), key=lambda item: item[1])

    top_5 = [f"{item[0]}: {item[1]}" for item in list(reversed(sorted_items))[:5]]
    bottom_5 = [f"{item[0]}: {item[1]}" for item in sorted_items[:5]]

    plt.clf()
    plt.axis("off")
    plt.text(
        0.5,
        1,
        f"{metric}",
        ha="center",
        va="center",
        fontsize="large",
        wrap=True,
    )
    nl = "\n"
    plt.text(
        0.5,
        0.5,
        f"Top 5 datasets:{nl}"
        f"{nl.join(top_5)}"
        f"{nl*3}"
        f"Bottom 5 datasets:{nl}"
        f"{nl.join(bottom_5)}"
        f"{nl*3}",
        ha="center",
        va="center",
        fontsize="medium",
        wrap=True,
    )
    pdf_writer.savefig()


@click.command()
@click.option("--output_file", type=str, default="qc_plots.pdf")
@click.option("--y2ave_metadata", type=str)
@click.argument("stat_files", nargs=-1, type=click.Path(exists=True))
def main(output_file, stat_files, y2ave_metadata):
    if y2ave_metadata:
        metadata_df = pd.read_csv(y2ave_metadata, sep="\t")
    stat_dfs = load_stat_files(stat_files)
    with PdfPages(output_file) as pdf_writer:
        add_intro_page(pdf_writer)
        for metric in METRICS:
            pdf_writer.savefig(plot_violin(stat_dfs, metric))
            pdf_writer.savefig(plot_swarm_box(stat_dfs, metric))
            save_outlier_stats(stat_dfs, metric, pdf_writer)
            if y2ave_metadata:
                pdf_writer.savefig(plot_scatter(stat_dfs, metadata_df, metric))


if __name__ == "__main__":
    main()
