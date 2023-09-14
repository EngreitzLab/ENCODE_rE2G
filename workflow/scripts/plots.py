import click
import glob
import os
import pandas as pd
from typing import List
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

STATS_SUFFIX = "_stats.tsv"
VALUE_KEY = "Value"
METRICS = [
    "num_enh",
    "num_genes",
    "num_enh_gene_links",
    "num_genes_with_1_enh_min",
    "mean_num_genes_per_enh",
    "mean_num_enh_per_gene",
    "mean_log10_dist_to_tss",
    "mean_enh_region_size"
]

def load_stat_files(stat_files) -> List[pd.DataFrame]:
    results = []
    for f in stat_files:
        results.append(pd.read_csv(f, sep="\t").set_index("Metric"))
    return results

def snake_to_title(snake_case):
    words = snake_case.split('_')
    # Capitalize the first letter of each word except the first one
    title = ' '.join(word.capitalize() for word in words)
    return title

def plot_violin(stat_dfs, metric):
    plt.clf()
    points = [df.loc[metric, VALUE_KEY] for df in stat_dfs]
    mean, median = np.mean(points), np.median(points)
    ax = sns.violinplot(y=points)
    ax.set_title(snake_to_title(metric))
    plt.scatter([], [], label=f"n={len(points)}\nMean={mean:.2f}\nMedian={median:.2f}")
    plt.legend()
    return ax.get_figure()
    

@click.command()
@click.option("--results_dir", type=str, required=True)
@click.option("--output_file", type=str, default="qc_plots.pdf")
def main(results_dir, output_file):
    stat_files = glob.glob(os.path.join(results_dir, '*', f"*{STATS_SUFFIX}"))
    stat_dfs = load_stat_files(stat_files)
    with PdfPages(output_file) as pdf_writer:
        for metric in METRICS:
            pdf_writer.savefig(plot_violin(stat_dfs, metric))

if __name__ == "__main__":
    main()