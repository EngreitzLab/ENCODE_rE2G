CircleCI [![CircleCI](https://dl.circleci.com/status-badge/img/gh/EngreitzLab/ENCODE_rE2G.svg?style=svg)](https://app.circleci.com/pipelines/github/EngreitzLab/ENCODE_rE2G?branch=main)

# ENCODE-rE2G
> :memo: **Note:** This repo is currently undergoing development. To access the version using for the encode_re2g paper, go to this [version](https://github.com/EngreitzLab/ENCODE_rE2G/tree/1906b6dcd97269374778e67592168c9da2dc455a). There are currently no clear instructions for stitching together the outputs from ABC, e2g features, and e2g, so use at your own discretion. We are working on creating 1 clean pipeline for the future

ENCODE-rE2G is a logistic regression pipeline built on top of [ABC](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction). Given a chromatin accessibility input file, it will generate a list of enhancer-gene predictions. You can read the preprint paper [here](https://www.biorxiv.org/content/10.1101/2023.11.09.563812v1)

![image](https://github.com/EngreitzLab/ENCODE_rE2G/assets/10254642/ce6d33b5-2c5f-49cc-8a09-8142f7ac9b62)

## Set up

Clone the repo and set it up for submodule usage
```
git clone --recurse-submodules git@github.com:EngreitzLab/ENCODE_rE2G.git
git config --global submodule.recurse true
```
We use `ABC` as a submodule, so this command will initialize it and set up your git config to automatically keep the submodule up to date. Depending on download speeds this step can take some time, but should be done within 5-10 minutes.

## Apply a pretrained model
You'll need to use a certain model based on your input. (e.g DNase-seq or ATAC-seq? Do you have H3K27ac data?) We've pretrained all the models and determined the right thresholding to get E-G links at 70% recall of a [CRISPR-validated E-G links](https://github.com/EngreitzLab/CRISPR_comparison/tree/main). 

Modify the `ABC_BIOSAMPLES` field in `config/config.yaml` to point to your ABC config. Read more about ABC config [here](https://abc-enhancer-gene-prediction.readthedocs.io/en/latest/usage/getting_started.html#configuring-abc).
- We have not trained powerlaw models. If you don't have cell specific hic, opt to use megamap hic instead: https://s3.us-central-1.wasabisys.com/aiden-encode-hic-mirror/bifocals_iter2/tissues.hic. Megamap is a culmination of hic across many different tissues. 
- [Advanced] If applying a model that includes external features, you must define an `external_features_config` in your `biosamples_config.` See "Train model" section for details on this file.

Activate a conda environment that has [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) installed. Note that building conda environments might take some time depending on your system (typically less than 30 minutes).

```
mamba env create -f workflow/envs/encode_re2g.yml
conda activate encode_re2g
snakemake -j1 --use-conda
```

Based on your biosample config, we will find the right model to use for you. If we haven't trained that model before, an exception will get raised. 

Output will show up in the `results/` directory
- Binarized predictions will be located at `results/{biosample_name}/{model_name}/encode_e2g_predictions_threshold.{threshold}.tsv.gz`
- Non-thresholded models will be located at `results/{biosample_name}/{model_name}/encode_e2g_predictions.tsv.gz` with the score in a column named "ENCODE-rE2G.Score".

### Supported Models
We have pre-trained ENCODE-rE2G on certain model types. You can find them in the `models` directory.
Each model must have the following:
1. model pickle file (`model.pkl` corresponding to `model_full.pkl` from the model training workflow)
2. feature table file (`feature_table.tsv`, the corresponding feature table file from model training)
3. threshold file (`threshold_0.XXX` where predictions with a score greater than 0.XXX are binarized as true links.

The way we choose the model depends on the biosamples input. The code for model selection can be found [here](https://github.com/EngreitzLab/ENCODE_rE2G/blob/main/workflow/rules/utils.smk#L42).
 
 To override default model selection and specify a different model (either one you've trained yourself or the extended model), add a column called `model_dir` to your biosample config. Multiple model directories can be specified as a comma-separated list. NOTE: The genome-wide feature tables to reproduce the ENCODE-rE2G_Extended model included in the prediction files on Synapse.org for [K562](https://www.synapse.org/#!Synapse:syn59478344) and [GM12878](https://www.synapse.org/#!Synapse:syn59478343). To use these feature tables, download the feature tables and remove the ".Feature" suffix from feature name columns.

## Train model

**Important: Only train models for biosamples matching the corresponding CRISPR data (in this case, K562)**
- Much of the the model training code was adapted from Alireza Karbalayghareh's [original implementation](https://github.com/karbalayghareh/ENCODE-E2G).

Modify `config/config_training.yaml` with your model and dataset configs
- `model_config` has columns:  model, dataset, ABC_directory, feature_table, polynomial (do you want to use polynomial features?), and override_params (are there model training parameters you would like to change from the default logistic regression settings specfied in `config/config_training.yaml`?)
    - See [this example](https://pastebin.com/zt1868R3) `model_config` for how to specfiy override parameters. If there are no override_params, leave the column blank but still include the header.
    - Feature tables must be specified for each model (example: `resources/feature_tables`) with columns: feature (name in final table), input_col (name in ABC output), second_input (multiplied by input_col if provided), aggregate_function (how to combine feature values when a CRISPR element overlaps more than one ABC element), fill_value (how to replace NAs), nice_name (used when plotting)
    - Note that trained models generated using polynomial features cannot directly be used in the **Apply model** workflow
- `dataset_config` is an ABC biosamples config to generate ABC predictions for datasets without an existing ABC directory. 
- Each dataset must correspond to a unique ABC_directory, with "biosample" in `dataset_config` equals "dataset" in `model_config`. If no ABC_directory is indicated in `model_config`, it must have an entry in `dataset_config`.
- If you are including features in addition to those generated within the pipeline (e.g. a value in input_col or second_input of a feature table is not included in `reference_features` in `config/config_training/yaml`), you must also define how to add these values with an external_features_config, which you include in `dataset_config` in the optional column external_features_config:
    - An `external_features_config` has columns feature (corresponding to the missing input_col or second_input value), source_col (column name in the source file), aggregate_function (how to combine values when merging different element definitions), join_by, and source_file
    - join_by must be either "TargetGene" (feature is defined per gene) or "overlap" (feature is defined per element-gene pair)
    - If join_by is "TargetGene," source_file must be a .tsv with, at minimum, columns source_col and TargetGene. If join_by is "overlap," source_file must be a .tsv with, at minimum, columns chr, start, end, TargetGene, source_col.

Activate a conda environment that has [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) installed. 

```
mamba env create -f workflow/envs/encode_re2g.yml 
conda activate encode_re2g 
snakemake -s workflow/Snakefile_training -j1 --use-conda
```

### Output
`results/{biosample_name}/{model_name}/model_name/model_full.pkl`: full model trained on all chromosomes
`results/{biosample_name}/{model_name}/model/training_predictions.tsv`: rE2G predictions on CRISPR training data, using leave 1 chromosome out models
