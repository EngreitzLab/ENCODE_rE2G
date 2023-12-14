# ENCODE E2G
> :memo: **Note:** This repo is currently undergoing development. To access the version using for the encode_re2g paper, go to this [version](https://github.com/EngreitzLab/ENCODE_rE2G/tree/1906b6dcd97269374778e67592168c9da2dc455a). There are currently no clear instructions for stitching together the outputs from ABC, e2g features, and e2g, so use at your own discretion. We are working on creating 1 clean pipeline for the future

Generate [Encode E2G](https://github.com/karbalayghareh/ENCODE-E2G) input features based on [ABC](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction) output

Supports running ENCODE E2G with a pre-trained model. Put the pretrained model in the config, under `models_dir`

## Usage

Clone the repo and set it up for submdoule usage
```
git clone --recurse-submodules https://github.com/EngreitzLab/ENCODE_rE2G.git
git config --global submodule.recurse true
```
We use `ABC` as a submodule so this command will initialize it and sets up your git config to automatically keep the submodule up to date.

Modify `config/config.yml` with your ABC biosamples config

```
mamba env create -f workflow/envs/encode_e2g.yml
conda activate encode_e2g
snakemake -j1 --use-conda
```

Output will show up in the `results/` directory





## Supported Models

We have pre-trained ENCODE-rE2G on certain model types. You can find them in the `models` directory.
Each model must have the following:
1. model pickle file
2. feature table file
3. threshold file

The way we choose the model depends on the biosamples input. The code for model selection is in
the utils.smk file, under the `_get_biosample_model_dir` function.
