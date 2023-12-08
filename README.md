# ENCODE E2G
Generate [Encode E2G](https://github.com/karbalayghareh/ENCODE-E2G) input features based on [ABC](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction) output

Supports running ENCODE E2G with a pre-trained model. Put the pretrained model in the config, under `models_dir`

## Usage

Clone the repo

Load `ABC` submodule: `git submodule update --init`

Modify `config/config.yml` with your ABC biosamples config

```
mamba env create -f workflow/envs/encode_e2g.yml
conda activate encode_e2g
snakemake -j1 --use-conda
```

Output will show up in the `results/` directory
