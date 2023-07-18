# ENCODE E2G Features
Generate [Encode E2G](https://github.com/karbalayghareh/ENCODE-E2G) input features based on [ABC](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction) output

## Usage

Modify `config/config.yml` to point to your ABC prediction file

```
mamba env create -f workflow/envs/encode_e2g_features.yml
conda activate encode_e2g_features
snakemake -j1
```

Output will show up in the `results/` directory
