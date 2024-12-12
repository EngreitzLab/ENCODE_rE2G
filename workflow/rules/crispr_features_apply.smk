# merge features with crispr data
rule overlap_features_crispr_apply:
    input:
        prediction_file = os.path.join(RESULTS_DIR, "{biosample}", "{model_name}", "encode_e2g_predictions.tsv.gz"),
        crispr = config['crispr_dataset'],
        feature_table_file = os.path.join(RESULTS_DIR, "{biosample}", "feature_table.tsv"),
        tss = config['gene_TSS500'],
    params:
        scripts_dir = SCRIPTS_DIR
    output: 
        features = os.path.join(RESULTS_DIR, "{biosample}", "{model_name}", "EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_{nafill}.tsv.gz"),
    conda:
        "../envs/encode_re2g.yml" 
    resources:
        mem_mb=32*1000
    script:
        "../scripts/model_application/merge_features_with_crispr_data_apply.R"

# calculate performance metrics 
rule crispr_benchmarking:
    input:
        crispr_features = expand(os.path.join(RESULTS_DIR, "{biosample}", "{model_name}", "EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_NAfilled.tsv.gz"), zip, biosample=BIOSAMPLE_DF["biosample"], model_name=BIOSAMPLE_DF["model_dir_base"]),
    output:
        comp_table = os.path.join(RESULTS_DIR, "crispr_benchmarking_performance_summary.tsv")
    params:
        model_names = BIOSAMPLE_DF["model_dir_base"].tolist(),
        model_thresh = BIOSAMPLE_DF["model_threshold"].tolist(),
        scripts_dir = SCRIPTS_DIR
    conda:
        "../envs/encode_re2g.yml" 
    resources:
        mem_mb=64*1000
    shell: 
        """ 
        python {params.scripts_dir}/model_training/benchmark_performance.py \
            --crispr_features "{input.crispr_features}" \
            --output_file {output.comp_table} \
            --model_thresholds "{params.model_thresh}" \
            --model_names "{params.model_names}" \
        """
