# overlap activity-only feature table table with K562 CRISPR data
rule overlap_activity_only_features_crispr:
	input:
		features = os.path.join(RESULTS_DIR, "{biosample}/ActivityOnly_features.tsv.gz"),
		crispr = config['crispr_dataset'],
		feature_table_file = lambda wildcards: get_feature_table_file(wildcards.biosample),
		tss = config['gene_TSS500']
	output: 
		os.path.join(RESULTS_DIR, "{biosample}", "EPCrisprBenchmark_ensemble_data_GRCh38.K562_ActivityOnly_features_{nafill}.tsv.gz")
	params:
		filter_genes = "none",
		activity = lambda wildcards: BIOSAMPLE_ACTIVITES[wildcards.biosample]
	conda:
		"../envs/encode_re2g.yml" 
	resources:
		mem_mb=32*1000
	script:
		"../scripts/overlap_features_with_crispr_data.R"