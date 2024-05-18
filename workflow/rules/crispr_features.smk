rule make_dataset_feature_table:
	input:
		model_config = config["model_config"],
	output:
		dataset_features = os.path.join(RESULTS_DIR, "{dataset}", "feature_table.tsv")
	params:
		model_config = config["model_config"],
		e2g_path = config["E2G_DIR_PATH"]
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=8*1000
	script:
		"../scripts/feature_tables/combine_feature_tables.R"

# reformat external features config
rule format_external_features_config:
	input:
		dataset_config = config["dataset_config"]
	output:
		external_features_config = os.path.join(RESULTS_DIR, "{dataset}", "external_features_config.tsv"),
	params:
		e2g_path = config["E2G_DIR_PATH"]
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=8*1000
	script:
		"../scripts/feature_tables/format_external_features_config.R"

# overlap feature table  with K562 CRISPR data
rule overlap_features_crispr:
	input:
		features = os.path.join(RESULTS_DIR, "{dataset}", "genomewide_features.tsv.gz"),
		crispr = config['crispr_dataset'],
		feature_table_file = os.path.join(RESULTS_DIR, "{dataset}", "feature_table.tsv"),
		tss = config['gene_TSS500']
	output: 
		features = os.path.join(RESULTS_DIR, "{dataset}", "EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_{nafill}.tsv.gz"),
		missing = os.path.join(RESULTS_DIR, "{dataset}",  "missing.EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_{nafill}.tsv.gz")
	conda:
		"../envs/encode_re2g.yml" 
	resources:
		mem_mb=32*1000
	script:
		"../scripts/feature_tables/overlap_features_with_crispr_data.R"

# process data for model training: rename columns, apply filter features, filter to gene list
# note: we use the NAfilled CRISPR feature data here!
rule process_crispr_data:
	input:
		crispr_features = os.path.join(RESULTS_DIR, "{dataset}", "EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_{nafill}.tsv.gz")
	params:
		genes = config["gene_TSS500"]
	output:
		processed = os.path.join(RESULTS_DIR, "{dataset}",  "for_training.EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_{nafill}.tsv.gz")
	conda:
		"../envs/encode_re2g.yml" 
	resources:
		mem_mb=32*1000
	script:
		"../scripts/feature_tables/process_crispr_data.R"
