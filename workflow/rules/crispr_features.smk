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

# overlap feature table  with CRISPR data for this dataset
rule overlap_features_crispr_for_dataset:
	input:
		features = os.path.join(RESULTS_DIR, "{dataset}", "genomewide_features.tsv.gz"),
		feature_table_file = os.path.join(RESULTS_DIR, "{dataset}", "feature_table.tsv"),
		crispr = lambda wildcards: config['crispr_dataset'][wildcards.crispr_dataset],
		tss = config['gene_TSS500']
	params:
		crispr_cell_type = lambda wildcards: dataset_config.loc[wildcards.dataset, "crispr_cell_type"]
	output: 
		features = os.path.join(RESULTS_DIR, "{dataset}", "CRISPR_dataset_{crispr_dataset}.overlapping_features.{nafill}.tsv.gz"),
		missing = os.path.join(RESULTS_DIR, "{dataset}", "CRISPR_dataset_{crispr_dataset}.missing_from_features.{nafill}.tsv.gz")
	conda:
		"../envs/encode_re2g.yml" 
	resources:
		mem_mb=32*1000
	script:
		"../scripts/feature_tables/overlap_features_with_crispr_data.R"

# process data for model training: rename columns, apply filter features, filter to gene list, combine across datasets
# note: we use the NAfilled CRISPR feature data here!
def get_crispr_files_for_model(wildcards, file_type):
	datasets = model_dataset_dict[wildcards.model].values()
	crispr_dataset = model_config.loc[wildcards.model, "crispr_dataset"]
	files = [os.path.join(RESULTS_DIR, ds, f"CRISPR_dataset_{crispr_dataset}.{file_type}.{wildcards.nafill}.tsv.gz") for ds in datasets]
	return files

rule process_crispr_data:
	input:
		crispr_features = lambda wildcards: get_crispr_files_for_model(wildcards, "overlapping_features"),
		crispr_missing = lambda wildcards: get_crispr_files_for_model(wildcards, "missing_from_features"),
	params:
		genes = config["gene_TSS500"],
	output:
		processed = os.path.join(MODELS_RESULTS_DIR, "{model}",  "for_training.merged_CRISPR_dataset.overlapping_features.{nafill}.tsv.gz"),
		missing = os.path.join(MODELS_RESULTS_DIR, "{model}",  "merged_CRISPR_dataset.missing_from_features.{nafill}.tsv.gz")
	conda:
		"../envs/encode_re2g.yml" 
	resources:
		mem_mb=32*1000
	script:
		"../scripts/feature_tables/process_crispr_data.R"
