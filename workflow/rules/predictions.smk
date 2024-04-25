
rule write_feature_table: # write feature table to output directory for future reference
    input:
        feature_table_file = lambda wildcards: get_feature_table_file(wildcards.biosample)
    output:
        feature_table_here = os.path.join(RESULTS_DIR, "{biosample}", "feature_table.tsv")
    run:
        shell("cat {input.feature_table_file} > {output.feature_table_here}")

# reformat external features config
rule format_external_features_config:
	input:
		dataset_config = config["ABC_BIOSAMPLES"]
	output:
		external_features_config = os.path.join(RESULTS_DIR, "{dataset}", "external_features_config.tsv"),
	params:
		e2g_path = config["E2G_DIR_PATH"]
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=8*1000
	script:
		"../scripts/format_external_features_config.R"

rule generate_e2g_predictions:
	input:
		final_features = os.path.join(RESULTS_DIR, "{biosample}", "final_features.tsv.gz"),
	params:
		epsilon = config["epsilon"],
		feature_table_file = lambda wildcards: get_feature_table_file(wildcards.biosample),
		trained_model = lambda wildcards: get_trained_model(wildcards.biosample),
		scripts_dir = SCRIPTS_DIR
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=determine_mem_mb
	output: 
		prediction_file = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", "encode_e2g_predictions.tsv.gz")
	shell: 
		""" 
		python {params.scripts_dir}/run_e2g.py \
			--predictions {input.final_features} \
			--feature_table_file {params.feature_table_file} \
			--epsilon {params.epsilon} \
			--trained_model {params.trained_model} \
			--output_file {output.prediction_file}
		"""

rule filter_e2g_predictions:
	input:
		prediction_file = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", "encode_e2g_predictions.tsv.gz")
	params:
		threshold = lambda wildcards: get_threshold(wildcards.biosample),
		include_self_promoter = config["include_self_promoter"],
		scripts_dir = SCRIPTS_DIR
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=determine_mem_mb
	output:
		thresholded = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", "encode_e2g_predictions_threshold{threshold}.tsv.gz")
	shell:
		"""
		python {params.scripts_dir}/threshold_e2g_predictions.py \
			--all_predictions_file {input.prediction_file} \
			--threshold {params.threshold} \
			--include_self_promoter {params.include_self_promoter} \
			--output_file {output.thresholded}
		"""