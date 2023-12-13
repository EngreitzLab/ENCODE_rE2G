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
		mem_mb=32*1000
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
		mem_mb=32*1000
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