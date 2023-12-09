rule generate_e2g_predictions:
	input:
		predictions_extended = os.path.join(RESULTS_DIR, "{biosample}", "ActivityOnly_features.tsv.gz"),
	params:
		feature_table_file = lambda wildcards: config["feature_table"][BIOSAMPLE_ACTIVITES[wildcards.biosample]],
		epsilon = config["epsilon"],
		models_dir = lambda wildcards: config["models_dir"][BIOSAMPLE_ACTIVITES[wildcards.biosample]],
		scripts_dir = SCRIPTS_DIR
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=32*1000
	output: 
		prediction_file = os.path.join(RESULTS_DIR, "{biosample}", "encode_e2g_predictions.tsv.gz")
	shell: 
		""" 
		python {params.scripts_dir}/run_e2g.py \
			--predictions {input.predictions_extended} \
			--feature_table_file {params.feature_table_file} \
			--epsilon {params.epsilon} \
			--models_dir {params.models_dir} \
			--output_file {output.prediction_file}
		"""

rule filter_e2g_predictions:
	input:
		prediction_file = os.path.join(RESULTS_DIR, "{biosample}", "encode_e2g_predictions.tsv.gz")
	params:
		threshold = config["threshold"],
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=4*1000
	output:
		thresholded = os.path.join(RESULTS_DIR, "{biosample}", "encode_e2g_predictions_threshold{threshold}.tsv.gz")
	shell:
		"""
		zcat {input.prediction_file} | awk -F'\t' '$NF >= {params.threshold}' | gzip > {output.thresholded}
		"""