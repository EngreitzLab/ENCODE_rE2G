def get_accessibility_files(wildcards):
	# Inputs have been validated so only DHS or ATAC is provided
	biosample = BIOSAMPLE_DF[BIOSAMPLE_DF["biosample"] == wildcards.biosample].iloc[0]
	files = biosample[biosample["default_accessibility_feature"]]
	return files.split(",")

rule get_stats:
	input:
		thresholded = os.path.join(RESULTS_DIR, "{biosample}", "{model_name}", "encode_e2g_predictions_threshold{threshold}.tsv.gz"),
		accessibility = get_accessibility_files
	params:
		scripts_dir = SCRIPTS_DIR
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=4*1000
	output:
		stats = os.path.join(RESULTS_DIR, "{biosample}", "{model_name}", "encode_e2g_predictions_threshold{threshold}_stats.tsv")
	shell:
		"""
		python {params.scripts_dir}/model_application/get_stats.py --predictions {input.thresholded} --accessibility "{input.accessibility}" --output_file {output.stats}
		"""

rule generate_plots:
	input:
		stat_files = expand(
            os.path.join(RESULTS_DIR, "{biosample}", "{model_name}", "encode_e2g_predictions_threshold{threshold}_stats.tsv"), zip, biosample=BIOSAMPLE_DF["biosample"],
            model_name=BIOSAMPLE_DF['model_dir_base'], threshold=BIOSAMPLE_DF["model_threshold"]
        )
	params:
		scripts_dir = SCRIPTS_DIR
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=determine_mem_mb
	output:
		plots = os.path.join(RESULTS_DIR, "qc_plots.pdf")
	shell:
		"""
		python {params.scripts_dir}/model_application/generate_plots.py \
			--output_file {output.plots} \
			{input.stat_files}
		"""
