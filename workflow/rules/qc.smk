rule get_stats:
	input:
		thresholded = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", "encode_e2g_predictions_threshold{threshold}.tsv.gz")
	params:
		scripts_dir = SCRIPTS_DIR
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=4*1000
	output:
		stats = os.path.join(RESULTS_DIR, "{biosample}", "Metrics", "encode_e2g_predictions_threshold{threshold}_stats.tsv")
	shell:
		"""
		python {params.scripts_dir}/get_stats.py --predictions {input.thresholded} --output_file {output.stats}
		"""

rule generate_plots:
	input:
		[os.path.join(RESULTS_DIR, f"{biosample}", "Predictions", f"encode_e2g_predictions_threshold{get_threshold(biosample)}_stats.tsv") for biosample in BIOSAMPLES]
	params:
		scripts_dir = SCRIPTS_DIR
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=4*1000
	output:
		plots = os.path.join(RESULTS_DIR, "qc_plots.pdf")
	shell:
		"""
		python {params.scripts_dir}/generate_plots.py \
			--results_dir {RESULTS_DIR} \
			--output_file {output.plots} \
			# --y2ave_metadata /oak/stanford/groups/engreitz/Users/atan5133/igvf_dataset_processing/Y2AVE_SingleCellDatasets.CellClusterTable.tsv
		"""