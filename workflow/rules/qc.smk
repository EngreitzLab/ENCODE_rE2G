rule get_stats:
	input:
		thresholded = os.path.join(RESULTS_DIR, "{biosample}", "encode_e2g_predictions_threshold{threshold}.tsv.gz")
	params:
		scripts_dir = SCRIPTS_DIR
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=4*1000
	output:
		stats = os.path.join(RESULTS_DIR, "{biosample}", "encode_e2g_predictions_threshold{threshold}_stats.tsv")
	shell:
		"""
		python {params.scripts_dir}/get_stats.py --predictions {input.thresholded} --output_file {output.stats}
		"""

rule generate_plots:
	input:
		expand(
			os.path.join(RESULTS_DIR, "{biosample}", f"encode_e2g_predictions_threshold{config['threshold']}_stats.tsv"), biosample=BIOSAMPLES
		)
	params:
		scripts_dir = SCRIPTS_DIR
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=4*1000
	output:
		plots = os.path.join(RESULTS_DIR, "qc_plots_threshold{threshold}.pdf")
	shell:
		"""
		python {params.scripts_dir}/generate_plots.py \
			--results_dir {RESULTS_DIR} \
			--output_file {output.plots} \
			# --y2ave_metadata /oak/stanford/groups/engreitz/Users/atan5133/igvf_dataset_processing/Y2AVE_SingleCellDatasets.CellClusterTable.tsv
		"""