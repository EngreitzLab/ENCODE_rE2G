import re
import glob
re_pattern = re.compile("/oak/stanford/groups/engreitz/Users/agschwin/distal_regulation_paper/ENCODE_portal_submission/encode_e2g/results/dnase_only/thresholded_predictions/encode_e2g_predictions_(.*)_DNaseOnly_thresholded_predictions.tsv.gz")
files = glob.glob('/oak/stanford/groups/engreitz/Users/agschwin/distal_regulation_paper/ENCODE_portal_submission/encode_e2g/results/dnase_only/thresholded_predictions/encode_e2g_predictions*thresholded_predictions.tsv.gz')

BIOSAMPLES = [re_pattern.match(file).group(1) for file in files]

rule get_stats:
	input:
		thresholded = "/oak/stanford/groups/engreitz/Users/agschwin/distal_regulation_paper/ENCODE_portal_submission/encode_e2g/results/dnase_only/thresholded_predictions/encode_e2g_predictions_{biosample}_DNaseOnly_thresholded_predictions.tsv.gz"
	params:
		scripts_dir = SCRIPTS_DIR
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=determine_mem_mb
	output:
		stats = os.path.join(RESULTS_DIR, "{biosample}", "Metrics", "encode_e2g_predictions_threshold_stats.tsv")
	shell:
		"""
		python {params.scripts_dir}/get_stats.py --predictions {input.thresholded} --output_file {output.stats}
		"""

rule generate_plots:
	input:
		stat_files = [os.path.join(RESULTS_DIR, f"{biosample}", "Metrics", f"encode_e2g_predictions_threshold_stats.tsv") for biosample in BIOSAMPLES]
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
		python {params.scripts_dir}/generate_plots.py \
			--output_file {output.plots} \
			{input.stat_files}
		"""