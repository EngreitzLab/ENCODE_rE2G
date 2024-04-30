
# compare cv-performance on training data across all models (note, this is not the true benchmarking performance CRISPR elements not overlapping prediction elements aren't considered)  
rule gather_model_performances:
	input:
		all_predictions = expand(os.path.join(RESULTS_DIR, "{dataset}", "{model}", "model", "training_predictions.tsv"), zip, dataset=model_config["dataset"], model=model_config["model"])
	output:
		comp_table = os.path.join(RESULTS_DIR, "performance_across_models.tsv")
	params:
		scripts_dir = SCRIPTS_DIR,
		out_dir = RESULTS_DIR,
		model_config_file = config["model_config"],
		crispr_dataset = config["crispr_dataset"]
	conda:
		"../envs/encode_re2g.yml" 
	resources:
		mem_mb=64*1000
	shell: 
		""" 
		python {params.scripts_dir}/model_training/compare_all_models.py \
			--model_config_file {params.model_config_file} \
			--output_file {output.comp_table}  \
			--crispr_data {params.crispr_dataset} \
			--out_dir {params.out_dir}
		"""

rule plot_model_performances:
	input:
		comp_table = os.path.join(RESULTS_DIR, "performance_across_models.tsv")
	output:
		comp_plot_auprc = os.path.join(RESULTS_DIR, "performance_across_models_auprc.pdf"),
		comp_plot_prec = os.path.join(RESULTS_DIR, "performance_across_models_precision.pdf")
	conda:
		"../envs/encode_re2g.yml" 
	resources:
		mem_mb=determine_mem_mb
	script:
		"../scripts/model_training/plot_model_comparison.R"
