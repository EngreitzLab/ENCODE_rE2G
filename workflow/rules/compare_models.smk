
# compare cv-performance on training data across all models (note, this is not the true benchmarking performance CRISPR elements not overlapping prediction elements aren't considered)  
rule gather_model_performances:
	input:
		all_predictions = expand(os.path.join(MODELS_RESULTS_DIR, "{model}", "model", "training_predictions.tsv"), model=model_config["model"]),
		all_missing = expand(os.path.join(MODELS_RESULTS_DIR, "{model}", "merged_CRISPR_dataset.missing_from_features.NAfilled.tsv.gz"), model=model_config["model"]),
	output:
		comp_table = os.path.join(MODELS_RESULTS_DIR, "performance_across_models.tsv")
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
			--all_pred "{input.all_predictions}" \
			--all_missing "{input.all_missing}" \
			--model_config_file {params.model_config_file} \
			--output_file {output.comp_table}  \
			--crispr_data {params.crispr_dataset} \
			--out_dir {params.out_dir}
		"""

rule plot_model_performances:
	input:
		comp_table = os.path.join(MODELS_RESULTS_DIR, "performance_across_models.tsv")
	output:
		comp_plot_auprc = os.path.join(MODELS_RESULTS_DIR, "performance_across_models_auprc.pdf"),
		comp_plot_prec = os.path.join(MODELS_RESULTS_DIR, "performance_across_models_precision.pdf")
	conda:
		"../envs/encode_re2g.yml" 
	resources:
		mem_mb=ABC.determine_mem_mb
	script:
		"../scripts/model_training/plot_model_comparison.R"
