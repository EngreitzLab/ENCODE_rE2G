
# integrate default params with any overriding params for model training
rule generate_model_params:
	input:
	params:
		override_params = lambda wildcards: model_config.loc[wildcards.model, "override_params"],
		default_params = config["default_params"],
		scripts_dir = SCRIPTS_DIR
	output:
		final_params = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "model", "training_params.pkl")
	conda:
		"../envs/encode_re2g.yml" 
	shell:
		""" 
		python {params.scripts_dir}/model_training/get_params.py \
			--default_params "{params.default_params}" \
			--override_params "{params.override_params}" \
			--output_file {output.final_params} 
		"""

# generate trained model and cross-validated predictions on CRISPR data
rule train_model:
	input:
		crispr_features_processed = os.path.join(RESULTS_DIR, "{dataset}", "for_training.EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_NAfilled.tsv.gz"),
		feature_table = lambda wildcards: model_config.loc[wildcards.model, 'feature_table'],
		model_params = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "model", "training_params.pkl")
	params:
		epsilon = config["epsilon"],
		scripts_dir = SCRIPTS_DIR,
		polynomial = lambda wildcards: model_config.loc[wildcards.model, 'polynomial'],
		out_dir = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "model")
	output:
		trained_model = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "model", "model_full.pkl"),
		pred = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "model", "training_predictions.tsv"),
		#in_order = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "model", "training_data_in_order.tsv"),
		#shap = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "model", "shap_scores.tsv")
	conda:
		"../envs/encode_re2g.yml" 
	resources:
		mem_mb=64*1000
	shell: 
		""" 
		python {params.scripts_dir}/model_training/train_model.py \
			--crispr_features_file {input.crispr_features_processed} \
			--feature_table_file {input.feature_table} \
			--out_dir {params.out_dir} \
			--polynomial {params.polynomial} \
			--epsilon {params.epsilon} \
			--params_file {input.model_params}
		"""
