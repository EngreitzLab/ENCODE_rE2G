# forward sequential feature selection
rule calculate_forward_feature_selection:
	input:
		crispr_features_processed = os.path.join(RESULTS_DIR, "{dataset}", "for_training.EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_NAfilled.tsv.gz"),
		feature_table = lambda wildcards: model_config.loc[wildcards.model, 'feature_table'],
		model_params = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "model", "training_params.pkl")
	params:
		epsilon = config["epsilon"],
		scripts_dir = SCRIPTS_DIR,
		polynomial = lambda wildcards: model_config.loc[wildcards.model, 'polynomial'],
		out_dir = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "feature_analysis")
	output:
		results = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "feature_analysis", "forward_feature_selection.tsv"),
	conda:
		"../envs/encode_re2g.yml" 
	resources:
		mem_mb=determine_mem_mb
	shell: 
		""" 
		python {params.scripts_dir}/model_training/forward_sequential_feature_selection.py \
			--crispr_features_file {input.crispr_features_processed} \
			--feature_table_file {input.feature_table} \
			--out_dir {params.out_dir} \
			--polynomial {params.polynomial} \
			--epsilon {params.epsilon} \
			--params_file {input.model_params}
		"""

rule plot_forward_feature_selection:
	input:
		results = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "feature_analysis", "forward_feature_selection.tsv"),
		feature_table = lambda wildcards: model_config.loc[wildcards.model, 'feature_table'],
	params:
		polynomial = lambda wildcards: model_config.loc[wildcards.model, 'polynomial']
	output:
		out_auprc = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "feature_analysis", "forward_feature_selection_auprc.pdf"),
		out_prec = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "feature_analysis", "forward_feature_selection_precision.pdf")
	conda:
		"../envs/encode_re2g.yml" 
	resources:
		mem_mb=determine_mem_mb
	script:
		"../scripts/model_training/plot_sffs.R"

# backward sequential feature selection
rule calculate_backward_feature_selection:
	input:
		crispr_features_processed = os.path.join(RESULTS_DIR, "{dataset}", "for_training.EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_NAfilled.tsv.gz"),
		feature_table = lambda wildcards: model_config.loc[wildcards.model, 'feature_table'],
		model_params = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "model", "training_params.pkl")
	params:
		epsilon = config["epsilon"],
		scripts_dir = SCRIPTS_DIR,
		polynomial = lambda wildcards: model_config.loc[wildcards.model, 'polynomial'],
		out_dir = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "feature_analysis")
	output:
		results = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "feature_analysis", "backward_feature_selection.tsv"),
	conda:
		"../envs/encode_re2g.yml" 
	resources:
		mem_mb=determine_mem_mb
	shell: 
		""" 
		python {params.scripts_dir}/model_training/backward_sequential_feature_selection.py \
			--crispr_features_file {input.crispr_features_processed} \
			--feature_table_file {input.feature_table} \
			--out_dir {params.out_dir} \
			--polynomial {params.polynomial} \
			--epsilon {params.epsilon} \
			--params_file {input.model_params}
		"""

rule plot_backward_feature_selection:
	input:
		results = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "feature_analysis", "backward_feature_selection.tsv"),
		feature_table = lambda wildcards: model_config.loc[wildcards.model, 'feature_table']
	params:
		polynomial = lambda wildcards: model_config.loc[wildcards.model, 'polynomial']
	output:
		out_auprc = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "feature_analysis", "backward_feature_selection_auprc.pdf"),
		out_prec = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "feature_analysis", "backward_feature_selection_precision.pdf")
	conda:
		"../envs/encode_re2g.yml" 
	resources:
		mem_mb=determine_mem_mb
	script:
		"../scripts/model_training/plot_sbfs.R"

# compare all features sets
rule compare_all_feature_sets:
	input:
		crispr_features_processed = os.path.join(RESULTS_DIR, "{dataset}", "for_training.EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_NAfilled.tsv.gz"),
		feature_table = lambda wildcards: model_config.loc[wildcards.model, 'feature_table'],
		model_params = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "model", "training_params.pkl")
	params:
		epsilon = config["epsilon"],
		scripts_dir = SCRIPTS_DIR,
		out_dir = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "feature_analysis")
	output:
		results = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "feature_analysis", "all_feature_sets.tsv"),
	conda:
		"../envs/encode_re2g.yml" 
	resources:
		mem_mb=determine_mem_mb,
		runtime='24h'
	shell: 
		""" 
		python {params.scripts_dir}/model_training/compare_all_feature_sets.py \
			--crispr_features_file {input.crispr_features_processed} \
			--feature_table_file {input.feature_table} \
			--out_dir {params.out_dir} \
			--epsilon {params.epsilon} \
			--params_file {input.model_params}
		"""

# permuation feature importance
rule calculate_permutation_feature_importance:
	input:
		crispr_features_processed = os.path.join(RESULTS_DIR, "{dataset}", "for_training.EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_NAfilled.tsv.gz"),
		feature_table = lambda wildcards: model_config.loc[wildcards.model, 'feature_table'],
		model_params = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "model", "training_params.pkl")
	params:
		epsilon = config["epsilon"],
		n_repeats = 20,
		scripts_dir = SCRIPTS_DIR,
		polynomial = lambda wildcards: model_config.loc[wildcards.model, 'polynomial'],
		out_dir = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "feature_analysis")
	output:
		results = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "feature_analysis", "permutation_feature_importance.tsv"),
	conda:
		"../envs/encode_re2g.yml" 
	resources:
		mem_mb=determine_mem_mb
	shell: 
		""" 
		python {params.scripts_dir}/model_training/permutation_feature_importance.py \
			--crispr_features_file {input.crispr_features_processed} \
			--feature_table_file {input.feature_table} \
			--out_dir {params.out_dir} \
			--polynomial {params.polynomial} \
			--epsilon {params.epsilon} \
			--n_repeats {params.n_repeats} \
			--params_file {input.model_params}
		"""

rule plot_permutation_feature_importance:
	input:
		results = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "feature_analysis", "permutation_feature_importance.tsv"),
		feature_table = lambda wildcards: model_config.loc[wildcards.model, 'feature_table']
	params:
		polynomial = lambda wildcards: model_config.loc[wildcards.model, 'polynomial'],
		n_repeats = 20
	output:
		out_auprc = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "feature_analysis", "permutation_feature_importance_auprc.pdf"),
		out_prec = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "feature_analysis", "permutation_feature_importance_precision.pdf")
	conda:
		"../envs/encode_re2g.yml" 
	resources:
		mem_mb=determine_mem_mb
	script:
		"../scripts/model_training/plot_pfi.R"
