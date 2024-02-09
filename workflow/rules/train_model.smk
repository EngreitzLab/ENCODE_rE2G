
# add activity only features and retain columns from ABC for future features 
rule activity_only_features:
	input:
		feature_table_file = lambda wildcards: model_config.loc[wildcards.model, 'feature_table'],
		abc = lambda wildcards: os.path.join(ABC_BIOSAMPLES_DIR[wildcards.dataset], "Predictions", "EnhancerPredictionsAllPutative.tsv.gz"),
		NumCandidateEnhGene = os.path.join(RESULTS_DIR, "{dataset}", "NumCandidateEnhGene.tsv"),
		NumTSSEnhGene = os.path.join(RESULTS_DIR, "{dataset}", "NumTSSEnhGene.tsv"),
		NumEnhancersEG5kb = os.path.join(RESULTS_DIR, "{dataset}", "NumEnhancersEG5kb.txt"),
		SumEnhancersEG5kb = os.path.join(RESULTS_DIR, "{dataset}", "SumEnhancersEG5kb.txt"),
		ubiqExprGenes = config["ubiq_expr_genes"]
	output: 
		predictions_extended = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "ActivityOnly_features.tsv.gz")
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=32*1000
	script:
		"../scripts/activity_only_features.R"

# add external features (modify this rule in future pipelines- note, outputs should have name INPUT COL and will be renamed later)
rule add_external_features:
	input:
		predictions_extended = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "ActivityOnly_features.tsv.gz")
	output:
		plus_external_features = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "ActivityOnly_plus_external_features.tsv.gz")
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=4*1000
	shell:
		"""
		cp {input.predictions_extended} {output.plus_external_features}
		"""

# compute interaction or squared terms, fill NAs, rename features to finals, fill nas
rule gen_final_features:
	input:
		plus_external_features = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "ActivityOnly_plus_external_features.tsv.gz"),
		feature_table_file = lambda wildcards: model_config.loc[wildcards.model, 'feature_table']
	output:
		final_features = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "final_features.tsv.gz")
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=32*1000
	script:
		"../scripts/gen_final_features.R"
	
# overlap activity-only feature table table with K562 CRISPR data
rule overlap_features_crispr:
	input:
		features = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "final_features.tsv.gz"),
		crispr = config['crispr_dataset'],
		feature_table_file = lambda wildcards: model_config.loc[wildcards.model, 'feature_table'],
		tss = config['gene_TSS500']
	output: 
		os.path.join(RESULTS_DIR, "{dataset}", "{model}", "EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_{nafill}.tsv.gz")
	params:
		filter_genes = "none",
	conda:
		"../envs/encode_re2g.yml" 
	resources:
		mem_mb=32*1000
	script:
		"../scripts/overlap_features_with_crispr_data.R"

# process data for model training: rename columns, apply filter features, filter to gene list
# note: we use the NAfilled CRISPR feature data here!
rule process_crispr_data:
	input:
		crispr_features = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_{nafill}.tsv.gz")
	params:
		genes = config["gene_TSS500"]
	output:
		processed = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "for_training.EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_{nafill}.tsv.gz")
	conda:
		"../envs/encode_re2g.yml" 
	resources:
		mem_mb=32*1000
	script:
		"../scripts/process_crispr_data.R"

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
		python {params.scripts_dir}/get_params.py \
			--default_params "{params.default_params}" \
			--override_params "{params.override_params}" \
			--output_file {output.final_params} 
		"""

# generate trained model and cross-validated predictions on CRISPR data
rule train_model:
	input:
		crispr_features_processed = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "for_training.EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_NAfilled.tsv.gz"),
		feature_table = lambda wildcards: model_config.loc[wildcards.model, 'feature_table'],
		model_params = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "model", "training_params.pkl")
	params:
		epsilon = config["epsilon"],
		scripts_dir = SCRIPTS_DIR,
		polynomial = lambda wildcards: model_config.loc[wildcards.model, 'polynomial'],
		out_dir = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "model")
	output:
		trained_model = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "model", "model_full.pkl"),
		pred = os.path.join(RESULTS_DIR, "{dataset}", "{model}", "model", "training_predictions.tsv")
	conda:
		"../envs/encode_re2g.yml" 
	resources:
		mem_mb=64*1000
	shell: 
		""" 
		python {params.scripts_dir}/train_model.py \
			--crispr_features_file {input.crispr_features_processed} \
			--feature_table_file {input.feature_table} \
			--out_dir {params.out_dir} \
			--polynomial {params.polynomial} \
			--epsilon {params.epsilon} \
			--params_file {input.model_params}
		"""
		
# compare cv-performance on training data across all models (note, this is not the true benchmarking performance CRISPR elements not overlapping prediction elements aren't considered)  
rule gather_model_performances:
	input:
		all_predictions = expand(os.path.join(RESULTS_DIR, "{dataset}", "{model}", "model", "training_predictions.tsv"), zip, dataset=model_config["dataset"], model=model_config["model"])
	output:
		comp_table = os.path.join(RESULTS_DIR, "relative_performance_across_models.tsv")
	params:
		scripts_dir = SCRIPTS_DIR,
		out_dir = RESULTS_DIR,
		model_config_file = config["model_config"]
	conda:
		"../envs/encode_re2g.yml" 
	resources:
		mem_mb=64*1000
	shell: 
		""" 
		python {params.scripts_dir}/compare_all_models.py \
			--model_config_file {params.model_config_file} \
			--output_file {output.comp_table}  \
			--out_dir {params.out_dir}
		"""

rule plot_model_performances:
	input:
		comp_table = os.path.join(RESULTS_DIR, "relative_performance_across_models.tsv")
	output:
		comp_plot = os.path.join(RESULTS_DIR, "relative_performance_across_models.pdf")
	conda:
		"../envs/encode_re2g.yml" 
	resources:
		mem_mb=64*1000
	script:
		"../scripts/plot_model_comparison.R"