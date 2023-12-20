# create activity-only feature table
rule activity_only_features:
	input:
		feature_table_file = lambda wildcards: get_feature_table_file(wildcards.biosample),
		abc = lambda wildcards: os.path.join(ABC_dirs[wildcards.biosample], "Predictions", "EnhancerPredictionsAllPutative.tsv.gz"),
		NumCandidateEnhGene = os.path.join(RESULTS_DIR, "{biosample}", "NumCandidateEnhGene.tsv"),
		NumTSSEnhGene = os.path.join(RESULTS_DIR, "{biosample}", "NumTSSEnhGene.tsv"),
		NumEnhancersEG5kb = os.path.join(RESULTS_DIR, "{biosample}", "NumEnhancersEG5kb.txt"),
		SumEnhancersEG5kb = os.path.join(RESULTS_DIR, "{biosample}", "SumEnhancersEG5kb.txt"),
		ubiqExprGenes = config["ubiq_expr_genes"]
	output: 
		predictions_extended = os.path.join(RESULTS_DIR, "{biosample}", "ActivityOnly_features.tsv.gz")
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=determine_mem_mb
	script:
		"../scripts/activity_only_features.R"

# add external features (modify this rule in future pipelines- note, outputs should have name INPUT COL and will be renamed later)
rule add_external_features:
	input:
		predictions_extended = os.path.join(RESULTS_DIR, "{biosample}", "ActivityOnly_features.tsv.gz")
	output:
		plus_external_features = os.path.join(RESULTS_DIR, "{biosample}", "ActivityOnly_plus_external_features.tsv.gz")
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=1*1000
	shell:
		"""
		mv {input.predictions_extended} {output.plus_external_features}
		"""

# compute interaction or squared terms, fill NAs, rename features to finals, fill nas
rule gen_final_features:
	input:
		plus_external_features = os.path.join(RESULTS_DIR, "{biosample}", "ActivityOnly_external_features.tsv.gz"),
		feature_table_file = lambda wildcards: get_feature_table_file(wildcards.biosample)
	output:
		final_features = os.path.join(RESULTS_DIR, "{biosample}", "final_features.tsv.gz")
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=determine_mem_mb
	script:
		"../scripts/gen_final_features.R"
