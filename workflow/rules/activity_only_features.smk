# create activity-only feature table
rule activity_only_features:
	input:
		feature_table_file = lambda wildcards: get_feature_table_file(wildcards.biosample),
		abc = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", "EnhancerPredictionsAllPutative.tsv.gz"),
		NumCandidateEnhGene = os.path.join(RESULTS_DIR, "{biosample}", "NumCandidateEnhGene.tsv"),
		NumTSSEnhGene = os.path.join(RESULTS_DIR, "{biosample}", "NumTSSEnhGene.tsv"),
		NumEnhancersEG5kb = os.path.join(RESULTS_DIR, "{biosample}", "NumEnhancersEG5kb.txt"),
		SumEnhancersEG5kb = os.path.join(RESULTS_DIR, "{biosample}", "SumEnhancersEG5kb.txt"),
		ubiqExprGenes = config["ubiq_expr_genes"]
	params:
		activity = lambda wildcards: BIOSAMPLE_ACTIVITES[wildcards.biosample]
	output: 
		predictions_extended = os.path.join(RESULTS_DIR, "{biosample}", "ActivityOnly_features.tsv.gz")
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=32*1000
	script:
		"../scripts/activity_only_features.R"

rule gen_final_features:
	# We really only do a file renaming here
	# This step should be replaced by other modules if they wish to add more features
	input:
		predictions_extended = os.path.join(RESULTS_DIR, "{biosample}", "ActivityOnly_features.tsv.gz"),
	output:
		final_features = os.path.join(RESULTS_DIR, "{biosample}", "final_features.tsv.gz")
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=4*1000
	shell:
		"""
		mv {input.predictions_extended} {output.final_features}
		"""