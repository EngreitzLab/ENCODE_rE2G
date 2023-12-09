# create activity-only feature table
rule activity_only_features:
	input:
		feature_config = "config/feature_config.tsv",
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