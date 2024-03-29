# reformat external features config
rule format_external_features_config:
	input:
		dataset_config = config["ABC_BIOSAMPLES"]
	output:
		external_features_config = os.path.join(RESULTS_DIR, "{dataset}", "external_features_config.tsv"),
	params:
		e2g_path = config["E2G_DIR_PATH"]
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=8*1000
	script:
		"../scripts/format_external_features_config.R"

# write feature table to output directory for future reference
rule write_feature_table:
    input:
        feature_table_file = lambda wildcards: get_feature_table_file(wildcards.biosample)
    output:
        feature_table_here = os.path.join(RESULTS_DIR, "{biosample}", "feature_table.tsv")
    run:
        shell("cat {input.feature_table_file} > {output.feature_table_here}")

# create activity-only feature table
rule activity_only_features:
	input:
		feature_table_file = lambda wildcards: get_feature_table_file(wildcards.biosample),
		abc = lambda wildcards: os.path.join(ABC_BIOSAMPLES_DIR[wildcards.biosample], "Predictions", "EnhancerPredictionsAllPutative.tsv.gz"),
		NumCandidateEnhGene = os.path.join(RESULTS_DIR, "{biosample}", "NumCandidateEnhGene.tsv"),
		NumTSSEnhGene = os.path.join(RESULTS_DIR, "{biosample}", "NumTSSEnhGene.tsv"),
		NumEnhancersEG5kb = os.path.join(RESULTS_DIR, "{biosample}", "NumEnhancersEG5kb.txt"),
		SumEnhancersEG5kb = os.path.join(RESULTS_DIR, "{biosample}", "SumEnhancersEG5kb.txt"),
		geneClasses = config["gene_classes"]
	output: 
		predictions_extended = os.path.join(RESULTS_DIR, "{biosample}", "ActivityOnly_features.tsv.gz")
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=determine_mem_mb
	script:
		"../scripts/activity_only_features.R"

# add external features
rule add_external_features:
	input:
		predictions_extended = os.path.join(RESULTS_DIR, "{biosample}", "ActivityOnly_features.tsv.gz"),
		feature_table_file = lambda wildcards: get_feature_table_file(wildcards.biosample),
		external_features_config = os.path.join(RESULTS_DIR, "{biosample}", "external_features_config.tsv"),
	output:
		plus_external_features = os.path.join(RESULTS_DIR, "{biosample}",  "ActivityOnly_plus_external_features.tsv.gz")
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=128*1000
	script:
		"../scripts/merge_external_features.R"

# compute interaction or squared terms, fill NAs, rename features to finals, fill nas
rule gen_final_features:
	input:
		plus_external_features = os.path.join(RESULTS_DIR, "{biosample}", "ActivityOnly_plus_external_features.tsv.gz"),
		feature_table_file = lambda wildcards: get_feature_table_file(wildcards.biosample)
	output:
		final_features = os.path.join(RESULTS_DIR, "{biosample}", "final_features.tsv.gz")
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=determine_mem_mb
	script:
		"../scripts/gen_final_features.R"
