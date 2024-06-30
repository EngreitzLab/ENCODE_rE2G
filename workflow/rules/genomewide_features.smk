from functools import partial

rule gen_new_features: 
	input:
		abc_predictions = lambda wildcards: os.path.join(ABC_BIOSAMPLES_DIR[wildcards.biosample], "Predictions", "EnhancerPredictionsAllPutative.tsv.gz"),
		enhancer_list = lambda wildcards: os.path.join(ABC_BIOSAMPLES_DIR[wildcards.biosample], "Neighborhoods", "EnhancerList.txt"),
	params:
		gene_TSS500 = config['gene_TSS500'],
		chr_sizes = config['chr_sizes'],
		scripts_dir = SCRIPTS_DIR
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=partial(determine_mem_mb, min_gb=32)
	output: 
		NumCandidateEnhGene = os.path.join(RESULTS_DIR, "{biosample}", "NumCandidateEnhGene.tsv"),
		NumTSSEnhGene = os.path.join(RESULTS_DIR, "{biosample}", "NumTSSEnhGene.tsv"),
		NumEnhancersEG5kb = os.path.join(RESULTS_DIR, "{biosample}", "NumEnhancersEG5kb.txt"),
		SumEnhancersEG5kb = os.path.join(RESULTS_DIR, "{biosample}", "SumEnhancersEG5kb.txt"),
	shell: 
		""" 
		python {params.scripts_dir}/feature_tables/gen_new_features.py \
			--enhancer_list {input.enhancer_list} \
			--abc_predictions {input.abc_predictions} \
			--ref_gene_tss {params.gene_TSS500} \
			--chr_sizes {params.chr_sizes} \
			--results_dir {RESULTS_DIR}/{wildcards.biosample}
		"""
		
# create activity-only feature table
rule activity_only_features:
	input:
		feature_table_file = os.path.join(RESULTS_DIR, "{biosample}", "feature_table.tsv"),
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
		"../scripts/feature_tables/activity_only_features.R"

# add external features
if config["final_score_col"] == "ENCODE-rE2G.Score.qnorm": # if sc-E2G
	min_mem = 32
else:
	min_mem = 8

rule add_external_features:
	input:
		predictions_extended = os.path.join(RESULTS_DIR, "{biosample}", "ActivityOnly_features.tsv.gz"),
		feature_table_file = os.path.join(RESULTS_DIR, "{biosample}", "feature_table.tsv"),
		external_features_config = os.path.join(RESULTS_DIR, "{biosample}", "external_features_config.tsv"),
	output:
		plus_external_features = os.path.join(RESULTS_DIR, "{biosample}",  "ActivityOnly_plus_external_features.tsv.gz")
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=partial(determine_mem_mb, min_gb=min_mem)  
	script:
		"../scripts/feature_tables/merge_external_features.R"

# compute interaction or squared terms, fill NAs, rename features to finals, fill nas
rule gen_final_features:
	input:
		plus_external_features = os.path.join(RESULTS_DIR, "{biosample}", "ActivityOnly_plus_external_features.tsv.gz"),
		feature_table_file = os.path.join(RESULTS_DIR, "{biosample}", "feature_table.tsv")
	output:
		final_features = os.path.join(RESULTS_DIR, "{biosample}", "genomewide_features.tsv.gz")
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=determine_mem_mb
	script:
		"../scripts/feature_tables/gen_final_features.R"
