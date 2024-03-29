rule make_dataset_feature_table:
	input:
		model_config = config["model_config"],
	output:
		dataset_features = os.path.join(RESULTS_DIR, "{dataset}", "feature_table.tsv")
	params:
		e2g_path = config["E2G_DIR_PATH"]
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=8*1000
	script:
		"../scripts/combine_feature_tables.R"

# reformat external features config
rule format_external_features_config:
	input:
		dataset_config = config["dataset_config"]
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

# add activity only features and retain columns from ABC for future features 
rule activity_only_features:
	input:
		abc = lambda wildcards: os.path.join(ABC_BIOSAMPLES_DIR[wildcards.dataset], "Predictions", "EnhancerPredictionsAllPutative.tsv.gz"),
		NumCandidateEnhGene = os.path.join(RESULTS_DIR, "{dataset}", "NumCandidateEnhGene.tsv"),
		NumTSSEnhGene = os.path.join(RESULTS_DIR, "{dataset}", "NumTSSEnhGene.tsv"),
		NumEnhancersEG5kb = os.path.join(RESULTS_DIR, "{dataset}", "NumEnhancersEG5kb.txt"),
		SumEnhancersEG5kb = os.path.join(RESULTS_DIR, "{dataset}", "SumEnhancersEG5kb.txt"),
		feature_table_file = os.path.join(RESULTS_DIR, "{dataset}", "feature_table.tsv"),
		geneClass = config["gene_classes"]
	output: 
		predictions_extended = os.path.join(RESULTS_DIR, "{dataset}", "ActivityOnly_features.tsv.gz")
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=32*1000
	script:
		"../scripts/activity_only_features.R"

# add external features
rule add_external_features:
	input:
		predictions_extended = os.path.join(RESULTS_DIR, "{dataset}", "ActivityOnly_features.tsv.gz"),
		feature_table_file = os.path.join(RESULTS_DIR, "{dataset}", "feature_table.tsv"),
		external_features_config = os.path.join(RESULTS_DIR, "{dataset}", "external_features_config.tsv"),
	output:
		plus_external_features = os.path.join(RESULTS_DIR, "{dataset}",  "ActivityOnly_plus_external_features.tsv.gz")
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=128*1000
	script:
		"../scripts/merge_external_features.R"

# compute interaction or squared terms, fill NAs, rename features to finals, fill nas
rule gen_final_features:
	input:
		plus_external_features = os.path.join(RESULTS_DIR, "{dataset}",  "ActivityOnly_plus_external_features.tsv.gz"),
		feature_table_file = os.path.join(RESULTS_DIR, "{dataset}", "feature_table.tsv")
	output:
		final_features = os.path.join(RESULTS_DIR, "{dataset}", "genomewide_features.tsv.gz")
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=32*1000
	script:
		"../scripts/gen_final_features.R"
	
# overlap feature table  with K562 CRISPR data
rule overlap_features_crispr:
	input:
		features = os.path.join(RESULTS_DIR, "{dataset}", "genomewide_features.tsv.gz"),
		crispr = config['crispr_dataset'],
		feature_table_file = os.path.join(RESULTS_DIR, "{dataset}", "feature_table.tsv"),
		tss = config['gene_TSS500']
	output: 
		features = os.path.join(RESULTS_DIR, "{dataset}", "EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_{nafill}.tsv.gz"),
		missing = os.path.join(RESULTS_DIR, "{dataset}",  "missing.EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_{nafill}.tsv.gz")
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
		crispr_features = os.path.join(RESULTS_DIR, "{dataset}", "EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_{nafill}.tsv.gz")
	params:
		genes = config["gene_TSS500"]
	output:
		processed = os.path.join(RESULTS_DIR, "{dataset}",  "for_training.EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_{nafill}.tsv.gz")
	conda:
		"../envs/encode_re2g.yml" 
	resources:
		mem_mb=32*1000
	script:
		"../scripts/process_crispr_data.R"

# rule validate_feature_table:
#     input: 
#          dataset_features = os.path.join(RESULTS_DIR, "{dataset}", "merged_feature_table.tsv")
#     output: 

# checkpoint validation:

# # validate features
# ref_features = config["reference_features"]
# for index, row in model_config.iterrows():
# 	dataset_row = dataset_config.loc[dataset_config["dataset"] == row["dataset"]]].squeeze() # corresponding row of datsaet_config as a series
# 	validate_feature_table(row["model"], row["feature_table"], dataset_row, ref_features)

# confirm external feature data is specified
	# for model application workflow, "entry" = biosample in biosample config, and biosample_row is its corresponding row
	# for model training workflow, "entry" = model in model_config, and biosample_row is the corresponding row in dataset_config
# def validate_feature_table(entry, feature_table, biosample_row, ref_features):
# 	extra_features = _determine_extra_features(feature_table, ref_features)
# 	if len(extra_features) > 0:
# 		if "external_feature_config" not in biosample_row.index.tolist():
# 			raise Exception("No external_feature_config is provided for these features required for " + entry + ": " + exra_features)
# 		no_data = _validate_external_features(extra_features, biosample_row["external_feature_config"])
# 		if len(no_data) > 0:
# 			raise Exception("Data for these external features are not present for " + entry + ": " + no_data)

# # check feature table features against reference and return list of extra features
# def _determine_extra_features(feature_table, ref_features):
# 	feature_table = pd.read_table(feature_table))
# 	req_features  = features_table[['input_col', 'second_input']].stack().dropna().unique().tolist()
# 	extra_features = [ft for ft in req_features if ft not in ref_features] # return things in feature table that aren't automatically generated
# 	return extra_features

# # check dataset config to see if extra features are deliniated in external_feature_config
# def _validate_external_features(extra_features, external_feature_config):
# 	efc = pd.read_table(external_feature_config)
# 	no_data = [ft for ft in extra_features if ft not in efc["input_col"]] 
# 	return no_data
	
############