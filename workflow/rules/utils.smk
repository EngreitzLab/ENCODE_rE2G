import glob
import os
import pandas as pd
MAX_MEM_MB = 250 * 1000  # 250GB

def make_paths_absolute(obj, base_path):
	"""
	Use absolute paths to be compatible with github submodules
	Recursively go through the dictionary and convert relative paths to absolute paths.
	"""
	if isinstance(obj, dict):
		for key, value in obj.items():
			obj[key] = make_paths_absolute(value, base_path)
	elif isinstance(obj, str):
		# We assume all strings are paths. If converting the string
		# to an absolute path results in a valid file, then the str was a path
		new_file = os.path.join(base_path, obj)
		if os.path.exists(new_file):
			return new_file
	return obj

def determine_mem_mb(wildcards, input, attempt, min_gb=8):
	# Memory resource calculator for snakemake rules
	input_size_mb = input.size_mb
	if ".gz" in str(input):
		input_size_mb *= 8  # assume gz compressesed the file <= 8x
	attempt_multiplier = 2 ** (attempt - 1)  # Double memory for each retry
	mem_to_use_mb = attempt_multiplier *  max(4 * input_size_mb, min_gb * 1000)
	return min(mem_to_use_mb, MAX_MEM_MB)

def _get_biosample_model_dir(biosample):
	row = BIOSAMPLE_DF.loc[BIOSAMPLE_DF["biosample"] == biosample].iloc[0]
	access_type = row["default_accessibility_feature"].lower()
	hic_file = row["HiC_file"]
	if pd.isna(hic_file):
		raise Exception("No model found for powerlaw")
		# return os.path.join(MODEL_DIR, f"{access_type}_powerlaw")
	
	if pd.notna(row["H3K27ac"]):
		raise Exception("H3K27ac model not supported")
	
	if row["HiC_type"] == "avg":
		raise Exception("No model found for avg hic")
		# return os.path.join(MODEL_DIR, f"{access_type}_avg_hic")
	
	if hic_file == config["MEGAMAP_HIC_FILE"]:
		return os.path.join(MODEL_DIR, f"{access_type}_megamap")
	else:
		# assume intact hi-c
		if access_type=="atac":
			return os.path.join(MODEL_DIR, "multiome_intact_hic_test") # enable testing of sc-E2G pipeline
		return os.path.join(MODEL_DIR, f"{access_type}_intact_hic")

def get_feature_table_file(biosample):
	model_dir = _get_biosample_model_dir(biosample)
	return os.path.join(model_dir, "feature_table.tsv")

def get_trained_model(biosample):
	model_dir = _get_biosample_model_dir(biosample)
	return os.path.join(model_dir, "model.pkl")

def get_threshold(biosample):
	model_dir = _get_biosample_model_dir(biosample)
	print(model_dir)
	threshold_file = glob.glob(os.path.join(model_dir, 'threshold_*'))[0]
	threshold_file = os.path.basename(threshold_file)
	return threshold_file.split("_")[1]

# confirm external feature data is specified
	# for model application workflow, "entry" = biosample in biosample config, and biosample_row is its corresponding row
	# for model training workflow, "entry" = model in model_config, and biosample_row is the corresponding row in dataset_config
def validate_feature_table(entry, feature_table, biosample_row, ref_features):
	extra_features = _determine_extra_features(feature_table, ref_features)
	if len(extra_features) > 0:
		if "external_feature_config" not in biosample_row.index.tolist():
			raise Exception("No external_feature_config is provided for these features required for " + entry + ": " + exra_features)
		no_data = _validate_external_features(extra_features, biosample_row["external_feature_config"])
		if len(no_data) > 0:
			raise Exception("Data for these external features are not present for " + entry + ": " + no_data)

# check feature table features against reference and return list of extra features
def _determine_extra_features(feature_table, ref_features):
	feature_table = pd.read_table(feature_table))
	req_features  = features_table[['input_col', 'second_input']].stack().dropna().unique().tolist()
	extra_features = [ft for ft in req_features if ft not in ref_features] # return things in feature table that aren't automatically generated
	return extra_features

# check dataset config to see if extra features are deliniated in external_feature_config
def _validate_external_features(extra_features, external_feature_config):
	efc = pd.read_table(external_feature_config)
	no_data = [ft for ft in extra_features if ft not in efc["input_col"]] 
	return no_data



