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

def expand_biosample_df(biosample_df):
	# add new columns
	if "model_dir" not in biosample_df.columns:
		biosample_df['model_dir']  = np.nan
	biosample_df['model_dir_base'] = ''
	biosample_df['model_threshold'] = float(0)

	new_rows = []
	for index, row in biosample_df.iterrows():
		if pd.isna(row['model_dir']):
			biosample_df.loc[index, "model_dir"] = os.path.normpath(_get_biosample_model_dir_from_row(row)) # infer model directory if not specified
		
		for model_dir in biosample_df.loc[index, 'model_dir'].split(','):
			new_row = row.copy()
			new_row['model_dir'] = model_dir
			new_row['model_dir_base'] = os.path.basename(model_dir)
			new_rows.append(new_row)
	new_df = pd.DataFrame(new_rows)
	new_df['model_threshold'] = [float(get_model_threshold(this_biosample, this_model, new_df)) for this_biosample, this_model in zip(new_df['biosample'], new_df['model_dir_base'])]

	return(new_df)

def _validate_model_dir(potential_dir):
	files = os.listdir(potential_dir)
	if "model.pkl" not in files:
		raise Exception("model.pkl is provided in the specified model directory")
	if "feature_table.tsv" not in files:
		raise Exception("feature_table.tsv is not provided in the specified model directory") 
	if sum(s.startswith("threshold_") for s in files) == 0:
		raise Exception("A threshold is not provided in the specified model directory")
	elif sum(s.startswith("threshold_") for s in files) > 1:
		raise Exception("More than one threshold is provided in the specified model directory")

def _get_model_dir_from_wildcards(biosample, model_name, biosample_df=None):
	if biosample_df is None:
		df = BIOSAMPLE_DF
	else:
		df = biosample_df
	model_dir = df.loc[(df["biosample"]==biosample) & (df["model_dir_base"]==model_name), "model_dir"]
	return (model_dir.values[0])

def _get_biosample_model_dir_from_row(row):
	# if user has explicitly defined model_dir, use by default
	if "model_dir" in BIOSAMPLE_DF.columns:
		if pd.notna(row["model_dir"]):
			_validate_model_dir(row["model_dir"])
			return row["model_dir"]

	access_type = row["default_accessibility_feature"].lower()
	input_params = [access_type]
	if not pd.isna(row["H3K27ac"]):
		input_params.append("h3k27ac")

	hic_file = row["HiC_file"]
	if pd.isna(hic_file):
		input_params.append("powerlaw")
	elif row["HiC_type"] == "avg":
		input_params.append("avg_hic")
	elif os.path.basename(hic_file) == os.path.basename(config["MEGAMAP_HIC_FILE"]):
		input_params.append("megamap")
	elif row["HiC_type"] == "hic":
		input_params.append("intact_hic")
	
	model_folder = os.path.join(MODEL_DIR, "_".join(input_params))
	if not os.path.exists(model_folder):
		raise Exception(f"{model_folder} not found. Model with input params not supported")
	return model_folder

def get_feature_table_file(biosample, model_name): 
	model_dir = _get_model_dir_from_wildcards(biosample, model_name)
	return os.path.join(model_dir, "feature_table.tsv")

def get_trained_model(biosample, model_name):
	model_dir = _get_model_dir_from_wildcards(biosample, model_name)
	return os.path.join(model_dir, "model.pkl")

def get_model_threshold(biosample, model_name, biosample_df=None):
	model_dir = _get_model_dir_from_wildcards(biosample, model_name, biosample_df)
	threshold_files = glob.glob(os.path.join(model_dir, 'threshold_*'))
	assert len(threshold_files) == 1, "Should have exactly 1 threshold file in directory"
	threshold_file = os.path.basename(threshold_files[0])
	return threshold_file.split("_")[1]
