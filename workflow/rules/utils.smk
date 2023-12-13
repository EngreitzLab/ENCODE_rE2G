import glob
import os
import pandas as pd

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

def _get_biosample_model_dir(biosample):
	row = BIOSAMPLE_DF.loc[BIOSAMPLE_DF["biosample"] == biosample].iloc[0]
	access_type = row["default_accessibility_feature"].lower()
	hic_file = row["HiC_file"]
	if pd.isna(hic_file):
		raise Exception("No model found for powerlaw")
		# return os.path.join(E2G_DIR_PATH, "models", f"{access_type}_powerlaw")
	
	if pd.notna(row["H3K27ac"]):
		raise Exception("H3K27ac model not supported")
	
	if row["HiC_type"] == "avg":
		raise Exception("No model found for avg hic")
		# return os.path.join(E2G_DIR_PATH, "models", f"{access_type}_avg_hic")
	
	if hic_file == config["MEGAMAP_HIC_FILE"]:
		return os.path.join(E2G_DIR_PATH, "models", f"{access_type}_megamap")
	else:
		# assume intact hi-c
		return os.path.join(E2G_DIR_PATH, "models", f"{access_type}_intact_hic")

def get_feature_table_file(biosample):
	model_dir = _get_biosample_model_dir(biosample)
	return os.path.join(model_dir, "feature_table.tsv")

def get_trained_model(biosample):
	model_dir = _get_biosample_model_dir(biosample)
	return os.path.join(model_dir, "model.pkl")

def get_threshold(biosample):
	model_dir = _get_biosample_model_dir(biosample)
	threshold_file = glob.glob(os.path.join(model_dir, 'threshold_*'))[0]
	threshold_file = os.path.basename(threshold_file)
	return threshold_file.split("_")[1]
	