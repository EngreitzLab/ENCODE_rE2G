rule make_biosample_feature_table:  # make feature table per biosample
	input:
		config["ABC_BIOSAMPLES"]
	output:
		biosample_features = os.path.join(RESULTS_DIR, "{biosample}", "feature_table.tsv")
	params:
		model_dirs = lambda wildcards: BIOSAMPLE_DF.loc[BIOSAMPLE_DF['biosample'] == wildcards.biosample]['model_dir'].to_list(),
		#biosample_config = config["ABC_BIOSAMPLES_MODELS"],
		e2g_path = config["E2G_DIR_PATH"]
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=4*1000
	script:
		"../scripts/feature_tables/combine_feature_tables_apply.R"

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
		mem_mb=4*1000
	script:
		"../scripts/feature_tables/format_external_features_config.R"

rule generate_e2g_predictions:
	input:
		final_features = os.path.join(RESULTS_DIR, "{biosample}", "genomewide_features.tsv.gz"),
	params:
		epsilon = config["epsilon"],
		feature_table_file = lambda wildcards: get_feature_table_file(wildcards.biosample, wildcards.model_name),
		trained_model = lambda wildcards: get_trained_model(wildcards.biosample, wildcards.model_name),
		scripts_dir = SCRIPTS_DIR
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=determine_mem_mb
	output: 
		prediction_file = os.path.join(RESULTS_DIR, "{biosample}", "{model_name}", "encode_e2g_predictions.tsv.gz")
	shell: 
		""" 
		python {params.scripts_dir}/model_application/run_e2g.py \
			--predictions {input.final_features} \
			--feature_table_file {params.feature_table_file} \
			--epsilon {params.epsilon} \
			--trained_model {params.trained_model} \
			--output_file {output.prediction_file}
		"""

rule filter_e2g_predictions:
	input:
		prediction_file = os.path.join(RESULTS_DIR, "{biosample}", "{model_name}", "encode_e2g_predictions.tsv.gz")
	params:
		threshold = lambda wildcards: get_model_threshold(wildcards.biosample, wildcards.model_name),
		include_self_promoter = config["include_self_promoter"],
		scripts_dir = SCRIPTS_DIR
	conda:
		"../envs/encode_re2g.yml"
	resources:
		mem_mb=determine_mem_mb
	output:
		thresholded = os.path.join(RESULTS_DIR, "{biosample}", "{model_name}", "encode_e2g_predictions_threshold{threshold}.tsv.gz")
	shell:
		"""
		python {params.scripts_dir}/model_application/threshold_e2g_predictions.py \
			--all_predictions_file {input.prediction_file} \
			--threshold {params.threshold} \
			--include_self_promoter {params.include_self_promoter} \
			--output_file {output.thresholded}
		"""


rule combine_filter_e2g_predictions:
	input:
		thresholded=expand(os.path.join(RESULTS_DIR, "{biosample}", "{model_name}", "encode_e2g_predictions_threshold{{threshold}}.tsv.gz"), zip,biosample=BIOSAMPLE_DF["biosample"],model_name=BIOSAMPLE_DF["model_dir_base"])
	output:
		combined_thresholded=os.path.join(RESULTS_DIR, "combined_Predictions_{model_name}", "encode_e2g_predictions_threshold{threshold}.tsv.gz")
	resources:
		mem_mb=8*1000
	conda:
		"../envs/encode_re2g.yml"
	threads: 8
	shell:
		"""
		i=0
        for sample in {input.thresholded}
        do 
            if [ $i -eq 0 ]
            then
                zcat $sample |pigz -p{threads} > {output.combined_thresholded}
            else
                zcat $sample |sed 1d | pigz -p{threads} >> {output.combined_thresholded}
            fi
            ((i=i+1))
        done
		"""

