rule e2g_variant_overlap:
	input:
		thresholded = os.path.join(RESULTS_DIR, "{biosample}", "{model_name}", "encode_e2g_predictions_threshold{threshold}.tsv.gz")
	params:
	 	score_column="ENCODE-rE2G.Score.qnorm",
		threshold="{threshold}",
		chr_sizes=config["chr_sizes"],
		scripts_dir = SCRIPTS_DIR
	output:
		BEDPE = os.path.join(
				RESULTS_DIR, 
				"{biosample}", 
				"{model_name}",
				"encode_e2g_predictions_threshold{threshold}.bedpe"),
		variantOverlap=os.path.join(
				RESULTS_DIR, 
				"{biosample}", 
				"{model_name}",
				"encode_e2g_predictions_threshold{threshold}.ForVariantOverlap.shrunk150bp.gz"),
	resources: 
		mem_mb=8*1000
	conda:
		"../envs/encode_re2g.yml"
	shell:
		"""
		python {params.scripts_dir}/model_application/process_e2g_output.py \
			--predictions_file {input.thresholded} \
			--score_column {params.score_column} \
			--chrom_sizes {params.chr_sizes} \
			--threshold {params.threshold} \
			--variant_overlap_output {output.variantOverlap} \
			--bedpe_output {output.BEDPE}
		"""	

rule combine_e2g_variant_overlap:
	input:
		variantOverlap=expand(os.path.join(
				RESULTS_DIR, 
				"{biosample}", 
				"{{model_name}}",
				"encode_e2g_predictions_threshold{{threshold}}.ForVariantOverlap.shrunk150bp.gz"), biosample=BIOSAMPLE_DF["biosample"].to_list())
	output:
		combined_variantOverlap = os.path.join(RESULTS_DIR, "combined_Predictions_{model_name}", "encode_e2g_predictions_threshold{threshold}.ForVariantOverlap.shrunk150bp.gz")
	resources:
		mem_mb=8*1000
	conda:
		"../envs/encode_re2g.yml"
	threads: 8
	shell:
		"""
		i=0
        for sample in {input.variantOverlap}
        do 
            if [ $i -eq 0 ]
            then
                zcat $sample |pigz -p{threads} > {output.combined_variantOverlap}
            else
                zcat $sample |sed 1d | pigz -p{threads} >> {output.combined_variantOverlap}
            fi
            ((i=i+1))
        done
		"""
