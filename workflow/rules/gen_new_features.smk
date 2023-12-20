from functools import partial

rule gen_new_features: 
	input:
		abc_predictions = lambda wildcards: os.path.join(ABC_dirs[wildcards.biosample], "Predictions", "EnhancerPredictionsAllPutative.tsv.gz"),
		enhancer_list = lambda wildcards: os.path.join(ABC_dirs[wildcards.biosample], "Neighborhoods", "EnhancerList.txt"),
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
		python {params.scripts_dir}/gen_new_features.py \
			--enhancer_list {input.enhancer_list} \
			--abc_predictions {input.abc_predictions} \
			--ref_gene_tss {params.gene_TSS500} \
			--chr_sizes {params.chr_sizes} \
			--results_dir {RESULTS_DIR}/{wildcards.biosample}
		"""