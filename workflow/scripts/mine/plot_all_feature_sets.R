# libraries
suppressMessages({
    library(plyr)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(data.table)
    library(ggpubr)
	library(tibble)
})

# input files
#input_file = (snakemake@input$results)
#feature_table_file = (snakemake@input$feature_table)
#output_file = (snakemake@output$out_file)

# plot auprc as a function of n_features
auprc_by_n_features <- function(df, n_highlight){
    df$col = ifelse(df$n_features==n_highlight, '#3c5488', '#3c548855')
    x = ggplot(df, aes(x=factor(n_features), y=AUPRC, fill=col)) +
        geom_violin(color="transparent") +
        labs(x="Number of features", y="AUPRC") +
        ylim(0, 1) +
        scale_fill_identity() +
        theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position='None')

    return(x)
}

# calcs and formatting for feature-weighted auprc
feature_weighted_auprc <- function(df){
    df_ft = dplyr::select(df, -c(features, AUPRC_95CI_low, AUPRC_95CI_high)) %>%
        pivot_longer(cols=-c(AUPRC, n_features), names_to="feature", values_to="present") %>%
        mutate(feature_wt  = AUPRC*present/n_features) %>%
        group_by(feature) %>%
        summarize(AUPRC_sum = sum(feature_wt), AUPRC_sd=sd(feature_wt))

    # add nice names
    ft_names = dplyr::select(feature_table, feature, nice_name)
    df_ft = left_join(df_ft, ft_names, by='feature')

    # plot feature by "value"
    y = ggplot(df_ft, aes(x=nice_name, y=AUPRC_sum)) + 
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=AUPRC_sum-3*AUPRC_sd, ymax=AUPRC_sum+3*AUPRC_sd), width=.2, position=position_dodge(.9)) +
    theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8)) +
    labs(x="Feature", y="AUPRC when included (3*s.d.)") +
    coord_flip()

    return(y)
}

# calculate "marginal contribution" to AUPRC (compare feature set with same one minus target)
marginal_contribution <- function(df, feature_table){
    all_features = feature_table$feature
    ft_names = dplyr::select(feature_table, feature, nice_name)
    
    df_diff = data.frame()
    for (i in 1:length(all_features)){
        target_feature = all_features[i]
        others = all_features[all_features != target_feature]

        df_plus = dplyr::select(df, -c(features, AUPRC_95CI_low, AUPRC_95CI_high)) %>%
            dplyr::filter(get({{target_feature}})==1) %>%
            rename(AUPRC_plus = AUPRC)

        df_minus = dplyr::select(df, -c(features, AUPRC_95CI_low, AUPRC_95CI_high)) %>%
            dplyr::filter(get({{target_feature}})==0) %>%
            rename(AUPRC_minus=AUPRC)

        df_merge = left_join(df_plus, df_minus, by=others) %>%
            mutate(AUPRC_diff = AUPRC_plus - AUPRC_minus, feature=target_feature) %>%
            drop_na()
        df_merge = dplyr::select(df_merge, c(feature, AUPRC_plus, AUPRC_minus, AUPRC_diff))    
        df_diff = rbind(df_diff, df_merge)
    }
    df_diff = left_join(df_diff, ft_names, by='feature')
    avg_diff = df_diff %>% group_by(nice_name) %>%
        summarize(avg = mean(AUPRC_diff)) %>%
        arrange(avg)
    df_diff$nice_name = factor(df_diff$nice_name, ordered=TRUE, levels=avg_diff$nice_name)
    color_diff = ifelse(df_diff$AUPRC_diff < 0, 'red', 'black')

    w = ggplot(df_diff, aes(x=nice_name, y=AUPRC_diff)) + 
    geom_point(color=color_diff, position=position_jitter(width=0.25, height=0), alpha=0.2, size=0.1) +
    theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position='none') +
    labs(x="", y="Marginal contribution to AUPRC") +
    coord_flip()

    return(w)
}

with_and_without <- function(df, n_highlight, feature_table){
    all_features = feature_table$feature
    ft_names = dplyr::select(feature_table, feature, nice_name)

    df = dplyr::filter(df, n_features==n_highlight)
    df_long = pivot_longer(df, cols=all_of(all_features), names_to='feature', values_to='present')
    df_long$col = ifelse(df_long$present==1, '#00A087FF', '#91D1C2FF')
    df_long = left_join(df_long, ft_names, by='feature')

    # compare distributions per feature
    stat_test = df_long %>%
        group_by(nice_name) %>%
        summarize(p_value = t.test(AUPRC ~ present)$p.value)

    # test for ubiq expression
    with_ubiq = dplyr::filter(df, ubiquitousExpressedGene==1)
    without_ubiq = dplyr::filter(df, ubiquitousExpressedGene==0)
    #print(t.test(with_ubiq$AUPRC, without_ubiq$AUPRC))

    stat_test$stars = case_when(stat_test$p_value<0.001 ~ '***',
        stat_test$p_value<0.01 ~ '**',
        stat_test$p_value<=0.05 ~ '*',
        TRUE ~ '')
    #print(stat_test)

    test = compare_means(AUPRC ~ present, data=df_long, group.by='nice_name', method='t.test')
    print(test)

    x = ggplot(df_long, aes(x=nice_name, y=AUPRC, fill=col)) +
        geom_violin(color="transparent") +
        stat_compare_means(label='p.signif', hide.ns=TRUE, method='t.test') +
        xlab('') +
        scale_fill_identity() +
        theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position='None') +
        coord_flip()

    return(x)
}

feature_enrichment <- function(df, n_highlight, feature_table, output_dir){
	all_features = feature_table$feature
    ft_names = dplyr::select(feature_table, feature, nice_name)

    #df = dplyr::filter(df, n_features==n_highlight)
	c = quantile(df$AUPRC, probs=c(0.25, 0.75))
	print(df$AUPRC[1])
	print(c)

	df = dplyr::select(df, all_of(c(all_features, "AUPRC")))
	df_high = data.frame((colSums(dplyr::filter(df, AUPRC>c[2])))) %>% rownames_to_column(var="feature")
	colnames(df_high) = c("feature", "sum_high")
	df_low = data.frame(colSums(dplyr::filter(df, AUPRC<c[1]))) %>% rownames_to_column(var="feature")
	colnames(df_low) = c("feature", "sum_low")

	df_plot = left_join(ft_names, df_high, by="feature") %>% left_join(df_low, by="feature")
	df_plot$enr = df_plot$sum_high/(df_plot$sum_low+0.01)
	df_plot = arrange(df_plot, enr)
	df_plot$nice_name = factor(df_plot$nice_name, ordered=TRUE, levels=df_plot$nice_name)

	df_plot = pivot_longer(df_plot, cols=c("sum_high", "sum_low", "enr"), names_to='metric', values_to='value')
	df_plot = mutate(df_plot, col = case_when(metric=="sum_high" ~ '#00A087FF', metric=="sum_low" ~ '#91D1C2FF', metric=="enr" ~ '#b778b3'))
	df_enr = dplyr::filter(df_plot, metric=="enr")
	df_sum = dplyr::filter(df_plot, metric!="enr")

	x = ggplot(df_sum, aes(x=nice_name, y=value, fill=col)) +
        geom_bar(stat="identity", position="dodge") +
        xlab('') +
		ylab('Number of times feature is present in top/bottom 25% of models') +
        scale_fill_identity() +
        theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position='None') +
        coord_flip()
	
	y = ggplot(df_enr, aes(x=nice_name, y=value, fill=col)) +
        geom_bar(stat="identity", position="dodge") +
        xlab('') +
		ylab('Enrichment of feature in top/bottom 25% of models') +
        scale_fill_identity() +
        theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position='None') +
        coord_flip()

	ggsave(paste0(output_dir, "sum_all_0.25.pdf"), x, height=4, width=4)
	ggsave(paste0(output_dir, "enrichment_all_0.25.pdf"), y, height=4, width=4)

	#return(c(x, y))

}
heatmap <- function(df, n_highlight, feature_table, output_dir){
	all_features = feature_table$feature
    ft_names = dplyr::select(feature_table, feature, nice_name)

	df = dplyr::filter(df, n_features==8)
	df$index = 1:nrow(df)
	df =  dplyr::select(df, all_of(c(all_features, "AUPRC", "index"))) %>%
        pivot_longer(cols=-c(AUPRC, index), names_to="feature", values_to="present") %>%
		left_join(ft_names, by="feature") %>%
		dplyr::filter(present==1)

	h =  ggplot(df, aes(x=index, y=nice_name, fill=AUPRC)) +
        geom_tile() +
        xlab('') +
		ylab('') +
		scale_fill_gradient(low="#c5373d", high="#006eae") +
        theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position='bottom') 

	ggsave(paste0(output_dir, "heatmap_n8_v1.pdf"), h, height=5, width=5)
}

models_ranked <- function(df, n_highlight, output_dir){
	df$index = rev(1:nrow(df))
	df$col = ifelse(df$n_features==n_highlight, '#3c5488', '#3c548811')

	# distance
	dist= dplyr::filter(df, features=="['distanceToTSS']")
	dist_auprc = dist$AUPRC

	# chosen model
	df$col[df$features=="['numTSSEnhGene', 'distanceToTSS', 'normalizedDNase_prom', 'sumNearbyEnhancers', 'ubiquitousExpressedGene', 'numCandidateEnhGene', '3DContact', 'ABC.Score']"] = "#c5373d"
	chosen = dplyr::filter(df, features=="['numTSSEnhGene', 'distanceToTSS', 'normalizedDNase_prom', 'sumNearbyEnhancers', 'ubiquitousExpressedGene', 'numCandidateEnhGene', '3DContact', 'ABC.Score']")

	all = ggplot(df, aes(x=index, y=AUPRC)) +
		geom_point(size=0.5, alpha=0.4) +
		geom_point(size=1, aes(x=chosen$index, y=chosen$AUPRC), color="#c5373d")+
		#scale_color_identity() +
		geom_hline(yintercept=dist_auprc, color="#96a0b3", linetype="dotted") +
		xlab('Models ranked by AUPRC') + ylab('AUPRC') +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position='None') 

	df_sub = dplyr::filter(df, n_features==n_highlight)
	df_sub$index = rev(1:nrow(df_sub))
	sub = ggplot(df_sub, aes(x=index, y=AUPRC, color=col)) +
		#geom_point(size=0.5) +
		geom_pointrange(aes(y=AUPRC,ymin=AUPRC_95CI_low, ymax=AUPRC_95CI_high)) +
		scale_color_identity() +
		geom_hline(yintercept=dist_auprc, color="#96a0b3", linetype="dotted") +
		xlab('8 feature models ranked by AUPRC') + ylab('AUPRC (95% CI)') +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position='None')

	AUPRC_cutoff = quantile(df$AUPRC, probs=0.95)
	df_best = dplyr::filter(df, AUPRC>AUPRC_cutoff)
	best = ggplot(df_best, aes(x=n_features, y=AUPRC, ymin=AUPRC_95CI_low, ymax=AUPRC_95CI_high)) +
		#geom_point(size=0.5) +
		geom_pointrange(position=position_jitter(width=0.2), alpha=0.5) +
		geom_pointrange(aes(x=chosen$n_features, y=chosen$AUPRC, ymin=chosen$AUPRC_95CI_low, ymax=chosen$AUPRC_95CI_high), color="#c5373d") +
		scale_color_identity() +
		xlab('Number of features') + ylab('AUPRC (95% CI)') +
		ggtitle('Top 5% of models') +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position='None')

	ggsave(paste0(output_dir, "ranked_models_all.pdf"), all, height=4, width=4)
	ggsave(paste0(output_dir, "ranked_models_best.pdf"), best, height=4, width=4)
	#ggsave(paste0(output_dir, "ranked_models_n8.pdf"), sub, height=4, width=4)

}

input_file = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G/results/2024_0209_feature_analysis_noPenalty/DNase_megamap/DNase_megamap_noReg/feature_analysis/all_feature_sets.tsv"
feature_table_file =  "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G/resources/feature_tables/all_features_DNase_hic.tsv"
output_dir = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G/results/2024_0209_feature_analysis_noPenalty/DNase_megamap/DNase_megamap_noReg/feature_analysis/"

df = fread(input_file, sep="\t")
feature_table = fread(feature_table_file)

r = models_ranked(df, 8, output_dir)

#h =  heatmap(df, 8, feature_table, output_dir)

#m = feature_enrichment(df, 8, feature_table, output_dir)

# df_8 = dplyr::filter(df, n_features %in% c(7, 8))
# w = marginal_contribution(df_8, feature_table)
# ggsave(paste0(output_dir, "marginal_contribution_n8.pdf"), w, height=4, width=4)

#x = auprc_by_n_features(df, 8)
#ggsave(paste0(output_dir, "auprc_by_n.pdf"), x, height=4, width=4)

#y = with_and_without(df, 8, feature_table)
#ggsave(paste0(output_dir, "with_and_without.pdf"), y, height=6, width=6)

#z = ggarrange(w, x, ncol=2, nrow=1)
#ggsave(output_file, z, height=4, width=8)


