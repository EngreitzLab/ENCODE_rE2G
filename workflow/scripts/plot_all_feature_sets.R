# libraries
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(ggpubr)

# input files
#input_file = (snakemake@input$results)
#feature_table_file = (snakemake@input$feature_table)
#output_file = (snakemake@output$out_file)

input_file = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G/results/2024_0209_feature_analysis_noPenalty/DNase_megamap/DNase_megamap_noReg/feature_analysis/all_feature_sets.tsv"
feature_table_file =  "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G/resources/feature_tables/all_features_DNase_hic.tsv"
output_file = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G/results/2024_0209_feature_analysis_noPenalty/DNase_megamap/DNase_megamap_noReg/feature_analysis/all_feature_sets.pdf"

df = fread(input_file, sep="\t")
feature_table = fread(feature_table_file)

# plot auprc as a function of n_features
x = ggplot(df, aes(x=n_features, y=AUPRC, ymin=AUPRC_95CI_low, ymax=AUPRC_95CI_high)) +
    geom_pointrange(size=0.1, alpha=0.2, position=position_jitter(width=0.25, height=0)) +
    labs(x="Number of features", y="AUPRC (95% CI)") +
    scale_x_continuous(breaks=1:nrow(feature_table)) +
    theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8))

# calcs and formatting for feature-weighted auprc
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

# calculate "marginal contribution" to AUPRC (compare feature set with same one minus target)
all_features = feature_table$feature
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

# save plots
z = ggarrange(w, x, ncol=2, nrow=1)
ggsave(output_file, z, height=4, width=8)


