# libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(ggpubr)

# input files
#input_file = (snakemake@input$shap_scores)
#output_file = (snakemake@output$shap_plot)
#output_file_prec = (snakemake@output$comp_plot_prec)

shap_scores = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/results/2024_0501_extended_regularization/DNase_both_H3K27ac_intact/DNase_extended_noReg/model/shap_scores.tsv"
ft_val = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/results/2024_0501_extended_regularization/DNase_both_H3K27ac_intact/DNase_extended_noReg/model/training_data_in_order.tsv"
ft_table = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/results/2024_0501_extended_regularization/DNase_both_H3K27ac_intact/feature_table.tsv"
out_file = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/results/2024_0501_extended_regularization/DNase_both_H3K27ac_intact/DNase_extended_noReg/test_shap_2.pdf"

df_shap = fread(shap_scores, sep="\t")
df_ft = fread(ft_val, sep="\t")
ft_names = fread(ft_table, sep="\t")

## reformat
n = nrow(df_shap)
df_shap = pivot_longer(df_shap, everything(), names_to="feature", values_to="shap_score")
df_shap$shap_abs = abs(df_shap$shap_score)
# add feature values (and standardize)
df_ft_norm = df_ft %>% mutate_all(function(x) x / max(abs(x)))
df_ft = pivot_longer(df_ft, everything(), names_to="feature", values_to="feature_value")
df_ft_norm = pivot_longer(df_ft_norm, everything(), names_to="feature", values_to="feature_value_norm")
df_shap$feature_value = df_ft$feature_value
df_shap$feature_value_norm = df_ft_norm$feature_value_norm
# add feature_names (if not polynomial)
#### ADD IF STATEMENT
ft_names = dplyr::select(ft_names, feature, nice_name)
df_shap = left_join(df_shap, ft_names, by="feature")
df_shap$feature=df_shap$nice_name


## mean abs shap 
# summarize with 95% CI
summ = group_by(df_shap, feature) %>%
    summarize(mean_shap = mean(shap_abs), sd_shap = sd(shap_abs))
summ$CI_width = 1.96 * summ$sd_shap / sqrt(n)
summ$CI_low = summ$mean_shap-summ$CI_width
summ$CI_high = summ$mean_shap+summ$CI_width
# plot
x = ggplot(summ, aes(x=reorder(feature, mean_shap), y=mean_shap)) +
    geom_bar(stat="identity", position=position_dodge(width=0.8), width=0.7)+
    geom_errorbar(aes(ymin=CI_low, ymax=CI_high), position=position_dodge(width=0.8), width=0.25) +
    labs(x='Feature', y='Mean(|SHAP value|) (95% CI)') +
    coord_flip() +
    theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8))
#ggsave(filename=out_file, plot=x, width=6, height=8)

## feature value by shap
# sort by abs mean shap,
summ = summ[order(summ$mean_shap), ]
df_shap$feature = factor(df_shap$feature, ordered = TRUE, levels = summ$feature)

y = ggplot(df_shap, aes(x=feature, y=shap_score, color=feature_value_norm)) +
    geom_point(position = position_jitter(width = 0.1), alpha = 0.5, size=0.3) +
    labs(x='Feature', y='SHAP value') +
    coord_flip() +
    scale_color_gradient2(low="#006eae", high="#c5373d") +
    theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position='right')
ggsave(filename=out_file, plot=y, width=6, height=8)

## feature value by shap (scatter)
z = ggplot(df_shap, aes(x=feature_value, y=shap_score)) +
    geom_point(size=0.2, alpha=0.1) +
    geom_hline(yintercept=0, color="grey", linetype='dashed') +
    facet_wrap(vars(feature), ncol=6, axes="all", scales="free_x") +
    xlab("Feature value") + ylab("SHAP value") +
    theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position='top', strip.background = element_blank(), aspect.ratio=1)
#ggsave(out_file, z, width = 20, height = 20)

