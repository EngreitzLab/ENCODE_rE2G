# libraries
suppressMessages({
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(ggpubr)
library(ggbeeswarm)
})

# input files
#input_file = (snakemake@input$shap_scores)
#output_file = (snakemake@output$shap_plot)
#output_file_prec = (snakemake@output$comp_plot_prec)

# extended
dataset_dir = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/results/2024_0501_extended_regularization/DNase_both_H3K27ac_intact"
model_file_name = "DNase_extended_noReg"

# DNase + megamap 
# dataset_dir = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/results/2024_0210_DNase_megamap_regularization/DNase_megamap/"
# model_file_name = "DNase_megamap_noReg"

shap_scores = file.path(dataset_dir, model_file_name, "model", "shap_scores.tsv")
#ft_table = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/results/2024_0210_DNase_megamap_regularization/DNase_megamap/feature_table.tsv"
ft_table = file.path(dataset_dir, "feature_table.tsv")
ft_val = file.path(dataset_dir, model_file_name, "model", "training_data_in_order.tsv")
out_dir = file.path(dataset_dir, model_file_name, "shap_scores_new")
dir.create(out_dir)

df_shap = fread(shap_scores, sep="\t")
df_ft = fread(ft_val, sep="\t")
ft_names = fread(ft_table, sep="\t")
colnames(df_shap) = ft_names$feature; #head(df_shap)
colnames(df_ft) = ft_names$feature;# head(df_ft)

## reformat and calculate things
n = nrow(df_shap)
n_features = nrow(ft_names)

# scale shap scores based on range of feature (non-norm)?
df_shap_norm = df_shap
for (i in 1:ncol(df_shap)){
	ft_range = max(df_ft[[i]]) - min(df_ft[[i]])
	df_shap_norm[[i]] = df_shap[[i]] / ft_range
}
df_shap$id = 1:nrow(df_shap)

df_shap = pivot_longer(df_shap, -id, names_to="feature", values_to="shap_score")
df_shap$shap_abs = abs(df_shap$shap_score)
df_shap_norm = pivot_longer(df_shap_norm, everything(), names_to="feature", values_to="shap_score")
df_shap$shap_score_norm = df_shap_norm$shap_score
df_shap$shap_norm_abs = abs(df_shap$shap_score_norm)

# add feature values and scale to 0-1
df_ft_norm = df_ft %>% mutate_all(function(x) (x-min(x))/ (max(x)-min(x)))
df_ft = pivot_longer(df_ft, everything(), names_to="feature", values_to="feature_value")
df_ft$feature_value_original = exp(df_ft$feature_value) - 0.01
df_ft_norm = pivot_longer(df_ft_norm, everything(), names_to="feature", values_to="feature_value_norm")
df_shap$feature_value = df_ft$feature_value
df_shap$feature_value_norm = df_ft_norm$feature_value_norm
df_shap$feature_value_original = df_ft$feature_value_original

# add feature_names (assume not polynomial)
ft_names = dplyr::select(ft_names, feature, nice_name)
df_shap = left_join(df_shap, ft_names, by="feature")
df_shap$feature=df_shap$nice_name

## mean abs shap 
# summarize with 95% CI
summ = group_by(df_shap, feature) %>%
    summarize(mean_shap = mean(shap_abs), sd_shap = sd(shap_abs),
						mean_shap_norm = mean(shap_norm_abs), sd_shap_norm=sd(shap_norm_abs)) %>%
	arrange(mean_shap, descending=TRUE)
summ$feature = factor(summ$feature, levels=summ$feature,ordered=TRUE)

summ$CI_width = 1.96 * summ$sd_shap / sqrt(n)
summ$CI_low = summ$mean_shap-summ$CI_width
summ$CI_high = summ$mean_shap+summ$CI_width

# plot
x = ggplot(summ, aes(x=feature, y=mean_shap)) +
    geom_bar(stat="identity", position=position_dodge(width=0.8), width=0.7)+
    geom_errorbar(aes(ymin=CI_low, ymax=CI_high), position=position_dodge(width=0.8), width=0.25) +
    labs(x='Feature', y='Mean(|SHAP value|) (95% CI)') +
    coord_flip() +
    theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8))
#ggsave(file=file.path(out_dir, "mean_shap_scores.pdf"), plot=x, width=5, height=ifelse(n_features>20, 6, 3), create.dir=TRUE)

## feature value by shap
# sort by abs mean shap
df_shap$feature = factor(df_shap$feature, levels=summ$feature, ordered=TRUE)
y = ggplot(df_shap, aes(x=feature, y=shap_score, color=feature_value_norm)) +
    #geom_point(position = position_jitter(width = 0.1), alpha = 0.5, size=0.3) +
	geom_quasirandom(size=0.8, alpha=0.8, method = "tukeyDense", shape=".") +
    labs(x='Feature', y='SHAP value', color = "Feature value") +
    coord_flip() +
    scale_color_gradient(low="#006eae", high="#c5373d", oob=scales::squish, na.value="#FFFFFF", limits=c(0.1, 0.9), breaks = c(0.1, 0.9), labels = c("Low", "High")) +
    theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position='bottom') + 
	guides(color = guide_colourbar(theme = theme(legend.key.width=unit(8, "lines"), legend.key.height = unit(0.8, "lines"), legend.ticks=element_blank()), title.position = "top"))
#ggsave(filename=file.path(out_dir, "shap_by_feature_value.pdf"), plot=y, width=5, height=ifelse(n_features>20, 9, 4))

## feature value by shap (scatter)
feature_to_color = "Promoter class (Bergman et al. 2022)"
color_key = dplyr::filter(df_shap, feature==feature_to_color) %>% dplyr::select(id, feature_value_norm) %>%
	mutate(feature_color= feature_value_norm) %>% dplyr::select(id, feature_color)
df_shap = left_join(df_shap, color_key, by="id")

z = ggplot(df_shap, aes(x=feature_value, y=shap_score, color=feature_color)) +
    geom_point(size=0.8, alpha=0.3, shape='.') +
    geom_hline(yintercept=0, color="grey", linetype='dashed') +
    facet_wrap(vars(feature), ncol=5, axes="all", scales="free_x") +
	scale_color_gradient(low="#006eae", high="#c5373d", oob=scales::squish, na.value="#FFFFFF", limits=c(0.1, 0.9), breaks = c(0.1, 0.9), labels = c("Low", "High")) +
    labs(x="Feature value", y="SHAP value", color=feature_to_color) +
    theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position='bottom', strip.background = element_blank(), aspect.ratio=1)
ggsave(file.path(out_dir, "shap_score_scatter_P2.pdf"), z, width = 10, height = 2*ceiling(nrow(ft_names)/5))

