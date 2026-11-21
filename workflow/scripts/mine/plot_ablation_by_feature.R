# libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

## input files

## DNase 12 features
pred_dir = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/results/2024_0210_DNase_megamap_regularization/DNase_megamap/DNase_megamap_noReg/ablation_analysis_by_feature/"
feature_table_file =  file.path(pred_dir, "DNase_hic_12features.tsv")
input_file = file.path(pred_dir, "ablation_by_feature.tsv")

## DNase 8 features
# pred_dir = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/results/2024_0329_final_models_two_inputs/DNase_ENCFF860XAE_megamap/DNase_ENCFF860XAE_megamap_8features//ablation_analysis_by_feature/"
# feature_table_file = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/results/2024_0329_final_models_two_inputs/DNase_ENCFF860XAE_megamap/feature_table.tsv"
# input_file = file.path(pred_dir, "ablation_by_feature_by_distance.tsv")

output_file_auprc =file.path(pred_dir, "ablation_auprc_sort_aupr.pdf")
output_file_prec =file.path(pred_dir, "ablation_precison_sort_prec.pdf")
output_file_combo=file.path(pred_dir, "ablation_combined_sort_prec.pdf")

df = fread(input_file, sep="\t")
feature_table = fread(feature_table_file, sep="\t")
feature_table$feature_ablated=feature_table$feature
key = dplyr::select(feature_table,feature_ablated, nice_name)

# name and order features for plotting
df = dplyr::left_join(df, key, by="feature_ablated") %>% arrange(delta_precision) # %>% arrange(delta_aupr) ## %>% dplyr::filter(distance_bin=="all")
df$nice_name = factor(df$nice_name, levels=rev(df$nice_name), ordered=TRUE)
ht = ifelse(nrow(df)<20, 3, 7) # set height of figure based on number of features

# add some columns...
df <- df %>%
  mutate(
	annot_aupr = case_when(pval_aupr < 0.001 ~ '***',
	  pval_aupr >= 0.001 & pval_aupr < 0.01 ~ '**',
	  pval_aupr >= 0.01 & pval_aupr <= 0.05 ~ '*',
	  TRUE ~ "  "),
	color_aupr = ifelse(delta_aupr < 0 | pval_aupr > 0.05, 'black', 'red'),
	loc_aupr = delta_aupr_low-0.008,
	annot_prec = case_when(pval_precision < 0.001 ~ '***',
						   pval_precision >= 0.001 & pval_precision < 0.01 ~ '**',
						   pval_precision >= 0.01 & pval_precision <= 0.05 ~ '*',
						   TRUE ~ "  "),
	color_prec = ifelse(delta_precision < 0 | pval_aupr > 0.05, 'black', 'red'),
	loc_prec = delta_precision_low-0.008
  )

# plot
x = ggplot(df, aes(x=nice_name, y=delta_aupr)) +
  geom_bar(stat="identity", , fill="#96a0b3") +
  geom_errorbar(aes(ymin=delta_aupr_low, ymax=delta_aupr_high), width=0.5) +
  geom_text(aes(label=annot_aupr, y=loc_aupr, colour=color_aupr),size=2) +
  scale_colour_identity() +
  xlab('Feature ablated') + ylab('Delta AUPRC from all feature model') +
  coord_flip() +
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position="none")
#ggsave(filename=output_file_auprc, plot=x, width=4, height=ht)

y = ggplot(df, aes(x=nice_name, y=delta_precision)) +
  geom_bar(stat="identity",  fill = "#435369") +
  geom_errorbar(aes(ymin=delta_precision_low, ymax=delta_precision_high), width=0.5) +
  geom_text(aes(label=annot_prec, y=loc_prec, colour=color_prec),size=2) +
  scale_colour_identity() +
  xlab('Feature ablated') + ylab('Delta precision at 70% recall\nfrom all feature model') +
  coord_flip() +
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position="none")
ggsave(filename=output_file_prec, plot=y, width=4, height=ht)

# pivot longer...
col_names = c("nice_name", "delta_x", "delta_x_low", "delta_x_high", "annot_x", "loc_x", "color_x")
df_auprc = dplyr::select(df, nice_name, delta_aupr, delta_aupr_low, delta_aupr_high, annot_aupr, loc_aupr, color_aupr)
colnames(df_auprc) = col_names
df_auprc$metric = "AUPRC"
df_prec = dplyr::select(df, nice_name, delta_precision, delta_precision_low, delta_precision_high, annot_prec, loc_prec, color_prec)
colnames(df_prec) = col_names
df_prec$metric = "Precision at 70% recall"
df_long = rbind(df_auprc, df_prec)

z = ggplot(df_long, aes(x=nice_name, y=delta_x, fill=metric)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9)) +
  geom_errorbar(aes(ymin=delta_x_low, ymax=delta_x_high), linewidth=0.25, width=0.2, position=position_dodge(width=0.9)) +
  geom_text(aes(label=annot_x, y=loc_x, colour=color_x),size=2, position=position_dodge(width=0.9)) +
  scale_colour_identity() +
  scale_fill_manual(values=c("#96a0b3",	"#435369"), name="Metric") +
  xlab('Feature ablated') + ylab('Delta from all feature model') +
  coord_flip() +
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
  	legend.position="bottom", legend.direction="vertical", legend.text=element_text(size=7), legend.title=element_text(size=8))
ggsave(filename=output_file_combo, plot=z, width=5, height=ifelse(nrow(df)<20, 5, 9) )