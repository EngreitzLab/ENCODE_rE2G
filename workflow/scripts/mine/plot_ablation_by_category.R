# libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

# input files
base_dir = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/results/2024_0210_DNase_megamap_regularization/DNase_megamap/DNase_megamap_noReg/ablation_analysis_by_category_v2/"
#base_dir = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/results/2024_0501_extended_regularization/DNase_both_H3K27ac_intact/DNase_extended_noReg/ablation_analysis_by_category/"
#base_dir = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/results/2024_0329_final_models_two_inputs/DNase_ENCFF860XAE_megamap/DNase_ENCFF860XAE_megamap_8features/ablation_analysis_by_category/"

input_file = paste0(base_dir,"ablation_by_category_by_distance.tsv")
output_file_auprc =  paste0(base_dir,"ablation_auprc_distance.pdf")
output_file_prec =  paste0(base_dir,"ablation_precision_distance.pdf")
output_file_combo =  paste0(base_dir,"ablation_combined_all_distance.pdf")
output_file_combo_dist =  paste0(base_dir,"ablation_combined_by_distance.pdf")


df = fread(input_file, sep="\t")
only_plot_all = FALSE

# orders for plotting
df_all = dplyr::filter(df, distance_bin=="all") %>% arrange(delta_aupr)
if (only_plot_all){
	df = df_all
} else {
	df$distance_bin = factor(df$distance_bin, levels=c('all', '0-10 Kb', '10-100 Kb', '>100 Kb'), ordered=TRUE)
}
df$ablation_category = factor(df$ablation_category, levels=rev(df_all$ablation_category), ordered=TRUE)

ht = ifelse(nrow(df_all)<20, 3, 7) # set height of figure based on number of features
wd = 3+ length(levels(df$distance_bin)) # set width based on facets
 
# add some columns...
df <- df %>%
  mutate(
	annot_aupr = case_when(pval_aupr < 0.001 ~ '***',
	  pval_aupr >= 0.001 & pval_aupr < 0.01 ~ '**',
	  pval_aupr >= 0.01 & pval_aupr <= 0.05 ~ '*',
	  TRUE ~ "  "),
	color_aupr = ifelse(delta_aupr < 0 | pval_aupr>0.05, 'black', 'red'),
	loc_aupr = delta_aupr_low-0.03,
	annot_prec = case_when(pval_precision < 0.001 ~ '***',
						   pval_precision >= 0.001 & pval_precision < 0.01 ~ '**',
						   pval_precision >= 0.01 & pval_precision <= 0.05 ~ '*',
						   TRUE ~ '  '),
	color_prec = ifelse(delta_precision < 0 | pval_precision>0.05, 'black', 'red'),
	loc_prec = ifelse(delta_precision_high<0, delta_precision_low-0.05, delta_precision_high+0.05)
  )

# plot
x = ggplot(df, aes(x=ablation_category, y=delta_aupr)) +
facet_wrap(facet=vars(distance_bin), nrow=1) +
  geom_bar(stat="identity", fill="#96a0b3") +
  geom_errorbar(aes(ymin=delta_aupr_low, ymax=delta_aupr_high), width=0.5) +
  geom_text(aes(label=annot_aupr, y=loc_aupr, colour=color_aupr),size=2) +
  scale_colour_identity() +
  xlab('Category of features ablated') + ylab('Delta AUPRC from all feature model') +
  coord_flip() +
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position="none", strip.background = element_blank(), strip.placement = "outside")
ggsave(filename=output_file_auprc, plot=x, width=wd, height=ht)

y = ggplot(df, aes(x=ablation_category, y=delta_precision)) +
  facet_wrap(facet=vars(distance_bin), nrow=1) +
  geom_bar(stat="identity", fill="#435369") +
  geom_errorbar(aes(ymin=delta_precision_low, ymax=delta_precision_high), width=0.5) +
  geom_text(aes(label=annot_prec, y=loc_prec, colour=color_prec),size=2) +
  scale_colour_identity() +
  xlab('Category of features ablated') + ylab('Delta precision at 70% recall from all feature model') +
  coord_flip() +
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position="none", strip.background = element_blank(), strip.placement = "outside")
ggsave(filename=output_file_prec, plot=y, width=wd, height=ht)



col_names = c("ablation_category", "delta_x", "delta_x_low", "delta_x_high", "annot_x", "loc_x", "color_x", "distance_bin")
df_auprc = dplyr::select(df, ablation_category, delta_aupr, delta_aupr_low, delta_aupr_high, annot_aupr, loc_aupr, color_aupr, distance_bin)
colnames(df_auprc) = col_names
df_auprc$metric = "AUPRC"
df_prec = dplyr::select(df, ablation_category, delta_precision, delta_precision_low, delta_precision_high, annot_prec, loc_prec, color_prec, distance_bin)
colnames(df_prec) = col_names
df_prec$metric = "Precision at 70% recall"
df_long = rbind(df_auprc, df_prec)
df_long$nice_name = factor(df_long$ablation_category, levels=df_all$ablation_category, ordered=TRUE)

# plot combined across dist
w = ggplot(df_long, aes(x=ablation_category, y=delta_x, fill=metric)) +
  facet_wrap(facet=vars(distance_bin), nrow=1) +
  geom_bar(stat="identity", position=position_dodge(width=0.9)) +
  geom_errorbar(aes(ymin=delta_x_low, ymax=delta_x_high), linewidth=0.25, width=0.2, position=position_dodge(width=0.9)) +
  geom_text(aes(label=annot_x, y=loc_x, colour=color_x),size=2, position=position_dodge(width=0.9)) +
  scale_colour_identity() +
  scale_fill_manual(values=c("#96a0b3",	"#435369"), name="Metric") +
  xlab('Category of features ablated') + ylab('Delta from all feature model') +
  coord_flip() +
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
  legend.position="bottom", strip.background = element_blank(), strip.placement = "outside")
ggsave(filename=output_file_combo_dist, plot=w, width=wd, height=ht+1)

# combined, "all" distance bin only
df_order = dplyr::filter(df, distance_bin=="all") %>% arrange(delta_aupr)
df_all = dplyr::filter(df_long, distance_bin=="all")
df_all$ablation_category = factor(df_all$ablation_category, levels = rev(df_order$ablation_category), ordered=TRUE)

z = ggplot(df_all, aes(x=ablation_category, y=delta_x, fill=metric)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9)) +
  geom_errorbar(aes(ymin=delta_x_low, ymax=delta_x_high), linewidth=0.25, width=0.2, position=position_dodge(width=0.9)) +
  geom_text(aes(label=annot_x, y=loc_x, colour=color_x),size=2, position=position_dodge(width=0.9)) +
  scale_colour_identity() +
  scale_fill_manual(values=c("#96a0b3",	"#435369"), name="Metric") +
  xlab('Feature ablated') + ylab('Delta from all feature model') +
  coord_flip() +
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
  	legend.position="bottom", legend.direction="vertical", legend.text=element_text(size=7), legend.title=element_text(size=8))
ggsave(filename=output_file_combo, plot=z, width=5, height=ht+1 )