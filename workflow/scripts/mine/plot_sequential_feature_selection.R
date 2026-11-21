# libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

## input files
#base_dir = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G/results/2024_0209_feature_analysis_noPenalty/DNase_megamap/DNase_megamap_noReg/feature_analysis/sfs/"
#base_dir = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/results/2024_0421_extended_no_conserv/DNase_both_H3K27ac_intact/DNase_both_extended/feature_analysis/sfs/"
#dir.create(base_dir)
#ffs_file = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/results/2024_0421_extended_no_conserv/DNase_both_H3K27ac_intact/DNase_both_extended/feature_analysis/forward_feature_selection.tsv"
#feature_table_file = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/results/2024_0421_extended_no_conserv/DNase_both_H3K27ac_intact/feature_table.tsv"

base_dir = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G/results/2024_0209_feature_analysis_noPenalty/DNase_megamap/DNase_megamap_noReg/feature_analysis/sfs/"
ffs_file = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G/results/2024_0209_feature_analysis_noPenalty/DNase_megamap/DNase_megamap_noReg/feature_analysis/sfs/forward_feature_selection.tsv"
feature_table_file = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/resources/feature_tables/12features_DNase_hic.tsv"

output_file_auprc =  paste0(base_dir,"sfs_auprc.pdf")
output_file_prec =  paste0(base_dir,"sfs_precision.pdf")
output_file_combo =  paste0(base_dir,"sfs_combined.pdf")

df = fread(ffs_file, sep="\t")
ymin_aupr = min(-0.05, min(df$delta_aupr_low)-0.05)
ymin_prec = min(-0.05, min(df$delta_precision_low)-0.05)
max_aupr = max(df$aupr)
max_precision = max(df$precision)
aupr_null = df$aupr[1]
precision_null = df$precision[1]

ft = fread(feature_table_file, sep="\t") %>% dplyr::select(feature, nice_name)
df = inner_join(df, ft, by=c("feature_added"="feature"))
df$feature_added = df$nice_name

# order features for plotting
features = df$feature_added
df$feature_added = factor(df$feature_added, levels=rev(features), ordered=TRUE)
ht = ifelse(nrow(df)<20, 3, 7) # set height of figure based on number of features
 
# add some columns...
df <- df %>%
  mutate(
	annot_aupr = case_when(pval_aupr < 0.001 ~ '***',
	  pval_aupr >= 0.001 & pval_aupr < 0.01 ~ '**',
	  pval_aupr >= 0.01 & pval_aupr <= 0.05 ~ '*',
	  TRUE ~ "  "),
	color_aupr = ifelse(delta_aupr > 0 | pval_aupr>0.05, 'black', 'red'),
	loc_aupr = ifelse(delta_aupr_low<0, delta_aupr_low-0.03, delta_aupr_high+0.03),
	annot_prec = case_when(pval_precision < 0.001 ~ '***',
						   pval_precision >= 0.001 & pval_precision < 0.01 ~ '**',
						   pval_precision >= 0.01 & pval_precision <= 0.05 ~ '*',
						   TRUE ~ '  '),
	color_prec = ifelse(delta_precision > 0 | pval_precision>0.05, 'black', 'red'),
	loc_prec = ifelse(delta_precision_high<0, delta_precision_low-0.03, delta_precision_high+0.03)
  )

# plot individuals
ht = ifelse(nrow(df)<20, 3, 7) # set height of figure based on number of features
wd = 4.5

x = ggplot(df, aes(x=feature_added, y=delta_aupr)) +
  geom_hline(yintercept=max_aupr-aupr_null, linewidth=0.5, color='#c5cad7', linetype='dashed') +
  geom_bar(stat="identity", fill="#96a0b3") +
  geom_point(aes(x=feature_added, y=aupr-aupr_null), size=0.8, color="#96a0b3") +
  geom_errorbar(aes(ymin=delta_aupr_low, ymax=delta_aupr_high), width=0.5) +
  geom_text(aes(label=annot_aupr, y=loc_aupr, colour=color_aupr), size=2) +
  scale_colour_identity() +
  xlab('Feature added') + ylab('Bars: delta AUPRC from previous model\nDots: cumulative delta AUPRC') +
  coord_flip() +
  ylim(c(ymin_aupr, max_aupr-aupr_null+0.1)) +
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position="none")
ggsave(filename=output_file_auprc, plot=x, width=wd, height=ht)

y = ggplot(df, aes(x=feature_added, y=delta_precision)) +
  geom_hline(yintercept= max_precision-precision_null, linewidth=0.5, color='#c5cad7', linetype='dashed') +
  geom_bar(stat="identity") +
  geom_point(aes(x=feature_added, y=precision-precision_null), size=0.8, color="#435369") +
  geom_errorbar(aes(ymin=delta_precision_low, ymax=delta_precision_high), width=0.5) +
  geom_text(aes(label=annot_prec, y=loc_prec, colour=color_prec), size=2) +
  scale_colour_identity() +
  xlab('Feature added') + ylab('Bars: delta precision at 70% recall from previous model \nDots: cumulative delta precision at 70% recall') +
  coord_flip() +
  ylim(c(ymin_prec,max_precision-precision_null+0.1)) +
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position="none")
ggsave(filename=output_file_prec, plot=y, width=wd, height=ht)

# convert to long
col_names = c("feature_added", "x", "delta_x", "delta_x_low", "delta_x_high", "annot_x", "loc_x", "color_x")
df_auprc = dplyr::select(df, feature_added, aupr, delta_aupr, delta_aupr_low, delta_aupr_high, annot_aupr, loc_aupr, color_aupr)
colnames(df_auprc) = col_names
df_auprc = mutate(df_auprc, metric = "AUPRC", max_x = max_aupr, x_null = aupr_null)

df_prec = dplyr::select(df, feature_added, precision, delta_precision, delta_precision_low, delta_precision_high, annot_prec, loc_prec, color_prec)
colnames(df_prec) = col_names
df_prec = mutate(df_prec, metric =  "Precision at 70% recall", max_x = max_precision, x_null = precision_null)
df_long = rbind(df_auprc, df_prec)

df_long$feature_added = factor(df_long$feature_added, levels=rev(features), ordered=TRUE)

z = ggplot(df_long, aes(x=feature_added, y=delta_x, fill=metric)) +
  geom_hline(yintercept= max_precision-precision_null, linewidth=0.5, alpha=0.7, color='#435369', linetype='dashed') +  
  geom_hline(yintercept= max_aupr-aupr_null, linewidth=0.5, alpha=0.7, color='#96a0b3', linetype='dashed') +  
  geom_bar(stat="identity", position=position_dodge(width=0.9)) +
  geom_point(aes(x=feature_added, y=x-x_null, color=metric), size=0.8, position=position_dodge(width=0.9)) +
  geom_errorbar(aes(ymin=delta_x_low, ymax=delta_x_high), linewidth=0.25, width=0.2, position=position_dodge(width=0.9)) +
  geom_text(aes(label=annot_x, y=loc_x),size=2, position=position_dodge(width=0.9)) +
  scale_fill_manual(values=c("#96a0b3",	"#435369"), name="Metric") +
  scale_color_manual(values=c("#96a0b3",	"#435369"), name="Metric") +
  xlab('Feature added') + ylab('Bars: delta from previous model\nDots: cumulative delta') +
  ylim(c(min(ymin_prec, ymin_aupr), max(max_precision, max_aupr))) +
  coord_flip() +
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
  	legend.position="bottom", legend.direction="vertical", legend.text=element_text(size=7), legend.title=element_text(size=8))
ggsave(filename=output_file_combo, plot=z, width=5, height=ht+1 )