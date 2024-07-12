# libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

# input files
input_file = (snakemake@input$results)
feature_table_file = (snakemake@input$feature_table)
polynomial = (snakemake@params$polynomial)
output_file_auprc = (snakemake@output$out_auprc)
output_file_prec = (snakemake@output$out_prec)

df = fread(input_file, sep="\t")
feature_table = fread(feature_table_file)
polynomial = gsub(" ", "", polynomial)

# if polynomial==False, use nice_name for plotting
if (polynomial=="FALSE" || polynomial == "False"){
  ft_names = dplyr::select(feature_table, feature, nice_name)
  df = left_join(df, ft_names, by=c("feature_added"="feature"))
    for (i in 1:nrow(df)){
		if (!is.na(df$nice_name[i]) || df$nice_name[i]!=""){
			df$feature_added[i] = df$nice_name[i]
		}
	}
}

# order features for plotting
features = df$feature_added
df$feature_added = factor(df$feature_added, levels=rev(features), ordered=TRUE)
ht = ifelse(nrow(df)<20, 3, 7) # set height of figure based on number of features
ymin_aupr = min(-0.05, min(df$delta_aupr_low)-0.01)
ymin_prec = min(-0.05, min(df$delta_precision_low)-0.01)
full_aupr = df$aupr[nrow(df)]
full_prec = df$precision[nrow(df)]

# add some columns...
df <- df %>%
  mutate(
    annot_aupr = case_when(pval_aupr < 0.001 ~ '***',
                           pval_aupr >= 0.001 & pval_aupr < 0.01 ~ '**',
                           pval_aupr >= 0.01 & pval_aupr <= 0.05 ~ '*',
                           TRUE ~ ''),
    color_aupr = ifelse(delta_aupr < 0, 'red', 'black'),
    loc_aupr = ifelse(delta_aupr_high < 0, 0.05, delta_aupr_high+0.05),
    annot_prec = case_when(pval_precision < 0.001 ~ '***',
                           pval_precision >= 0.001 & pval_precision < 0.01 ~ '**',
                           pval_precision >= 0.01 & pval_precision <= 0.05 ~ '*',
                           TRUE ~ ''),
    color_prec = ifelse(delta_precision < 0, 'red', 'black'),
    loc_prec = ifelse(delta_precision_high < 0, 0.05, delta_precision_high+0.05),
  )

# plot
x = ggplot(df, aes(x=feature_added, y=delta_aupr)) +
  geom_hline(yintercept=full_aupr, linewidth=0.5, color='gray', linetype='dotted') +
  geom_bar(stat="identity") +
  geom_point(aes(x=feature_added, y=aupr), size=0.8) +
  geom_errorbar(aes(ymin=delta_aupr_low, ymax=delta_aupr_high), width=0.5) +
  geom_text(aes(label=annot_aupr, y=loc_aupr, colour=color_aupr), size=2) +
  scale_colour_identity() +
  xlab('Feature added') + ylab('Delta AUPRC') +
  coord_flip() +
  ylim(c(ymin_aupr, 1)) +
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position="none")
ggsave(filename=output_file_auprc, plot=x, width=4, height=ht)

y = ggplot(df, aes(x=feature_added, y=delta_precision)) +
  geom_hline(yintercept=full_prec, linewidth=0.5, color='gray', linetype='dotted') +
  geom_bar(stat="identity") +
  geom_point(aes(x=feature_added, y=precision), size=0.8) +
  geom_errorbar(aes(ymin=delta_precision_low, ymax=delta_precision_high), width=0.5) +
  geom_text(aes(label=annot_prec, y=loc_prec, colour=color_prec), size=2) +
  scale_colour_identity() +
  xlab('Feature added') + ylab('Delta precision at 70% recall') +
  coord_flip() +
  ylim(c(ymin_prec, 1)) +
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position="none")
ggsave(filename=output_file_prec, plot=y, width=4, height=ht)


