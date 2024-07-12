# libraries
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

# input files
input_file = (snakemake@input$results)
feature_table_file = (snakemake@input$feature_table)
polynomial = (snakemake@params$polynomial)
n_repeats = (snakemake@params$n_repeats)
output_file_auprc = (snakemake@output$out_auprc)
output_file_prec = (snakemake@output$out_prec)

df = fread(input_file, sep="\t")
feature_table = fread(feature_table_file)
polynomial = gsub(" ", "", polynomial)

# if polynomial==False, use nice_name for plotting
if (polynomial=="FALSE" || polynomial == "False"){
  ft_names = dplyr::select(feature_table, feature, nice_name)
  ft_names %>% add_row(feature='None', nice_name='None')
  df = left_join(df, ft_names, by=c("feature_permuted"="feature"))
	for (i in 1:nrow(df)){
		if (!is.na(df$nice_name[i]) || df$nice_name[i]!=""){
			df$feature_permuted[i] = df$nice_name[i]
		}
	}
}

# define y-label using n_repeats
y_label = paste0('Feature permuted (N=', n_repeats, ')')

calculate_summary_stats <- function(input_df, feature_column, value_column) {
  summary_df <- input_df %>%
    group_by({{feature_column}}) %>%
    summarize(
      mean_value = mean({{value_column}}, na.rm = TRUE),
      lower_ci = mean({{value_column}}, na.rm = TRUE) - qt(0.975, n() - 1) * sd({{value_column}}, na.rm = TRUE) / sqrt(n()),
      upper_ci = mean({{value_column}}, na.rm = TRUE) + qt(0.975, n() - 1) * sd({{value_column}}, na.rm = TRUE) / sqrt(n())
    ) 
  
  return(summary_df)
}
# plot delta auprc
delta_aupr_summ <- calculate_summary_stats(df, feature_column=feature_permuted, value_column=delta_aupr)
delta_aupr_summ <- delta_aupr_summ[order(delta_aupr_summ$mean_value,decreasing=TRUE),]
delta_aupr_summ$feature_permuted <- factor(delta_aupr_summ$feature_permuted, levels = delta_aupr_summ$feature_permuted, ordered=TRUE)
ht = ifelse(nrow(delta_aupr_summ)<20, 3, 7) # set height of figure based on number of features

aupr = ggplot(delta_aupr_summ, aes(x=feature_permuted, y=mean_value)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=.2, position=position_dodge(.9)) +
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8)) +
  xlab('Delta AUPRC') + ylab(y_label) +
  coord_flip()

# plot delta precision
delta_precision_summ <- calculate_summary_stats(df, feature_column=feature_permuted, value_column=delta_precision)
delta_precision_summ <- delta_precision_summ[order(delta_precision_summ$mean_value,decreasing=TRUE),]
delta_precision_summ$feature_permuted <- factor(delta_precision_summ$feature_permuted, levels = delta_precision_summ$feature_permuted, ordered=TRUE)

# plot delta precision
prec = ggplot(delta_precision_summ, aes(x=feature_permuted, y=mean_value)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=.2, position=position_dodge(.9)) +
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8)) +
  xlab('Delta precision at 70% recall') + ylab(y_label) +
  coord_flip()

ggsave(output_file_auprc, aupr, height=ht, width=4)
ggsave(output_file_prec, prec, height=ht, width=4)


