# libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(ggpubr)

# input files
input_file = (snakemake@input$comp_table)
output_file_auprc = (snakemake@output$comp_plot_auprc)
output_file_prec = (snakemake@output$comp_plot_prec)

df = fread(input_file, sep="\t")
use_colors = dplyr::n_distinct(df$dataset) <= 10 # only color by dataset if 10 or fewer datasets
if (!use_colors){
	df$model=paste0(df$dataset,"\n",df$model)
	df$dataset="filler"
}

# plot auprc
auprc = ggplot(df, aes(x = reorder(model, -AUPRC), y = AUPRC, fill = dataset)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = AUPRC_95CI_low, ymax = AUPRC_95CI_high), position = position_dodge(width = 0.8), width = 0.25) +
  labs(x = "Model", y = "AUPRC") +
  ylim(c(0, 1)) +
  guides(fill = guide_legend(title = "Dataset")) +
  coord_flip() +
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position='bottom')

# plot precision
prec = ggplot(df, aes(x = reorder(model, -precision), y = precision, fill = dataset)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = precision_95CI_low, ymax = precision_95CI_high), position = position_dodge(width = 0.8), width = 0.25) +
  labs(x = "Model", y = "Precision at 70% recall") +
  ylim(c(0, 1)) +
  guides(fill = guide_legend(title = "Dataset")) +
  coord_flip() +
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position='bottom')

if (!use_colors) {
    auprc = auprc + theme(legend.position='None')
    prec = prec + theme(legend.position='None')
}
auprc = set_palette(auprc, 'npg')
prec = set_palette(prec, 'npg')

ggsave(filename=output_file_auprc, plot=auprc, width=6, height=4)
ggsave(filename=output_file_prec, plot=prec, width=6, height=4)

