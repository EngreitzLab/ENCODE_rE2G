# libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(ggpubr)

# input files
input_file = (snakemake@input$comp_table)
output_file_auprc = (snakemake@output$comp_plot)

df = fread(input_file, sep="\t")

# plot
x = ggplot(df, aes(x = reorder(model, -AUPRC), y = AUPRC, fill = dataset)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = AUPRC_95CI_low, ymax = AUPRC_95CI_high), position = position_dodge(width = 0.8), width = 0.25) +
  labs(x = "Model", y = "AUPRC") +
  ylim(c(0, 1)) +
  guides(fill = guide_legend(title = "Dataset")) +
  coord_flip() +
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position='bottom')

# set colors
x = set_palette(x, 'npg')

ggsave(filename=output_file_auprc, plot=x, width=6, height=4)

