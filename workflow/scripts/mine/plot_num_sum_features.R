# libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(ggpubr)

# input files
new_num = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G/tests/test_output/generic/K562_chr22/new_features/NumEnhancersEG5kb.txt"
new_sum = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G/tests/test_output/generic/K562_chr22/new_features/SumEnhancersEG5kb.txt"
old_num = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G/temp_results/NumEnhancersEG5kb.txt"
old_sum = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G/temp_results/SumEnhancersEG5kb.txt"
out_num = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G/temp_results/compare_numNearbyEnhancers.pdf"
out_sum = "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G/temp_results/compare_sumNearbyEnhancers.pdf"


# cp
cp = c("#c5373d", "#006eae", "#429130")
names(cp) = c("genic", "intergenic", "promoter")

# numNearbyEnhancers
new_df = fread(new_num, sep="\t", header=FALSE) %>% setNames(c("name", "new_value"))
old_df = fread(old_num, sep="\t", header=FALSE) %>% setNames(c("name", "old_value"))
df = left_join(new_df, old_df, by="name") %>%
	separate_wider_delim(name, delim = "|", names = c("class", "location"))

n = ggplot(df, aes(x=old_value, y=new_value, color=class)) +
	geom_point(alpha=0.5) +
	xlab("Old numNearbyEnhancers") + ylab("New numNearbyEnhancers") +
	stat_cor(aes(label=..rr.label..)) +
	scale_colour_manual(values=cp) +
	theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position="bottom")
ggsave(filename=out_num, plot=n, width=5, height=5.2)

# sumNearbyEnhancers
new_df = fread(new_sum, sep="\t", header=FALSE) %>% setNames(c("name", "new_value"))
old_df = fread(old_sum, sep="\t", header=FALSE) %>% setNames(c("name", "old_value"))
df = left_join(new_df, old_df, by="name") %>%
	separate_wider_delim(name, delim = "|", names = c("class", "location"))

s = ggplot(df, aes(x=old_value, y=new_value, color=class)) +
	geom_point(alpha=0.5) +
	xlab("Old sumNearbyEnhancers") + ylab("New sumNearbyEnhancers") +
	stat_cor(aes(label=..rr.label..)) +
	scale_colour_manual(values=cp) +
	theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position="bottom")
ggsave(filename=out_sum, plot=s, width=5, height=5.2)