# libraries
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

# if polynomial==False, use nice_name for plotting
if (polynomial==FALSE){
  ft_names = dplyr::select(feature_table, feature, nice_name)
  df = left_join(df, ft_names, by=c("feature_permuted"="feature"))
  df$feature_permuted = df$nice_name
}

# define y-label using n_repeats
y_label = paste0('Feature permuted (N=', n_repeats, ')')

## function that summarizes data nicely that I found online somewhere
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# plot delta auprc
delta_aupr_summ <- summarySE(df, measurevar="delta_aupr", groupvars=c("feature_permuted"))
delta_aupr_summ <- delta_aupr_summ[order(delta_aupr_summ$delta_aupr,decreasing=TRUE),]
delta_aupr_summ$feature_permuted <- factor(delta_aupr_summ$feature_permuted, levels = delta_aupr_summ$feature_permuted, ordered=TRUE)
ht = ifelse(nrow(delta_aupr_summ)<20, 3, 7) # set height of figure based on number of features

aupr = ggplot(delta_aupr_summ, aes(x=feature_permuted, y=delta_aupr)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=delta_aupr-ci, ymax=delta_aupr+ci), width=.2, position=position_dodge(.9)) +
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8)) +
  xlab('Delta AUPRC') + ylab(y_label) +
  coord_flip()

# plot delta precision
delta_precision_summ <- summarySE(df, measurevar="delta_precision", groupvars=c("feature_permuted"))
delta_precision_summ <- delta_precision_summ[order(delta_precision_summ$delta_precision,decreasing=TRUE),]
delta_precision_summ$feature_permuted <- factor(delta_precision_summ$feature_permuted, levels = delta_precision_summ$feature_permuted, ordered=TRUE)

# plot delta precision
prec = ggplot(delta_precision_summ, aes(x=feature_permuted, y=delta_precision)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=delta_precision-ci, ymax=delta_precision+ci), width=.2, position=position_dodge(.9)) +
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8)) +
  xlab('Delta precision at 70% recall') + ylab(y_label) +
  coord_flip()

ggsave(output_file_auprc, aupr, height=ht, width=4)
ggsave(output_file_prec, prec, height=ht, width=4)


