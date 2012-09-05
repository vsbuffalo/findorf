library(knitr)
library(ggplot2)
source("orf-analysis.R")

d <- importGTF("/Users/vinceb/Projects/wheat-transcriptome/orf-prediction/Tdurum/predictions/out.gtf")

## plots for poster
d$frameshift.or.internal.stop <- as.factor(d$majority_frameshift | d$internal_stop)
levels(d$frameshift.or.internal.stop) <- c("False", "True")
p <- ggplot(d) + geom_density(aes(x=orf.length, fill=frameshift.or.internal.stop), alpha=0.4)
p <- p + scale_x_continuous("ORF Length", limits=c(0, 5000))
p <- p + scale_fill_discrete(name = "Frameshift or \ninternal stop codon")
p <- p + scale_y_continuous("Density")
ggsave("orf-lengths.pdf", plot=p, height=8, width=8)

p <- ggplot(subset(d, !is.na(closest_relative))) + geom_bar(aes(x=closest_relative, fill=closest_relative))

p <- p + scale_x_discrete("Closest Relative") + scale_y_continuous("Count") + opts(legend.position="none")
ggsave("relatives.pdf", plot=p, height=8, width=8)
