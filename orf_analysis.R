## orf_analysis.R -- look at GTF from findorf
library(ggplot2)

GTF.FIELDS <- c("seqname", "source", "feature", "start",
                "end", "score", "strand", "frame", "group")

ANNOTATION.TERMS <- c("missing_stop", "has_orf", "inconsistent_strand", "internal_stop", 
                      "hsp_orf_overlap", "has_relatives", "missing_start", "majority_frameshift", 
                      "missing_5prime", "num_relatives", "closest_relative", "num_orf_candidates")

makeUnragged <- function(x) {
  group <- do.call(rbind, strsplit(x, " "))
  out <- group[match(ANNOTATION.TERMS, group[, 1]), 2]
  names(out) <- ANNOTATION.TERMS
  return(out)
}

importGTF <- function(file) {
  d <- read.table(file, sep="\t", header=FALSE, col.names=GTF.FIELDS,
                  na.strings=".",
                  colClasses=c("character", "factor", "factor", "numeric",
                    "numeric", "numeric", "factor", "factor", "character"))
  
  splits <- strsplit(as.character(d$group), "; ?", perl=TRUE)
  group <- do.call(rbind, lapply(splits, makeUnragged))
  tmp <- cbind(d[, -9], group)
  tmp$orf.length <- with(tmp, end-start)
  return(tmp)
}

d <- importGTF("src/findorf/nrev.gtf")

## verify some aspects of annotation
stopifnot(sum(table(d$has_orf)) == )

ggplot(d) + geom_bar(aes())
