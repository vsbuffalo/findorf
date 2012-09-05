## orf_analysis.R -- look at GTF from findorf
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
  tmp$missing_start <- as.logical(tmp$missing_start)
  tmp$missing_stop <- as.logical(tmp$missing_stop)
  tmp$has_relatives <- as.logical(tmp$has_relative)
  tmp$hsp_orf_overlap <- as.logical(tmp$hsp_orf_overlap)
  tmp$inconsistent_strand <- as.logical(tmp$inconsistent_strand)
  tmp$majority_frameshift <- as.logical(tmp$majority_frameshift)
  tmp$missing_5prime <- as.logical(tmp$missing_5prime)
  tmp$internal_stop <- as.logical(tmp$internal_stop)
  return(tmp)
}
