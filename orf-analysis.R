## orf_analysis.R -- look at GTF from findorf
GTF.FIELDS <- c("seqname", "source", "feature", "start",
                "end", "score", "strand", "frame", "group")

ANNOTATION.TERMS <- c("missing_stop", "has_orf", "inconsistent_strand", "internal_stop", 
                      "hsp_orf_overlap", "has_relatives", "missing_start", "majority_frameshift",
                      "missing_5prime", "num_relatives", "closest_relative", "num_orf_candidates",
                      "contig_length")

unfactor <- function(f) as.numeric(levels(f))[f]

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
  tmp$orf_length <- with(tmp, end-start)
  tmp$missing_start <- as.logical(tmp$missing_start)
  tmp$missing_stop <- as.logical(tmp$missing_stop)
  tmp$has_relatives <- as.logical(tmp$has_relative)
  tmp$hsp_orf_overlap <- as.logical(tmp$hsp_orf_overlap)
  tmp$inconsistent_strand <- as.logical(tmp$inconsistent_strand)
  tmp$majority_frameshift <- as.logical(tmp$majority_frameshift)
  tmp$missing_5prime <- as.logical(tmp$missing_5prime)
  tmp$internal_stop <- as.logical(tmp$internal_stop)
  tmp$contig_length <- unfactor(tmp$contig_length)
  
  # append a column annotating fullness
  tmp$status <- NA
  tmp$status[tmp$missing_start & !tmp$missing_stop] <- "missing start only"
  tmp$status[tmp$missing_stop & !tmp$missing_start] <- "missing stop only"
  tmp$status[tmp$missing_start & tmp$missing_stop] <- "missing start and stop"
  tmp$status[!tmp$missing_start & !tmp$missing_stop] <- "full length"
  # TODO check these cases here:
  tmp$status[is.na(tmp$status) & !tmp$hsp_orf_overlap] <- "no ORF/HSP overlap"
  tmp$status[is.na(tmp$status) & tmp$inconsistent_strand] <- "inconsistent strand"
  
  tmp$alt_status <- NA
  tmp$alt_status[tmp$missing_5prime] <- "missing 5'"
  tmp$alt_status[!tmp$missing_5prime] <- "not missing 5'"
  tmp$alt_status[tmp$missing_5prime & tmp$majority_frameshift] <- "missing 5' & majority frameshift"
  tmp$alt_status[!tmp$missing_5prime & tmp$majority_frameshift ] <- "majority frameshift"
  tmp$alt_status[is.na(tmp$alt_status) & tmp$internal_stop] <- "internal stop"
  tmp$alt_status[!tmp$missing_5prime & !tmp$majority_frameshift & !tmp$internal] <- "no missing 5', frameshift, or internal stop"

  # add imputed 3' and 5'-UTR lengths
  tmp$length.3utr <- NA
  tmp$length.5utr <- NA
  tmp$length.3utr <- tmp$contig_length - tmp$end
  tmp$length.5utr <- tmp$start
  
  return(tmp)
}
