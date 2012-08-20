## output_formatters.py -- strings for outputting objects

contig_str = """
Contig instance for ID: $id
length: $length
number of relatives: $num_relatives

# Frames - these values are in [-3, -2, -1, 1, 2, 3]. GTF/GFF uses [0, 1, 2]
majority frame: $majority_frame
inconsistent strand: $inconsistent_strand
majority frameshift: $majority_frameshift

# Closest Relative Anchor HSPs
relative: $relative
anchors:
$anchor_hsps
# ORF Integrity
missing start codon: $missing_start
missing stop codon: $missing_stop
5'-end likely missing: $missing_5prime
HSP coverage of ORF: $overlap

# Predicted ORF
$orf

# Predicted ORF Sequence
seq: $seq

"""
