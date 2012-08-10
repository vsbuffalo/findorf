## output_formatters.py -- strings for outputting objects

contig_sequence_repr = """
ContigSequence element for ID: $id
Length: $length
Number of relatives: $num_relatives

# Frames - these values are in [-3, -2, -1, 1, 2, 3]. GTF/GFF uses [0, 1, 2]
Majority frame: $majority_frame

# Frameshift
Majority frameshift: $majority_frameshift

# ORF Integrity
Missing start codon: $missing_start
Missing stop codon: $missing_stop
5'-end likely missing: $missing_5prime

# Predicted ORF - these values are 0-indexed
ORF is full length: $full_length_orf
ORF start: $orf_start
ORF stop: $orf_stop
ORF seq: $seq

# Relatives Start Sites

"""
