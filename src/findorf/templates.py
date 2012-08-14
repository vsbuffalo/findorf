## output_formatters.py -- strings for outputting objects

contig_sequence_repr = """
ContigSequence element for ID: $id
length: $length
number of relatives: $num_relatives

# Frames - these values are in [-3, -2, -1, 1, 2, 3]. GTF/GFF uses [0, 1, 2]
majority frame: $majority_frame

# Frameshift
majority frameshift: $majority_frameshift

# Anchor HSPs
$anchor_hsps

# ORF Integrity
missing start codon: $missing_start
missing stop codon: $missing_stop
5'-end likely missing: $missing_5prime
HSP coverage of ORF: $orf_hsp_coverage

# Predicted ORF
$orf

# Predicted ORF Sequence
seq: $seq

"""
    
orf_repr = """length: $length_bp bp, $length_aa aa
[start, end], in frame: [$start, $end]
[start, end], in query: [$query_start, $query_end]
frame: $frame
missing start/stop: $missing_start/$missing_stop"""

hsp_repr = """identities/length: $identities/$length
percent identity: $percent_identity
e-value: $e
frame: $frame
[query start, query end]: [$query_start - $query_end]
[subject start, start end]: [$sbjct_start, $sbjct_end]"""

anchor_hsps_repr = """strand: $strand
most 3':
$most_3prime
most 5':
$most_5prime"""

out = """
findorf summary
---------------------------

total: $total
full length orfs: $full_length
missing 5': $missing_5prime
missing start codon: $missing_start
missing stop codon: $missing_stop
closest relative frameshift: $cr_frameshift
majority frameshift: $majority_frameshift
any frameshift: $any_frameshift
contains stop codon: $contains_stop
"""