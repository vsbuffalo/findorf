## sequtils.py -- sequeunce utilities

from Bio import SeqIO
from Bio.Seq import Seq

STOP_CODONS = set(("TAG", "TGA", "TAA"))
START_CODONS = set(("ATG"))
GTF_FIELDS = ("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "group")

def get_codons(seq_in_frame):
    """
    Return a list of (codon, position) tuples.
    """
    seq = seq_in_frame
    return [(seq[pos:pos+3], pos) for pos in range(0, len(seq), 3)]

def find_all_stops(codons):
    """
    Find all stop codons from a list of (codon, position) tuples.
    """
    stop_codons = list()
    for codon_tuple in codons:
        if codon_tuple[0].upper() in STOP_CODONS:
            stop_codons.append(codon_tuple)
    return stop_codons

def find_all_starts(codons):
    """
    Find all stop codons from a list of (codon, position) tuples.
    """
    start_codons = list()
    for codon_tuple in codons:
        if codon_tuple[0].upper() in START_CODONS:
            start_codons.append(codon_tuple)
    return start_codons
