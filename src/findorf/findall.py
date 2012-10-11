"""
findall has methods for finding all ORFs, in all 6 reading frames, and
either outputting them or choosing the best based on KL divergence.
"""

import sys
from Bio import SeqIO
from collections import defaultdict
from utilities import get_codons, STOP_CODONS, START_CODONS
from RangedFeatures import ORF

def find_all_orfs(seq, min_length=30):
    """
    This function adds all ORFs to a dictionary organized by frame. It
    is very similar to get_all_orfs, except that it always starts in
    frame, and does this for every frame.
    """

    orfs = defaultdict(list)
    in_reading_frame = True
    query_length = len(seq)
    for frame in [-3, -2, -1, 1, 2, 3]:
        codons = get_codons(seq, frame)
        orf_start_pos = None
        query_start_pos = None

        for codon, orf_pos, query_pos in codons:
            codon = codon.upper()
            if codon in START_CODONS and not in_reading_frame:
                in_reading_frame = True
                orf_start_pos = orf_pos
                query_start_pos = query_pos
                continue
            if in_reading_frame and codon in STOP_CODONS:
                # a full reading frame, unless we haven't hit any stop
                # yet, then we're in the possible ORF from the start of
                # the query to the end.
                orf = ORF(query_start_pos, query_pos, query_length,
                          frame, no_start=query_start_pos is None,
                          no_stop=False)
                if len(orf) >= min_length:
                    orfs[frame].append(orf)

                in_reading_frame = False
                orf_start_pos = None
                query_start_pos = None

        # add any partial ORFs, and mark as having no stop codon
        if in_reading_frame:
            orf = ORF(query_start_pos, query_pos, query_length,
                      frame, no_start=query_start_pos is None,
                      no_stop=True)
            if len(orf) >= min_length:
                orfs[frame].append(orf)
    return orfs

 
def findall(seqfile, min_length=30, translate=False):
    with open(seqfile, 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            all_orfs = find_all_orfs(record.seq, min_length)
            # TODO add KL logic here
            for frame, orf_list in all_orfs.items():
                for orf in orf_list:
                    end = orf.end + 3 # include stop codon
                    start = orf.start
                    if orf.no_start:
                        start = abs(orf.frame) - 1

                    if frame < 0:
                        seq = record.seq.reverse_complement()[start:end]
                    else:
                        seq = record.seq[start:end]

                    if translate:
                        seq = seq.translate()
                    out = ">%s|frame=%s|start=%s\n%s\n" % (record.id, frame,
                                                              orf.abs_start(), seq)
                    sys.stdout.write(out)
