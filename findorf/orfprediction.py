"""
orfprediction.py - requisite functions for ORF finding and translation,
e.g.: getting codons, ORFs, etc.

Assume standard NCBI codon table.
"""

from Bio.Data import CodonTable
from Bio.Seq import Seq
from BioRanges.lightweight import Range, SeqRange, SeqRanges

CODON_TABLE = "Standard"
CODON_TABLE = CodonTable.unambiguous_dna_by_name[CODON_TABLE]
STOP_CODONS = set(CODON_TABLE.stop_codons)
START_CODONS = set(["ATG"])

def get_codons(seq, frame):
    """
    Return a list of (codon, position_in_orf, pos_in_forward_query)
    tuples, on the foward strand.
    """
    # for user friendliness with interactive shell, if seq is a
    # string, convert it to Bio.Seq.Seq
    if type(seq) is str:
        seq = Seq(seq)
        
    if frame < 0:
        seq = seq.reverse_complement()
        
    frame = abs(frame)
    tmp = [(str(seq[pos:pos+3]), pos-(frame-1), pos) for
           pos in range(frame-1, len(seq), 3)]
    
    # remove the last string if not a full codon.
    return [(c, po, pfq) for c, po, pfq in tmp if len(c) == 3]

def get_all_orfs(seqrecord, frame, in_reading_frame=False):
    """
    Generic ORF finder; it returns a list of all ORFs as they are
    found, given codons (a list if tuples in the form (codon,
    position)) from `get_codons`. This list is a list of SeqRange
    objects.

    If in_reading_frame is True, it ignores looking for a start codon
    first.

    Earlier versions did not allow for overlapping ORFs, since our
    approach was to always take the 5'-most start codon. However, now
    new versions will allow for other start codon prediction methods.

    For example:

        |-----M---M---M----------------*----|

    Earlier methods would just give a single ORF candidate:

                  |----------ORF-----------|
            |-----M---M---M----------------*----|

    Now, this new version will allow overlapping ORFs (via a queue):

                          |-----ORF3-------|
                      |---------ORF2-------|
                  |-------------ORF1-------|
            |-----M---M---M----------------*----|

    """
    seq = seqrecord.seq
    seqname = seqrecord.id
    seqlength = len(seq)

    # initialize ORF collector, and starting positions
    all_orfs = SeqRanges()
    orf_start_pos = None
    query_start_pos = None

    codons = get_codons(seq, frame)
    orf_frame = "+" if frame > 0 else "-"
    
    # Note that query_pos is forward strand.
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
            orf_data = {"no_start":query_start_pos is None, "no_stop":False}
            orf = SeqRange(Range(query_start_pos, query_pos), seqname,
                           orf_frame, seqlength=seqlength, data=orf_data)
            all_orfs.append(orf)
            
            # reset values
            in_reading_frame = False
            orf_start_pos = None
            query_start_pos = None
            
    # add any partial ORFs, and mark as having no stop codon
    if in_reading_frame:
        orf_data = {"no_start":query_start_pos is None, "no_stop":True}
        orf = SeqRange(Range(query_start_pos, query_pos), seqname,
                       orf_frame, seqlength=seqlength, data=orf_data)
        all_orfs.append(orf)
    return all_orfs

def predict_orf(seq_record, frame, method="1st-M", in_frame=False):
    """
    Prediction of ORF.

    method options are '1st-M' or 'conservative'.

     - '1st-M': Take the 5'-most M and use this

     - 'conservative': Take the 5'-most HSP and use the M 5' of this.
    """

    if method == "1st-M":
        pass
    elif method == "conservative":
        pass
    else:
        raise ValueError("method must be either '1st-M' or 'conservative'")
