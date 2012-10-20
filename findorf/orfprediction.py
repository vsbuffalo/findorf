"""
orfprediction.py - requisite functions for ORF finding and translation,
e.g.: getting codons, ORFs, etc.

Assume standard NCBI codon table.
"""

from collections import deque
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

def get_all_orfs(seqrecord, frame):
    """
    Generic ORF finder; it returns a list of all ORFs as they are
    found, given codons (a list if tuples in the form (codon,
    position)) from `get_codons`. This list is a list of SeqRange
    objects.

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

    Save for annotating them differently, this function will not
    ignore partial 5'-incomplete for 3'-incomplete cases. In other
    words, we assume we are in an open reading frame from the
    start. This allows us to enumerate every biologically possible
    ORF.
    
    """
    seq = seqrecord.seq
    seqname = seqrecord.id
    seqlength = len(seq)

    # initialize ORF collector, and starting positions
    all_orfs = SeqRanges() # for final ORFs

    # to handle keeping many reading possible reading frame candidates
    # open at once, we use a queue. Tuples maintain key data:
    # (start orf position, start query position, whether had start codon)
    orf_queue = deque(list())

    # get all codons in frame
    codons = get_codons(seq, frame)
    orf_frame = "+" if frame > 0 else "-"

    # push case that we're in reading frame from start onto queue, but
    # only if it's not a stop codon
    codon, orf_pos, query_pos = codons[0]
    if codon not in STOP_CODONS:
        # note what we're adding here: query_pos is query position *in
        # frame*. So even if the ORF is open-ended, we will still not
        # start from the beginning of the sequence, but rather the
        # beginning of the sequence in frame.
        orf_queue.append((orf_pos, query_pos, False))
    
    for codon, orf_pos, query_pos in codons:
        codon = codon.upper()

        if codon in START_CODONS:
            #print "adding start codon '%s' pos %d to queue" % (codon, query_pos)
            orf_queue.append((orf_pos, query_pos, True))
            continue
        if codon in STOP_CODONS:
            # pop everything off queue and make it an ORF to add to
            # the candidates list
            while True:
                try:
                    orf_start_pos, query_start_pos, had_start = orf_queue.popleft()
                except IndexError:
                    break
                orf_data = {"no_start":not had_start, "no_stop":False}
                orf = SeqRange(Range(query_start_pos, query_pos+2), seqname,
                               orf_frame, seqlength=seqlength, data=orf_data)
                all_orfs.append(orf)

    # iteration complete. If there are still items in the ORF queue,
    # pop them off and add them as incomplete.
    if len(orf_queue) > 0:
        while True:
            try:
                orf_start_pos, query_start_pos, had_start = orf_queue.pop()
            except IndexError:
                break
            orf_data = {"no_start":not had_start, "no_stop":True}
            orf = SeqRange(Range(query_start_pos, query_pos+2), seqname,
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
