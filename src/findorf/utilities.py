"""
utilities.py contains utilities that handle creating codons, puting
sequences in frame, etc.

"""

from Bio.Seq import Seq
from RangedFeatures import ORF
from collections import Counter

## Biological constants TODO get from BioPython
STOP_CODONS = set(["TAG", "TGA", "TAA"])
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


def get_all_orfs(seq, frame, query_length,
                 in_reading_frame=False, add_partial=False):
    """
    Generic ORF finder; it returns a list of all ORFs as they are
    found, given codons (a list if tuples in the form (codon,
    position)) from `get_codons`.
    """
    
    all_orfs = list()
    orf_start_pos = None
    query_start_pos = None

    codons = get_codons(seq, frame)
    
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
            all_orfs.append(ORF(query_start_pos, query_pos, query_length,
                                frame, no_start=query_start_pos is None,
                                no_stop=False))
            
            in_reading_frame = False
            orf_start_pos = None
            query_start_pos = None
            
    # add an partial ORFs, and mark as having no stop codon
    if in_reading_frame:
        all_orfs.append(ORF(query_start_pos, query_pos, query_length,
                            frame, no_start=query_start_pos is None,
                            no_stop=True))
    if add_partial:
        for codon, orf_pos, query_pos in codons:
            if codon in STOP_CODONS:
                all_orfs.insert(0, ORF(None, query_pos, query_length,
                                       frame, no_start=True,
                                       no_stop=True))
                break

    return all_orfs

def summarize_contigs(contigs):
    """
    Summarize annotation.
    """

    terms = ['majority_frameshift', 'missing_5prime', 'hsp_orf_overlap',
             'inconsistent_strand', 'has_orf', 'has_relatives',
             'missing_start', 'missing_stop', 'internal_stop',
             'num_relatives', 'num_orf_candidates', 'closest_relative']

    annotation_summary = dict([(t, Counter()) for t in terms])

    if type(contigs) is dict:
        contigs = contigs.values()
    
    for contig in contigs:
        all_anno = contig.get_annotation()
        for key, value in all_anno.items():
            annotation_summary[key][value] += 1

    return annotation_summary


def put_seq_in_frame(seq, frame):
    """
    Take a sequence (of Bio.Seq class) and transform it to into the
    correct frame.
    """
    if seq.__class__.__name__ != "Seq":
        seq = Seq(seq)
    if frame < 0:
        seq = seq.reverse_complement()
        frame = -1*frame

    if not frame in range(1, 4):
        raise Exception, "improper frame: frame must be in [1, 3]"
    return seq[(frame-1):]
