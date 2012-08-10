"""
Rules.py contains the rules that predict ORF. These are not class
methods in the ContigSequence class, even though they take these
objects because these should be a bit more modular. ContigSequence
provides all information that can be gathered from the ContigSequence;
these prediction and annotation rules are those that can be applied to
infer ORFs and annotation.


TODO: don't just grab the 5'-most HSP, grab the 5'-most that overlaps
anchor HSPs of closet relatives.
"""


try:
    from Bio.Seq import Seq
except ImportError, e:
    sys.exit("Cannot import BioPython modules; please install it.")
from collections import namedtuple
import pdb
from operator import attrgetter, itemgetter

## Biological constants
STOP_CODONS = set(["TAG", "TGA", "TAA"])
START_CODONS = set(["ATG"])

## Predefined tuple structures
orf_fields = ['start', 'stop', 'frame']
PredictedORF = namedtuple("PredictedORF", orf_fields)

HSP = namedtuple('HSP', ['e', 'identities', 'length',
                         'percent_identity', 'title',
                         'query_start', 'query_end',
                         'sbjct_start', 'sbjct_end',
                         'frame'])


def get_codons(seq_in_frame):
    """
    Return a list of (codon, position) tuples.
    """
    seq = seq_in_frame
    return [(seq[pos:pos+3], pos) for pos in range(0, len(seq), 3)]

def put_seq_in_frame(seq, frame):
    """
    Take a sequence (of Bio.Seq class) and transform it to into the
    correct frame.
    """
    if frame < 0:
        try:
            seq = seq.reverse_complement()
        except AttributeError:
            msg = "no reverse_complement() method; 'seq' must of class Bio.Seq.Seq"
            raise AttributeError, msg
        frame = -1*frame
    if not frame in range(1, 4):
        raise Exception, "improper frame: frame must be in [1, 3]"
    return seq[(frame-1):]

def any_overlap(range_a, range_b, reverse=False, closed=True):
    """
    Is there any overlap in two ranges? We could use interval trees,
    but in our application, we'll only be considering two ranges.

    forward strand:
    
          |---------------| a
                  |------------| b
    
         |---------------| a
                         |----| b

    or              
          |------------| a
    |----------| b

    on the forward strand, always: start < end.
    
    """
    a_start, a_end = range_a
    b_start, b_end = range_b

    if not reverse:
        if closed:
            return b_start <= a_end and a_start <= b_end
        if not closed:
            return b_start < a_end and a_start < b_end
    else:
        if closed:
            return b_start >= a_end and a_start >= b_end
        if not closed:
            return b_start > a_end and a_start > b_end



def overlaps_anchor_HSPs(start_pos, stop_pos, anchor_hsps, query_length):
    """
    Return True if start_pos and stop_pos have *any* overlap with the
    5'-most anchor HSP of the closest relative.
    """
    most_5prime_hsp = closest_relative_anchors(anchor_hsps)

    if anchor_hsp.frame < 0:
        anchor_hsp = put_HSP_on_foward_strand(most_5prime_hsp, query_length)

    return 



def put_HSP_on_foward_strand(hsp, query_length):
    """
    Take an HSP and return it on the forward strand.

    
    If a blastx HSP is reverse complemented to the subject sequence,
    the 5'-most amino acid is the query end and the 3'-most amino acid
    is the query start.

    TODO: unit test this and double check that query_length is correct
    units (bp vs aa).
    """

    if hsp.frame > 0:
        # already on forward strand
        return hsp

    hsp = HSP(e=hsp.e,
              identities=hsp.identities,
              length=hsp.length,
              percent_identity=hsp.percent_identity,
              title=hsp.title,
              query_start=query_length - hsp.query_end + 1,
              query_end=query_length - hsp.query_start + 1,
              sbjct_start=hsp.sbjct_start,
              sbjct_end=hsp.sbjct_end,
              frame=hsp.frame)

    return hsp

def anchor_hsp_attrgetter(which=None, key='e'):
    """
    Like operator.attrgetter, but for AnchorHSPs named tuples. Dense
    data structures need functions for digging into them cleanly.

    This may look strange, but recall attrgetter returns a function,
    so this is essentially currying.

    `x` are the key value tuples from a dictionary of AnchorHSPs.
    """
    
    return lambda x: attrgetter(key)(attrgetter(which)(x[1]))

def get_closest_relative_anchor_HSP(anchor_hsps, which='most_5prime', key='e'):
    """
    When looking for anchor HSPs (for use with finding internal stop
    codons and frameshifts), we will concern ourselves with the HSPs
    of the closest relative. This is a function to get the closet
    relative's AnchorHSPs, sorting by `key` (usually, idenitites,
    percent_identity, or e). `which` revers which anchor HSP to look
    at: 5' or 3'.
    """
    tmp = sorted(anchor_hsps.iteritems(),
                                  key=anchor_hsp_attrgetter(which, key))

    cr_anchor_hsps = tmp[0][1]

    return cr_anchor_hsps

def contains_internal_stop_codon(cs, e_value, pi_range, predicted_orf):
    """
    After we make an ORF prediction (the 5'-most ORF), we want to
    check that we don't have an HSP more 3' than our ORF stop
    position, which would indicate a disrupted protein (which could
    still be functional or a pseudogene).

    We look at all relative's anchor HSPs for maximum sensitivity, in
    case the closest relative's proteins are misannotated.

    TODO: unit test
    """

    anchor_hsps = cs.get_anchor_HSPs(e_value, pi_range)

    most_3prime = None
    for relative, ahsp in anchor_hsps.iteritems():
        this_hsp = ahsp.most_3prime
        if this_hsp.frame < 0:
            this_hsp = put_HSP_on_foward_strand(this_hsp, cs.len)
        if this_hsp.query_start > predicted_orf.stop:
            return True
    return False


def get_all_ORFs(codons, in_reading_frame=False):
    """
    Generic ORF finder; it returns a list of all ORFs as they are
    found, given codons (a list if tuples in the form (codon,
    position)) from `get_codons`.
    
    """
    
    all_orfs = list()
    start_pos = None
    for codon, position in codons:
        if codon in START_CODONS and not in_reading_frame:
            in_reading_frame = True
            start_pos = position
            continue
        if in_reading_frame and codon in STOP_CODONS:
            all_orfs.append((start_pos, position))
            in_reading_frame = False
            continue
    if in_reading_frame:
        all_orfs.append((start_pos, position))

    return all_orfs


def predict_ORF_frameshift(cs, e_value=None, pi_range=None, key='e'):
    """
    Predict an ORF in the case that we have a frameshift mutation. In
    this case, we can't rule out the possibility that our protein is
    real in the first frame, so we use the 5'-most frame and start
    finding from there. If there's a missing 5'-end we dispatch to
    `get_all_ORFs` as `predict_ORF_missing_5prime` would.

    `key` is what we should use to determine closest relative.
    
    """
    anchor_hsps = cs.get_anchor_HSPs(e_value, pi_range)

    
    # Get the frame of the closest relative's (by number of
    # identities) most 5'-end. Note that we look at relative closeness
    # by lowest e-value (by default `key='e'`) of which ever anchor
    # HSP we are talking about (in this case, 5'-most).
    closest_relative_anchors = get_closest_relative_anchor_HSP(anchor_hsps, key=key)
    frame_5prime_hsp = closest_relative_anchors.most_5prime.frame

    seq = cs.seq
    seq_in_frame = str(put_seq_in_frame(seq, frame_5prime_hsp))
    codons = get_codons(seq_in_frame)

    missing_5prime = cs.missing_5prime(anchor_hsps)

    # missing 5'-end so we're assuming we're already reading.
    all_orfs = get_all_ORFs(codons, missing_5prime)
    
    if not len(all_orfs):
        return PredictedORF(None, None, None)

    best_start, best_stop = all_orfs[0]
    return PredictedORF(best_start, best_stop, frame_5prime_hsp)

    

def predict_ORF_missing_5prime(cs):
    """
    Predict an ORF in the case that we have a missing 5'-end.

    We trust relatives' information here, using the 5'-most HSP. Since
    we assume the 5'-end is missing, we just read until we hit a stop
    codon.

    Note that there is an interesting exception here: if a stop codon
    exists because this is a pseudogene. We have to refer to
    information about the HSP's 3'-end for this.
    
    """

    seq = cs.seq
    frame = cs.majority_frame

    seq_in_frame = str(put_seq_in_frame(seq, frame))
    codons = get_codons(seq_in_frame)

    all_orfs = get_all_ORFs(codons, in_reading_frame=True)

    if not len(all_orfs):
        return PredictedORF(None, None, None)

    best_start, best_stop = all_orfs[0]
    return PredictedORF(best_start, best_stop, frame)    


def predict_ORF_vanilla(cs):
    """
    

    """
    seq = cs.seq
    frame = cs.majority_frame
    
    seq_in_frame = str(put_seq_in_frame(seq, frame))
    codons = get_codons(seq_in_frame)

    all_orfs = get_all_ORFs(codons, in_reading_frame=False)
    
    # first ORF is 5'-most
    if not len(all_orfs):
        return PredictedORF(None, None, None)
    
    best_start, best_stop = all_orfs[0]
    return PredictedORF(best_start, best_stop, frame)

def annotate_ORF(cs, e_value, pi_range, orf):
    """
    If we have an ORF (from PredictedORF named tuple), annotate some
    obvious characteristics about it.
    """
    annotation = dict()

    annotation['missing_start'] = orf.start is None
    annotation['missing_stop'] = orf.stop is None
    annotation['full_length'] = None not in (orf.start, orf.stop)

    if annotation['full_length']:
        cisc = contains_internal_stop_codon(cs, e_value, pi_range,  orf)
        annotation['contains_stop'] = cisc

    return annotation


def generic_predict_ORF(cs, e_value=None, pi_range=None):
    """
    The central dispatcher; logic function.

    In the future, it might be nice to do some sort of CLOS-style
    generic method dispatching.

    """

    if cs.has_relatives:
        # we have relatives; we can predict the ORF
        if cs.majority_frameshift:
            orf = predict_ORF_frameshift(cs, e_value, pi_range)
        elif cs.missing_5prime(cs.get_anchor_HSPs(e_value, pi_range)):
            orf = predict_ORF_missing_5prime(cs)
        else:        
            orf = predict_ORF_vanilla(cs)

        # with an ORF, we annotate it    
        orf_annotation = annotate_ORF(cs, e_value, pi_range, orf)
    else:
        orf = None
        orf_annotation = dict()
        
    cs.add_orf_prediction(orf)
    cs.add_annotation(orf_annotation)
    
