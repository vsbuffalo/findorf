"""
Rules.py contains the rules and functions that predict ORF.


ContigSequence provides all information that can be gathered from the
ContigSequence; these prediction and annotation rules are those that
can be applied to infer ORFs and annotation.

"""


try:
    from Bio.Seq import Seq
except ImportError, e:
    sys.exit("Cannot import BioPython modules; please install it.")
from collections import namedtuple
import pdb
from operator import attrgetter, itemgetter
import ContigSequence 

## Biological constants
STOP_CODONS = set(["TAG", "TGA", "TAA"])
START_CODONS = set(["ATG"])

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


def get_codons(seq, frame):
    """
    Return a list of (codon, position_in_orf, pos_in_forward_query) tuples.
    """

    if frame < 0:
        try:
            seq = seq.reverse_complement()
        except AttributeError:
            msg = "no reverse_complement() method; 'seq' must of class Bio.Seq.Seq"
            raise AttributeError, msg
        frame = abs(frame)

    tmp = [(str(seq[pos:pos+3]), pos-(frame-1), pos) for
           pos in range(frame-1, len(seq), 3)]

    # remove the last string if not a full codon.
    return [(c, po, pfq) for c, po, pfq in tmp if len(c) == 3]


def any_overlap(range_a, range_b, closed=True):
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

    assert(a_start <= a_end)
    assert(b_start <= b_end)

    if closed:
        return b_start <= a_end and a_start <= b_end
    if not closed:
        return b_start < a_end and a_start < b_end

def get_ORF_overlaps_5prime_HSP(orf_list, anchor_hsps, allow_one_side=True):
    """
    Return True if start_pos and stop_pos have *any* overlap with the
    5'-most anchor HSP of the closest relative.

    Note that this assumes start and stop are forward strand; we can
    safely assume this because orf prediction is done this way (via
    `put_seq_in_frame`.

    `allow_one_side` indicates whether we should allow one side
    overlap in the case of missing start/stop positions (i.e. those
    are None).
    """
    # we have to specify (via attribute) that we want the most_5prime,
    # as `get_closest_relative_anchor_HSP` will return an AnchorHSP
    # object, even with `which` specified, since `which` just
    # indicates which to sort by, not to return.
    m5p_hsp = get_closest_relative_anchor_HSP(anchor_hsps).most_5prime

    for orf in orf_list:
        if orf is None:
            continue
        qstart, qend = orf.query_start, orf.query_end
        if None in (qstart , qend):
            if not allow_one_side:
                continue
            # extend the HSP as far towards to the 5' or 3' end as so
            # to not prevent overlap.
            if qstart is None:
                qstart = 0
            if qend is None:
                qend = float('inf')            
        has_overlap = any_overlap((qstart, qend),
                                  (m5p_hsp.query_start, m5p_hsp.query_end))

        if has_overlap:
            return orf
    return None

def anchor_hsp_attrgetter(which=None, key='e'):
    """
    Like operator.attrgetter, but for AnchorHSPs object. Dense
    data structures need functions for digging into them cleanly.

    This may look strange, but recall attrgetter returns a function,
    so this is essentially currying.

    `x` are the key value tuples from a dictionary of AnchorHSPs.
    """
    
    return lambda x: attrgetter(key)(attrgetter(which)(x[1]))


def get_anchor_HSPs(relative_dict, is_reversed):
    strand = -1 if is_reversed else 1
    anchor_hsps = dict()
        
    for relative, hsps in relative_dict.items():
        # hsp_1 corresponds to the HSP with the latest end position
        hsp_1 = sorted(hsps, key=attrgetter('query_end'), reverse=True)[0]

        # hsp_2 corresponds to the HSP with the earliest start position
        hsp_2 = sorted(hsps, key=attrgetter('query_start'))[0]
        
        if is_reversed:
            # if reversed, the HSP closest to the protein N-terminus
            # is the one with the latest end position
            anchor_hsps[relative] = ContigSequence.AnchorHSPs(hsp_1, hsp_2, strand)
        else:
            # on the forward strand, the opposite is true.
            anchor_hsps[relative] = ContigSequence.AnchorHSPs(hsp_2, hsp_1, strand)

    return anchor_hsps


def get_closest_relative_anchor_HSP(anchor_hsps, which='most_5prime', key='e'):
    """
    When looking for anchor HSPs (for use with finding internal stop
    codons and frameshifts), we will concern ourselves with the HSPs
    of the closest relative. This is a function to get the closet
    relative's AnchorHSPs, sorting by `key` (usually, idenitites,
    percent_identity, or e). `which` revers which anchor HSP to look
    at: 5' or 3'.
    """
    if not len(anchor_hsps):
        return None
    
    tmp = sorted(anchor_hsps.iteritems(),
                                  key=anchor_hsp_attrgetter(which, key))

    cr_anchor_hsps = tmp[0][1]

    return cr_anchor_hsps

def contains_internal_stop_codon(anchor_hsps, predicted_orf):
    """
    After we make an ORF prediction (the 5'-most ORF), we want to
    check that we don't have an HSP more 3' than our ORF stop
    position, which would indicate a disrupted protein (which could
    still be functional or a pseudogene).

    We look at all relative's anchor HSPs for maximum sensitivity, in
    case the closest relative's proteins are misannotated.

    TODO: unit test
    """

    most_3prime = None
    for relative, ahsp in anchor_hsps.iteritems():
        this_hsp = ahsp.most_3prime
        if this_hsp.query_start > predicted_orf.query_end:
            return True
    return False


def get_all_ORFs(codons, frame, in_reading_frame=False):
    """
    Generic ORF finder; it returns a list of all ORFs as they are
    found, given codons (a list if tuples in the form (codon,
    position)) from `get_codons`.
    
    """
    
    all_orfs = list()
    orf_start_pos = None
    query_start_pos = None
    first_stop_hit = False
    # Note that query_pos is forward strand.
    for codon, orf_pos, query_pos in codons:
        if codon in START_CODONS and not in_reading_frame:
            in_reading_frame = True
            orf_start_pos = orf_pos
            query_start_pos = query_pos
            continue
        if in_reading_frame and codon in STOP_CODONS:
            # a full reading frame, unless we haven't hit any stop
            # yet, then we're in the possible ORF from the start of
            # the query to the end.
            all_orfs.append(ContigSequence.ORF(orf_start_pos, orf_pos, query_start_pos,
                                               query_pos, frame,
                                               no_start=(not first_stop_hit)))
            in_reading_frame = False
            orf_start_pos = None
            query_start_pos = None
            first_stop_hit = True
            
            continue
        
    if in_reading_frame:
        all_orfs.append(ContigSequence.ORF(orf_start_pos, orf_pos,
                                           query_start_pos, query_pos,
                                           frame, no_start=False,
                                           no_stop=True))

    return all_orfs


def predict_ORF_frameshift(seq, anchor_hsps, missing_5prime, key='e'):
    """
    Predict an ORF in the case that we have a frameshift mutation. In
    this case, we can't rule out the possibility that our protein is
    real in the first frame, so we use the 5'-most frame and start
    finding from there. If there's a missing 5'-end we dispatch to
    `get_all_ORFs` as `predict_ORF_missing_5prime` would.

    `key` is what we should use to determine closest relative.
    
    """
    # Get the frame of the closest relative's (by number of
    # identities) most 5'-end. Note that we look at relative closeness
    # by lowest e-value (by default `key='e'`) of which ever anchor
    # HSP we are talking about (in this case, 5'-most).
    closest_relative_anchors = get_closest_relative_anchor_HSP(anchor_hsps, key=key)
    frame_5prime_hsp = closest_relative_anchors.most_5prime.frame

    codons = get_codons(seq, frame_5prime_hsp)

    # missing 5'-end so we're assuming we're already reading.
    all_orfs = get_all_ORFs(codons, frame_5prime_hsp, missing_5prime)
    
    if not len(all_orfs):
        return None

    return all_orfs
    

def predict_ORF_missing_5prime(seq, frame):
    """
    Predict an ORF in the case that we have a missing 5'-end.

    We trust relatives' information here, using the 5'-most HSP. Since
    we assume the 5'-end is missing, we just read until we hit a stop
    codon.

    Note that there is an interesting exception here: if a stop codon
    exists because this is a pseudogene. We have to refer to
    information about the HSP's 3'-end for this.
    
    """

    codons = get_codons(seq, frame)

    all_orfs = get_all_ORFs(codons, frame, in_reading_frame=True)

    if not len(all_orfs):
        return None

    return all_orfs


def annotate_ORF(anchor_hsps, orf):
    """
    If we have an ORF (from ORF class), annotate some
    obvious characteristiself about it.
    """
    annotation = dict()

    annotation['missing_stop'] = orf.no_stop
    annotation['missing_start'] = orf.no_start
    annotation['full_length'] = None not in (orf.start, orf.end)

    if annotation['full_length']:
        cisc = contains_internal_stop_codon(anchor_hsps, orf)
        annotation['contains_stop'] = cisc
    return annotation


def predict_ORF_vanilla(seq, frame, assuming_missing_start=True):
    """
    The vanilla case: we have a sequence and a frame, full 5'-end, and
    no frameshift, so we create all ORFs in this frame, as a ribosome
    would.

    If `assuming_missing_start` is True, we start pretending we are in
    an ORF, reading until the first stop, and including this as an
    ORF.

    """
    codons = get_codons(seq, frame)

    all_orfs = get_all_ORFs(codons, frame, in_reading_frame=assuming_missing_start)
    
    # first ORF is 5'-most
    if not len(all_orfs):
        return None
    
    return all_orfs
