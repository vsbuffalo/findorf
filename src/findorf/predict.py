"""
predict.py contains the rules and functions that predict ORF.


ContigSequence provides all information that can be gathered from the
ContigSequence; these prediction and annotation rules are those that
can be applied to infer ORFs and annotation.

"""

from utilities import get_all_orfs, get_codons

## Biological constants TODO get from BioPython
STOP_CODONS = set(["TAG", "TGA", "TAA"])
START_CODONS = set(["ATG"])

def orf_with_frameshift(seq, closest_relative_anchors,
                        missing_5prime, key='e'):
    """
    Predict an ORF in the case that we have a frameshift mutation. In
    this case, we can't rule out the possibility that our protein is
    real in the first frame, so we use the 5'-most frame and start
    finding from there. If there's a missing 5'-end we dispatch to
    `get_all_ORFs` as `predict_ORF_missing_5prime` would.

    `key` is what we should use to determine closest relative.
    
    """
    frame_5prime_hsp = closest_relative_anchors.most_5prime.frame

    # missing 5'-end so we're assuming we're already reading.
    all_orfs = get_all_ORFs(seq, frame_5prime_hsp, len(seq), missing_5prime)
    
    return all_orfs
    

def orf_with_missing_5prime(seq, frame):
    """
    Predict an ORF in the case that we have a missing 5'-end.

    We trust relatives' information here, using the 5'-most HSP. Since
    we assume the 5'-end is missing, we just read until we hit a stop
    codon.

    Note that there is an interesting exception here: if a stop codon
    exists because this is a pseudogene. We have to refer to
    information about the HSP's 3'-end for this.
    
    """

    all_orfs = get_all_ORFs(seq, frame, len(seq), in_reading_frame=True)

    return all_orfs


def annotate_orf(anchor_hsps, orf):
    """
    If we have an ORF (from ORF class), annotate some
    obvious characteristiself about it.
    """
    annotation = dict()

    annotation['missing_stop'] = orf.no_stop
    annotation['missing_start'] = orf.no_start
    annotation['full_length'] = not any((orf.no_start, orf.no_stop))

    cisc = contains_internal_stop_codon(anchor_hsps, orf)
    annotation['contains_stop'] = cisc

    return annotation


def orf_vanilla(seq, frame):
    """
    The vanilla case: we have a sequence and a frame, full 5'-end, and
    no frameshift, so we create all ORFs in this frame, as a ribosome
    would.

    If `assuming_missing_start` is True, we start pretending we are in
    an ORF, reading until the first stop, and including this as an
    ORF. This is False, and should not be set to True until (if)
    `get_all_ORFs` is changed so that it restarts the loop when the
    first ORF with a stop codon but no start is found. This is tricky;
    a 5'-UTR could have a stop codon, so we would end up with
    overlapping cases. For now, we take the conservative approach.

    Note that the missing 5' logic will also catch these cases.
    """
    all_orfs = get_all_ORFs(seq, frame, len(seq),
                            in_reading_frame=assuming_missing_start)
    return all_orfs
