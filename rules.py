"""
Rules.py contains the rules that predict ORF. These are not class
methods in the ContigSequence class, even though they take these
objects because these should be a bit more modular.

There are also functions for annotating traits dependent upon the ORF,
i.e. if there's a pseudogene due to internal stop codon.
"""


try:
    from Bio.Seq import Seq
except ImportError, e:
    sys.exit("Cannot import BioPython modules; please install it.")
from collections import namedtuple
import pdb

## Biological constants
STOP_CODONS = set(["TAG", "TGA", "TAA"])
START_CODONS = set(["ATG"])

## Predefined tuple structures
orf_fields = ['start', 'stop', 'orf']
PredictedORF = namedtuple("PredictedORF", orf_fields)


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

def contains_internal_stop_codon(cs):
    """
    """
    pass

def predict_ORF_frameshift(cs):
    """
    Predict an ORF in the case that we have a frameshift mutation. In
    this case, we can't rule out the possibility that our protein is
    real in the first frame, so we use the 5'-most frame and start
    finding from there.
    
    """
    pass

def predict_ORF_missing_5prime(cs):
    """
    Predict an ORF in the case that we have a missing 5'-end.

    """
    pass

def predict_ORF_vanilla(cs):
    """
    

    """
    seq = cs.seq
    frame = cs.majority_frame
    
    seq_in_frame = str(put_seq_in_frame(seq, frame))
    codons = get_codons(seq_in_frame)

    all_orfs = list()
    start_pos = None
    position = None    
    in_reading_frame = False
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

    # first ORF is 5'-most
    if not len(all_orfs):
        return PredictedORF(None, None, None)
    
    best_start, best_stop = all_orfs[0]
    return PredictedORF(best_start, best_stop, None)

