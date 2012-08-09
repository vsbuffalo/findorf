"""
Rules.py contains the rules that predict ORF. These are not class
methods in the ContigSequence class, even though they take these
objects because these should be a bit more modular.

"""


try:
    from Bio.Seq import Seq
except ImportError, e:
    sys.exit("Cannot import BioPython modules; please install it.")


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


    """

def predict_ORF(cs):
    """
    

    """
    seq = cs.seq
    frame = cs.majority_frame
    
    seq_in_frame = put_seq_in_frame(seq, frame)
    codons = get_codons(seq_in_frame)
    
