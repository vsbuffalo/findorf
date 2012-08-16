"""
utilities.py contains utilities that handle creating codons, puting
sequences in frame, etc.

"""

from Bio.Seq import Seq

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

