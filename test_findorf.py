"""
test_findorf.py contains unit tests for high-impact functions,
i.e. those that have a high value of:

(number of calls) * P(error is silent)

"""
import unittest
from nose.tools import with_setup
from Bio.Seq import Seq
from rules import put_seq_in_frame
from rules import any_overlap, get_codons
from rules import HSP
from rules import get_anchor_HSPs, AnchorHSPs

def test_seq_frame_generator():
    seq = "AAGATGT"
    seq_rc = str(Seq(seq).reverse_complement())
    seq_in_frame = {1: seq, 2: seq[1:], 3:seq[2:],
                    -1: seq_rc, -2: seq_rc[1:], -3:seq_rc[2:]}
    for frame, seq_if in seq_in_frame.iteritems():
        yield assert_put_seq_in_frame, seq, frame, seq_if
    

def assert_put_seq_in_frame(seq, frame, seq_in_frame):
    assert(str(put_seq_in_frame(Seq(seq), frame)) == seq_in_frame)


def test_get_codons():
    seq = Seq("AGTGAGATCTGATCG")

    assert(get_codons(str(seq), 1) == [('AGT', 0, 0), ('GAG', 3, 3),
                                       ('ATC', 6, 6), ('TGA', 9, 9), ('TCG', 12, 12)])
    
    assert(get_codons(str(seq), 2) == [('GTG', 0, 1), ('AGA', 3, 4),
                                       ('TCT', 6, 7), ('GAT', 9, 10)])
    
    assert(get_codons(str(seq), 3) == [('TGA', 0, 2), ('GAT', 3, 5),
                                       ('CTG', 6, 8), ('ATC', 9, 11)])

    seq_rc = seq.reverse_complement()
    assert(get_codons(str(seq_rc), 1) == [('CGA', 0, 0), ('TCA', 3, 3),
                                          ('GAT', 6, 6), ('CTC', 9, 9),
                                          ('ACT', 12, 12)])

    assert(get_codons(str(seq_rc), 2) == [('GAT', 0, 1), ('CAG', 3, 4),
                                          ('ATC', 6, 7), ('TCA', 9, 10)])

    assert(get_codons(str(seq_rc), 3) == [('ATC', 0, 2), ('AGA', 3, 5),
                                          ('TCT', 6, 8), ('CAC', 9, 11)])


def test_any_overlap_generator():
    """
    `any_overlap` must work on both forward and reverse strands.

    These cases were tested in IRanges.
    """

    forward_a = (10, 20)
    forward_b = (20, 24)
    forward_c = (1, 15)
    forward_d = (1, 10)
    forward_e = (1, 3)
    forward_f = (16, 20)

    foward_cases = [(forward_a, forward_b, True, True),
                    (forward_a, forward_b, False, False),
                    (forward_a, forward_c, True, True),
                    (forward_a, forward_e, True, False), 
                    (forward_a, forward_e, False, False), 
                    (forward_f, forward_e, False, False), 
                    (forward_c, forward_d, False, True),
                    (forward_c, forward_d, True, True)]
    # forward cases
    for a, b, closed, truth in foward_cases:
        yield assert_any_overlap, a, b, closed, truth

def assert_any_overlap(a, b, closed, truth):
    assert(any_overlap(a, b, closed=closed) == truth)


def test_get_anchor_HSPs():
    """
    This is a very important function that must work on both forward
    and reverse strands.
    
    """
    # forward strand
    hsp_1 = HSP(e=0.0,identities=132, length=156,
                percent_identity=0.8461, title=u'Bradi2g60290.1',
                query_start=3, query_end=470, sbjct_start=632,
                sbjct_end=787, frame=2)

    hsp_2 = HSP(e=0.0, identities=53, length=89,
                percent_identity=0.5955, title=u'Bradi2g60290.1',
                query_start=1321, query_end=1572, sbjct_start=256,
                sbjct_end=344, frame=2)

    hsp_3 = HSP(e=2.89e-173, identities=182, length=268,
                percent_identity=0.6791, title=u'GRMZM2G303587_T01',
                query_start=488, query_end=1291, sbjct_start=198,
                sbjct_end=465, frame=2)
    hsp_4 = HSP(e=2.891e-173, identities=38, length=58,
                percent_identity=0.6551, title=u'GRMZM2G303587_T01',
                query_start=1324, query_end=1497, sbjct_start=125,
                sbjct_end=182, frame=2)
    
    # reverse strand
    hsp_5 = HSP(e=0.0, identities=216, length=278, percent_identity=0.7769,
                title=u'Bradi2g60290.1', query_start=488, query_end=1321,
                sbjct_start=348, sbjct_end=625, frame=-2)
    hsp_6 = HSP(e=0.0, identities=53, length=89,
                percent_identity=0.5955, title=u'Bradi2g60290.1',
                query_start=1323, query_end=1572, sbjct_start=256,
                sbjct_end=344, frame=-2)

    hsp_7 = HSP(e=2.89e-173, identities=182, length=268,
                percent_identity=0.6791, title=u'GRMZM2G303587_T01',
                query_start=489, query_end=1291, sbjct_start=198,
                sbjct_end=465, frame=-2)
    hsp_8 = HSP(e=2.891e-173, identities=38, length=58,
                percent_identity=0.6551, title=u'GRMZM2G303587_T01',
                query_start=1324, query_end=1497, sbjct_start=125,
                sbjct_end=182, frame=-2)

    
    forward_relatives = {
        'rel_a': [
            hsp_1,
            hsp_3,
            hsp_4],
        'rel_b': [
            hsp_1,
            hsp_2,
            hsp_3,
            hsp_4]
            }

    reverse_relatives = {
        'rel_a': [hsp_5, hsp_6, hsp_7, hsp_8]
        }

    # forward
    assert(get_anchor_HSPs(forward_relatives, False)['rel_a'] == AnchorHSPs(hsp_1, hsp_4, 1))
    assert(get_anchor_HSPs(forward_relatives, False)['rel_b'] == AnchorHSPs(hsp_1, hsp_2, 1))

    # reverse
    assert(get_anchor_HSPs(reverse_relatives, True)['rel_a'] == AnchorHSPs(hsp_6, hsp_5, -1))


