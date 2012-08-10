## test_findorf.py -- unit tests for findorf.py

import unittest
from nose.tools import with_setup
from Bio.Seq import Seq
from rules import put_seq_in_frame
from rules import any_overlap

def test_seq_frame_generator():
    seq = "AAGATGT"
    seq_rc = str(Seq(seq).reverse_complement())
    seq_in_frame = {1: seq, 2: seq[1:], 3:seq[2:],
                    -1: seq_rc, -2: seq_rc[1:], -3:seq_rc[2:]}
    for frame, seq_if in seq_in_frame.iteritems():
        yield assert_put_seq_in_frame, seq, frame, seq_if
    

def assert_put_seq_in_frame(seq, frame, seq_in_frame):
    assert(str(put_seq_in_frame(Seq(seq), frame)) == seq_in_frame)

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

    reverse_a = (20, 10)
    reverse_b = (24, 20)
    reverse_c = (15, 1)
    reverse_d = (10, 1)
    reverse_e = (3, 1)
    reverse_f = (30, 14)

    reverse_cases = [(reverse_a, reverse_b, True, True),
                    (reverse_a, reverse_b, False, False),
                    (reverse_a, reverse_c, True, True),
                    (reverse_a, reverse_e, True, False),
                    (reverse_a, reverse_e, False, False),
                    (reverse_f, reverse_e, False, False),
                    (reverse_c, reverse_d, False, True),
                    (reverse_c, reverse_d, True, True)]

    # forward cases
    for a, b, closed, truth in foward_cases:
        yield assert_any_overlap, a, b, False, closed, truth

    # reverse cases
    for a, b, closed,truth in reverse_cases:
        yield assert_any_overlap, a, b, True, closed, truth

def assert_any_overlap(a, b, reverse, closed, truth):
    assert(any_overlap(a, b, reverse, closed=closed) == truth)


