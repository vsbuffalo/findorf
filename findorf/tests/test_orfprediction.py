# test_orfprediction.py

from nose.tools import assert_equal

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

from findorf.orfprediction import get_all_orfs

def test_get_all_orfs_case_generator():
    cases = list()
    cases.append(dict(seq=SeqRecord(Seq("GGGGGGATGGGGGGGATGGGGATGGGGTAGGGGGGGATGGGG"), "case with frame=1"),
                      frame=1, start=[0, 6, 15, 21, 36], end=[27, 27, 27, 27, 41]))

    cases.append(dict(seq=SeqRecord(Seq("CCCCATCCCCCCCTACCCCATCCCCATCCCCCCCATCCCCCCA"), "case with frame=-2"),
                      frame=-2, start=[1, 7, 16, 22, 37], end=[28, 28, 28, 28, 42]))

    # corner cases:
    cases.append(dict(seq=SeqRecord(Seq("ATGTAG"), "case with frame=1, start stop"), frame=1,
                      start=[0], stop=[4]))
    for case in cases:
        yield check_get_all_orfs, case['seq'], case['frame'], case['start'], case['end']

def check_get_all_orfs(seq, frame, start, end):
    orfs = get_all_orfs(seq, frame)
    assert_equal(sorted(orfs.start), sorted(start))
    assert_equal(sorted(orfs.end), sorted(end))
