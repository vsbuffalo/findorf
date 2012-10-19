# test_orfprediction.py

from nose.tools import assert_equal

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

from findorf.orfprediction import get_all_orfs

def test_get_all_orfs_case_generator():
    case_1 = dict(seq=SeqRecord(Seq("GGGGGGATGGGGGGGATGGGGATGGGGTAGGGGGGGATGGGG"), "case_1"),
              frame=1, start=[0, 6, 15, 21, 36], end=[27, 27, 27, 27, 41])

    case_2 = dict(seq=SeqRecord(Seq("CCCCATCCCCCCCTACCCCATCCCCATCCCCCCCATCCCCCCA"), "case_2"),
                  frame=-2, start=[1, 7, 16, 22, 37], end=[28, 28, 28, 28, 42])
    cases = [case_1, case_2]
    for case in cases:
        yield check_get_all_orfs, case['seq'], case['frame'], case['start'], case['end']

def check_get_all_orfs(seq, frame, start, end):
    orfs = get_all_orfs(seq, frame)
    assert_equal(sorted(orfs.start), sorted(start))
    assert_equal(sorted(orfs.end), sorted(end))
