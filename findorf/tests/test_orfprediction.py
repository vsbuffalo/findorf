# test_orfprediction.py

from nose.tools import assert_equal

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

from findorf.orfprediction import get_all_orfs, count_5prime_ATG

def test_get_all_orfs_case_generator():
    cases = list()
    cases.append(dict(seq=SeqRecord(Seq("GGGGGGATGGGGGGGATGGGGATGGGGTAGGGGGGGATGGGG"), "case with frame=1"),
                      frame=1, start=[0, 6, 15, 21, 36], end=[29, 29, 29, 29, 41]))

    cases.append(dict(seq=SeqRecord(Seq("CCCCATCCCCCCCTACCCCATCCCCATCCCCCCCATCCCCCCA"), "case with frame=-2"),
                      frame=-2, start=[1, 7, 16, 22, 37], end=[30, 30, 30, 30, 42]))

    ## corner cases
    #
    # Start and stop codons next to each other. Here, we
    # could have two possible ORFs:
    #  1. open, containing unknown sequence, then an ATG, then stop
    #  2. not open, starting with first ATG, then stopping
    cases.append(dict(seq=SeqRecord(Seq("ATGTAG"), "case with frame=1, start stop"), frame=1,
                      start=[0, 0], end=[5, 5]))
    # if all we have is a stop, we won't even assume we have an open ORF
    cases.append(dict(seq=SeqRecord(Seq("TAG"), "case with frame=1, just stop"), frame=1,
                      start=[], end=[]))
    cases.append(dict(seq=SeqRecord(Seq("GCGGCGTAGGCGATGGCG"), "case with frame=1, mix"), frame=1,
                      start=[0, 12], end=[8, 17]))
    cases.append(dict(seq=SeqRecord(Seq("AAGCGGCGTAGGCGATGGCG"), "case with frame=3, mix"), frame=3,
                      start=[2, 14], end=[10, 19]))
    cases.append(dict(seq=SeqRecord(Seq("ATGGCGGCGTAGGCGATGGCG"), "case with frame=1, duplicate ranges"),
                      frame=1, start=[0, 0, 15], end=[11, 11, 20]))

    for case in cases:
        yield check_get_all_orfs, case['seq'], case['frame'], case['start'], case['end']

def check_get_all_orfs(seq, frame, start, end):
    orfs = get_all_orfs(seq, frame)
    assert_equal(orfs.start, start)
    assert_equal(orfs.end, end)

def test_count_5prime_ATG():
    s1 = Seq("ATGGGGATGGGGATG") # frame 1
    s2 = Seq("AATGGGGATGGGGATG").reverse_complement() # frame -2
    
    assert_equal(count_5prime_ATG(s1, 1, 0), 0)
    assert_equal(count_5prime_ATG(s1, 1, 100), 3)
    assert_equal(count_5prime_ATG(s1, 1, 4), 1)
    assert_equal(count_5prime_ATG(s1, 1, 8), 2)
    assert_equal(count_5prime_ATG(s1, 1, 6), 1)

    # same cases above, with frame adjusted. 
    assert_equal(count_5prime_ATG(s2, -2, 0), 0)
    assert_equal(count_5prime_ATG(s2, -2, 100), 3)
    assert_equal(count_5prime_ATG(s2, -2, 4), 1)
    assert_equal(count_5prime_ATG(s2, -2, 8), 2)
    assert_equal(count_5prime_ATG(s2, -2, 6), 1)
    
