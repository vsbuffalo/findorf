# test_contig.py
# many tests cases based off k36_contig_9886

from nose.tools import assert_equal

from BioRanges.lightweight import SeqRange, Range
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from findorf.orfprediction import ORFTypes
from findorf.contig import Contig

fake_seq = Seq("""GCCAAGCACATTGATGTTCTAAATAAAAGCCCTGTGGTGTCACATAAAAGATCAAGGTATCAACCAAAGACCAGACAGATACTGTGCTGCCTTCAGACAATACAAACGGCAGGCATGCTCAGCGACAAGCATATCTCTCTGCCTCTGTACACTGGATACAGACTGGTCTGTTGTTAACCTTACTTAATTATTCCCGCCATTACAGATCGGCGTCTATGTAGCAGTATAATACAAAGGTCATGCAATATATATAACACTATGATCCATGATCACCTACAACAGCAGCAGTTCTAAGCACTAGGGCAACAAAAAACAAGGGAACCATGCCATGCCACAACATGGAAGGTAGGCAGGCACATCATCACAACAAAAAGAAGCAAGGTGGTATTAGACAGGAAGGAGCTTTACATTTGATGACAGTAATGGCTTGGCTTGGCTTGGCTTACTCATCAGTCCTCTCATTTCATTCTGACTGAACATCGTCAGCTCATGGTGCCAATGCCGCCGCCGCTACCTACCATGGCTTCAGGTCGCCCTTGGGGCCCTGGGTCCAGTTCAGCGCCTTCTGCCCCGGGGTCAGGTACACGCTCTTGCTGTGCGGGAGCTTCCGCGCCAGGCCGTTCAGGATCTGGTGGAAGGCGTCCCCGGGCTGGGCCGGCTGCGGGAAAAGCATGGCCTGCTTCATCGACGGGTTCGCCAGCAGCCTCTGCTTCCGTAGGATCACCTCCATCGGGATGGAGGTCAGGATGGTGTCAAGCTTGGGGACATCGTCCTCCGCCACGAACACGCCGATCTCGTCCCACGGGATGGCGTCCGCGAAGGGCAGGACAATGTCGTCCGCGATGATCACCGGGATGCACCCGAACACCACGGCCTCAACGAGACGGGGGCTCCACGGGGCCCAGCCCAGCGGGCACAGGCAGAAGATGGCACGCTGCATGTCCTCATAGTAGGTGGGTGGGTGGTCCGTCGAGATGTCGAACAGAGGGTTGTTCTTGAAGTTCTCCCACACCGACGCACGGGCGCCTCTTGCATAGTAACCACCCTCGGGATCATTTGCCGTATCGTAGAACAAACCACGGAAATAGACAAAGATGGAGCGCGGGGTCTCTGGGGGTACAAGGTGAGTTTTCATTTTCTGGGGAGGAGCATATGGTGGGATGTTGATTGAGCCCTCCTTTAGGCAGACATGATCCTTCTGCCCAAATGTCTGGACAAGTGTAGCGCGGCGAAGCAATGGAAGGATCCCCCTTTCGATCGCCTTCTCTTCCTGATAATGGAAGCATGCCCCAAAGTCATGTGGCACGACGAAGAAATGATCAGCACCGGCCGTCCGGTTCCAATAGGGCCAGTGCGAGGAGATGAACTGAACCGCACTCCTCATGATCCGCGGGGATTTGAAGGGCAACGGATGACCCCATGGAGTAAGGTCACAGGTCGTGTACACCGGGGTGTAGAACCAGTCGGCCTCCTCCGGGTTCATCGTCCGGATCGCGCTCGACAGCAGGAACCGGTGCAT""")


def create_fake_contig():
    c1 = Contig(SeqRecord(fake_seq, id="fake contig 1", description="fake description"))

    # this example found using findorf
    d = {'frame': -1, 'no_start': True, 'most_5prime_hsp': {'relative':'fake-rel', 'title':'fake-hsp'},
         'no_stop': False}
    orf = SeqRange(Range(0, 1004), seqname="fake contig 1", strand="+", seqlength=len(fake_seq), data=d)
    c1.orf = orf
    c1.orf_type = ORFTypes(orf)
    return orf, c1

def test_contig_orf_seq():
    """
    Check that the correct ORF sequence given contig and ORF SeqRange
    is returned.
    """

    # ORFs are always on forward strand. The ORFs on this one were
    # found on the reverse strand, frame -1. To compare we must first
    # reverse complement.
    orf, c1 = create_fake_contig()
    assert_equal(str(c1.orf_seq.seq), str(fake_seq.reverse_complement()[orf.start:orf.end+1]))


def test_contig_orf_masked():
    """
    Check whether masking is working.
    """
    orf, c1 = create_fake_contig()
    masked_seq = str(c1.orf_masked.seq)
    should_be_masked_seq = "X"*1005 + str(fake_seq.reverse_complement())[1005:]
    assert_equal(masked_seq, should_be_masked_seq)

def test_conistency():
    """
    We should be able to rebuild a sequence from the ORF and masked
    regions - this tests that.

    This builds off our test case (in which the ORF starts from the
    first basepair, so we when reconstructing masked sequence and ORF,
    we just concatenate them)
    """
    orf, c1 = create_fake_contig()
    masked_seq = str(c1.orf_masked.seq)
    orf_seq = str(c1.orf_seq.seq)
    assert_equal(len(orf_seq), masked_seq.count("X"))
    assert_equal(orf_seq+masked_seq.replace("X", ""), str(fake_seq.reverse_complement()))

