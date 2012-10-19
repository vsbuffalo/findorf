# test_orfprediction.py

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

import findorf.orfprediction

TEST_FASTA = "cases.fasta"

def sequence_setup():
    """
    Setup some sequence fixtures using the cases.fasta file, which has
    IDs and fame frame in description.
    """
    case_1 = dict(seq="GGGGGGATGGGGGGGATGGGGATGGGGTAGGGGGGGATGGGG", frame=1,
                  start=[0, 6, 15, 21, 36], end=[27, 27, 27, 27, 41])
    


@withsetup(sequence_setup)
def test_get_all_orfs():
    orfs_1 = orfprediction.get_all_orfs(case_1['seq'], case_1['frame'])
    assert(sorted(orfs_1.start) == sorted(case_1['start']))
    
