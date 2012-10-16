"""
contig.py contains the Contig class. This new version is a depature
from the previous version in that (1) it relies upon BioRanges and (2)
it stores uses objects from BioPython directly.

"""

from BioRanges.lightweight import Range, SeqRange, SeqRanges
from collections import namedtuple

HSP = namedtuple("HSP", ['expect', 'identities', 'align_length',
                         'percent_identity', 'sbjct_start',
                         'sbjct_end', 'frame'])

def _HSP_object_to_Named_Tuple(hsp):
    """
    Convert an HSP object from BioPython's HSP class to a named
    tuple. This is a lightweight format, since we don't need
    everything from the full class. We don't store query start and end
    positions because this will be stored in the SeqRange object.

    We also do some sanity checking here.
    """
    # the BioPython parser doesn't give us a non-zero second
    # frame (which is for use with non-blastx parsers).
    assert(hsp.frame[1] is 0)
    
    # blastx has protein subjects, so this should always be the case
    assert(hsp.sbjct_start < hsp.sbjct_end)

    percent_identity = hsp.identities/float(hsp.align_length)

    hsp = HSP(expect=hsp.expect,
              identities=hsp.identities,
              align_length=hsp.align_length,
              percent_identity=percent_identity,
              sbjct_start=hsp.sbjct_start,
              sbjct_end=hsp.sbjct_end,
              frame=hsp.frame[0])
    return hsp

class Contig():
    """
    Contig represents a contig from the assembly, and has attributes
    and methods to add more information or make predictions about this
    conitg.

    """

    def __init__(self, record):
        self.record = record
        self.annotation = dict()
        self.hsps = SeqRanges()
        
    @property
    def seq(self):
        """
        Return the sequence of the contig.
        """
        return self.record.seq

    @property
    def id(self):
        """
        Return the sequence header ID.
        """
        return self.record.id

    @property
    def description(self):
        """
        Return the sequence header description.
        """

        return self.record.description

    def add_alignment(self, relative, blast_record):
        """
        Add a BLASTX alignment from a relative.
        """
        if len(blast_record.alignments) == 0:
            # no alignments, so we dont have any info to add for this
            # relative.
            return 

        best_alignment = blast_record.alignments[0]
        for hsp in best_alignment.hsps:

            # Negative strand alignments have query_start >
            # query_end. Ranges are require start < end.
            qstart = hsp.query_start
            qend = hsp.query_end
            strand = "+" if hsp.frame[0] < 0 else "+"
            assert(qstart < qend)

            nthsp = _HSP_object_to_Named_Tuple(hsp)
            seqrng = SeqRange(Range(qstart, qend),
                              seqname=self.record.id,
                              strand=strand,
                              data={'relative':relative,
                                    'title':best_alignment.title,
                                    'hsp':nthsp})
            self.hsps.append(seqrng)
            
    def most_5prime_hsp(self):
        """
        Get the 5'-most HSP or PFAM domain, handling the possibility
        the contig is in the reverse orientation.
        """
        pass

    def majority_frame(self):
        """
        Get the majority frame by looking at the ranges on the
        contig. We could also break this down by relative.
        """
        pass

    def majority_frame(self):
        pass
