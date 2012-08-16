"""
Contig

# BLASTX Output Notes

Alignments with a negative frame (that is the reverse complement of
the contig maps to the subject protein), have a `query_start` that
corresponds to the *end* of the subject protein, and the `query_end`
is the 5'-most to the subject protein. So for frame < 0, `query_end`
is the 5'-most.

The HSP class from BioPython is documented here:
http://biopython.org/DIST/docs/api/Bio.Blast.Record.HSP-class.html

"""

import sys
import pdb
from collections import Counter, defaultdict
from string import Template
from operator import itemgetter, attrgetter
from math import floor

try:
    from Bio.Blast import NCBIXML
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError, e:
    sys.exit("Cannot import BioPython modules; please install it.")

from RangedFeatures import HSP, AnchorHSPs, RelativeHSPs

GTF_FIELDS = ("seqname", "source", "feature", "start",
              "end", "score", "strand", "frame", "group")        

class Contig():
    """
    Contig represents a contig from the assembly, and has attributes
    and methods to add more information or make predictions about this
    conitg.

    """

    def __init__(self, record):
        """
        Initialize a Contig. The record.id must correspond to the same
        one used in the blastx results.

        Any ORF predictions will be stored as SeqFeatures on this
        Contig.

        We don't subclass SeqRecord here because we don't use their
        SeqFeatures with it - they are too weak in terms of look at
        overlap. It's easier to store particular features (anchorHSPs,
        ORF candidates) as dict or list attributes here.
        """
        # core data attributes
        self.id = record.id
        self.seq = record.seq

        # information added by blastx results
        self.relative_hsps = RelativeHSPs()

        # Annotation attributes
        self.orf = None ## This will be a BioPython SeqFeature
        self.all_orfs = None
        self.annotation = dict()

    def __len__(self):
        """
        Return the contig length.
        """
        return len(self.seq)
        
    def add_relative_alignment(self, relative, blast_record):
        """
        Given a relative and a BioPython BLAST alignment objects,
        extract and store the relevant parts of the _best_ alignment
        only.
        """
        if len(blast_record.alignments) == 0:
            # no alignments, so we dont have any info to add for this
            # relative.
            return 

        # TODO check: are these guaranteed in best first order?
        best_alignment = blast_record.alignments[0]
        for hsp in best_alignment.hsps:
            percent_identity = hsp.identities/float(hsp.align_length)

            # the BioPython parser doesn't give us a non-zero second
            # frame (which is for use with non-blastx parsers).
            assert(hsp.frame[1] is 0)

            # blastx has protein subjects, so this should always be the case
            assert(hsp.sbjct_start < hsp.sbjct_end)
            
            hsp = HSP(e=hsp.expect,
                      identities=hsp.identities,
                      length=hsp.align_length,
                      percent_identity=percent_identity,
                      title=best_alignment.title,
                      query_start=hsp.query_start,
                      query_end=hsp.query_end,
                      sbjct_start=hsp.sbjct_start,
                      sbjct_end=hsp.sbjct_end,
                      frame=hsp.frame[0])

            self.relative_hsps.add_relative_hsp(relative, hsp)
        
        
    def has_relatives(self):
        return len(self.relative_hsps) > 0

    def predict_orf(self, e_value=None, pi_range=None):
        """
        Predict an ORF, given some parameters.

        1. Frameshift? If so, check 5' missing,

        2. 5' missing?

        3. Vanilla
        """

        # filter by the parameters; we use these filtered cases from
        # now on.
        rel_hsps = self.relative_hsps.get_relatives(e_value, pi_range)
        self.filtered_relative_hsps = rel_hsps

        ## 1. Frameshift?
        if rel_hsps.majority_frame:

            ## 1.1 Frameshift and missing 5'?
            if rel_hsps.missing_5prime:
                pass
            else:
                pass
        # 2. Missing 5'?
        elif rel_hsps.missing_5prime:
            pass

        # 3. Vanilla prediction
        else:
            pass
            
            
        
