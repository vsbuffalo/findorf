"""
Contig

# Blastx Output Notes

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
import predict

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
        self.orf = None
        self.orf_candidates = list()
        self.annotation = dict()

    def __repr__(self):
        info = dict(id=self.id, len=len(self),
                    num_relatives=len(self.relative_hsps))
        msg = Template("""Contig instance for ID: $id
length: $len
number of relatives: $num_relatives
""").substitute(info)

        return msg


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

    def add_annotation(self, key, value):
        """
        Add annotation, keeping track of consistency.
        """
        self.annotation[key] = value
        
    def get_annotation(self, key=None):
        """
        Get annotation from the annotation dict, but return None if
        there is something is not defined.
        """
        if key is not None:
            return self.annotation.get(key, None)
        
        current_anno = self.annotation.items()
        attribute_anno = [('has_relatives', self.has_relatives),
                          ('num_relatives', len(self.relative_hsps)),
                           ('has_orf', self.orf is not None),
                           ('num_orf_candidates', len(self.orf_candidates))]

        return dict(current_anno + attribute_anno)
        

    @property
    def has_relatives(self):
        return len(self.relative_hsps) > 0

    def predict_orf(self, e_value=None, pi_range=None, qs_thresh=16, ss_thresh=40):
        """
        Predict an ORF, given some parameters.

        1. Frameshift?
          1.1. Missing 5' also?

        2. 5' missing?

        3. Vanilla
        """

        if not self.has_relatives:
            return None

        # Filter by the parameters by getting only the relative HSPs
        # that satisfy our e-value and percent identity requirements;
        # we use these filtered cases from now. This is where
        # self.relative_hsps and rel_hsp can differ.
        rel_hsps = self.relative_hsps.get_relatives(e_value, pi_range)
        self.filtered_relative_hsps = rel_hsps

        # it's possible that self.relative_hsps has relatives, but the
        # filtered rel_hsps does not. If this is the case, we return
        # None

        if not rel_hsps.has_relatives:
            return None

        # get the frame and closest relative anchor HSPs and whether
        # the HSPs lie on different strands (really degenerate case).
        frame = rel_hsps.majority_frame
        inconsistent_strand = rel_hsps.inconsistent_strand
        self.add_annotation('inconsistent_strand', inconsistent_strand)
        
        if inconsistent_strand:
            return None
        
        closest_relative, cr_ahsps = rel_hsps.closest_relative_anchor_hsps()
        self.cr_ahsps = cr_ahsps # for interactive debugging
        self.add_annotation('closest_relative', closest_relative)
        
        # get whether the relative HSPs indicate a missing 5' or frameshift 
        missing_5prime = rel_hsps.missing_5prime(len(self), qs_thresh, ss_thresh)
        self.add_annotation('missing_5prime', missing_5prime)
        majority_frameshift = rel_hsps.majority_frameshift
        self.add_annotation('majority_frameshift', majority_frameshift)
        
        ## 1. Frameshift?
        if majority_frameshift:
            # note that we're handling 1.1 (passing in missing 5' as
            # whether we're in reading frame).
            orf_candidates = predict.orf_with_frameshift(self.seq, cr_ahsps,
                                                         missing_5prime)          

        # 2. Missing 5'?
        elif missing_5prime:
            orf_candidates = predict.orf_with_missing_5prime(self.seq, frame)

        # 3. Vanilla prediction
        else:
            orf_candidates = predict.orf_vanilla(self.seq, frame)

        self.orf_candidates = orf_candidates

        # With ORF candidates, we chose the one with an overlap with
        # the 5'-most closest relative's anchor HSP. We use some
        # partial currying to allow this to be done in a list
        # comphrension. 
        f = lambda x: cr_ahsps.orf_overlaps_5prime(x, len(self))
        overlap_candidates = [orf for orf in orf_candidates if f(orf)]
        self.orfs_overlap = overlap_candidates

        if len(overlap_candidates) == 0:
            return None
        
        self.orf = overlap_candidates[0]

        return self.orf
