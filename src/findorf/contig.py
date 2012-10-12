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
    from Bio.SeqRecord import SeqRecord
except ImportError, e:
    sys.exit("Cannot import BioPython modules; please install it.")

from RangedFeatures import HSP, AnchorHSPs, RelativeHSPs, indent
import predict
from templates import contig_str


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
        self.cr_ahsps = None
        self.filtered_relative_hsps = None
        self.orf = None
        self.missing_start = None
        self.missing_stop = None
        self.orf_candidates = list()
        self.annotation = dict()

    def __repr__(self):
        """
        A string summary of the Contig object.
        """
        out = '-' * 80 + "\n"
        info = dict(id=self.id, length=len(self),
                    num_relatives=len(self.relative_hsps),
                    majority_frame=self.relative_hsps.majority_frame,
                    majority_frameshift=self.relative_hsps.majority_frameshift,
                    inconsistent_strand=self.relative_hsps.inconsistent_strand,
                    relative=self.get_annotation('closest_relative'),
                    missing_5prime=self.relative_hsps.missing_5prime(len(self)),
                    overlap=self.get_annotation('hsp_orf_overlap'),
                    orf=self.orf)
        
        if hasattr(self, 'cr_ahsps'):
            info.update(dict(anchor_hsps=self.cr_ahsps))
        else:
            info.update(dict(anchor_hsps=None))

        if self.orf is not None:
            tmp = dict(missing_start=self.orf.no_start,
                       missing_stop=self.orf.no_stop,
                       seq=self.orf.get_sequence(self))
        else:
            tmp = dict(missing_start=None,
                       missing_stop=None,
                       seq=None)
            
        info.update(tmp)
        out += Template(contig_str).substitute(info)
        return out

    def __str__(self):
        self.__repr__()

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
        there is something is not defined. These are updated first,
        before retrieving.
        """
        current_anno = self.annotation.items()
        attribute_anno = [('has_relatives', self.has_relatives),
                          ('num_relatives', len(self.relative_hsps)),
                          ('has_orf', self.orf is not None),
                          ('contig_length', len(self)),
                          ('num_orf_candidates', len(self.orf_candidates))]

        self.annotation =  dict(current_anno + attribute_anno)
        
        if key is not None:
            return self.annotation.get(key, None)

        return self.annotation
        

    @property
    def has_relatives(self):
        return len(self.relative_hsps) > 0

    def gff_dict(self):
        """
        Return a dictionary of some key attribute's values,
        corresponding to a GFF file's columns.
        
        Note that GFFs are 1-indexed, so we add one to positions.
        """
        out = dict()
        out["seqname"] = self.id
        out["source"] = "findorf"
        out["feature"] = "predicted_orf"
        # we increment the start because GTF is 1-indexed, but not for
        # the end, since we want the ORF to (but not including) the
        # stop codon.
        out["start"] = self.orf.abs_start() + 1 if self.orf is not None else '.'
        out["end"] = self.orf.abs_end(len(self)) + 1 if self.orf is not None else '.'
        out["score"] = "."

        mf = self.relative_hsps.majority_frame
        if mf is not None:
            out["strand"] = mf/abs(mf)
        else:
            out["strand"] = "."
            
        if mf is not None:
            # GFF uses frames in [0, 2]
            out["frame"] = abs(mf) - 1
        else:
            out["frame"] = "."
            out["group"] = "."
        return out
        
    def gtf_dict(self):
        """
        Return a dictionary corresponding to the columns of a GTF
        file.
        """
        anno = self.get_annotation()
        # a GTF's file's "group" column contains a merged set of
        # attributes, which in ContigSequence's case are those below
        group = ";".join(["%s %s" % (k, v) for k, v in anno.items()])
        out = self.gff_dict()
        out["group"] = group
        return out

    @property
    def protein(self):
        """
        Return a protein sequence.
        """
        if self.orf is not None:
            seq = self.orf.get_sequence(self)
            frame = self.relative_hsps.majority_frame
            return SeqRecord(seq=seq.translate(), id=self.id)
        return None

    @property
    def orf_seq(self):
        """
        Return the nucleotide sequence record
        """
        if self.orf is not None:
            seq = self.orf.get_sequence(self)
            return SeqRecord(seq=seq, id=self.id)
        return None

    def five_prime_utr(self):
        """
        Return the 5'-UTR sequences of contigs with ORFS (full ORFs
        and partials).
        """
        if self.orf is None:
            return None

        return self.orf.get_5prime_utr(self)
        
        
    def three_prime_utr(self):
        """
        Return the 3'-UTR sequences of contigs with ORFs (full ORFs
        and partials).
        """
        if self.orf is None:
            return None

        return self.orf.get_3prime_utr(self)
        

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

        if len(overlap_candidates) == 0:
            self.add_annotation('hsp_orf_overlap', False)
            return None
        else:
            self.add_annotation('hsp_orf_overlap', True)

        self.orf = overlap_candidates[0]
        self.orfs_overlap = overlap_candidates

        # Annotate missing start/stop for ORF
        self.add_annotation('missing_start', self.orf.no_start)
        self.add_annotation('missing_stop', self.orf.no_stop)
        self.add_annotation('full_length', all((not self.orf.no_start, not self.orf.no_start)))
        self.add_annotation('internal_stop', self.internal_stop)
        return self.orf

    @property
    def internal_stop(self):
        """
        If there's an ORF predicted, check that there are no (1) no
        other ORF candidats with overlaps (that's an internal stop)
        and (2) that the 3' ends are not after the stop position.
        """
        if self.orf is None:
            # we need an ORF for this - the ORF before the stop.
            return False

        if len(self.orfs_overlap) > 1:
            return True

        # getting the frame this way is necessary; there could be a
        # frameshift (majority frame is None then) and a stop codon
        frame = self.cr_ahsps.most_5prime.frame
        l = len(self)

        ahsps = self.relative_hsps.get_anchor_hsps()
        for relative, ahsps in ahsps:
            if frame < 0:
                if self.orf.end < self.cr_ahsps.most_3prime.put_on_forward_strand(l).start:
                    return True
            else:
                if self.orf.end < self.cr_ahsps.most_3prime.start:
                    return True

        return False
                
