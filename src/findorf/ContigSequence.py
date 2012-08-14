"""
ContigSequence.py contains the class declarations for ContigSequence,
HSP, ORF, and AnchorHSPs.

# ContigSequence

ContigSequence is first created by joining the FASTA file of contigs
with the separate XML BLASTX results against each relative's
databases. Then, in a prediction step, each ContigSequence's
`generic_predict_ORF` method is called, which dispatches the
appropriate information to separate function. These functions that
handle prediction are not class methods because:

1. Not all methods would be applicable, i.e. `predict_ORF_frameshift`
would not be useful if a ContigSequence's relatives all agree on the
frame (and then what, it returns None in this case)? This just seems
messy.

2. They're easier to unit test without having to make a bunch of
ContigSequence fixtures.

3. They're more modular, and their interfaces are generic (and
biologically rooted), so new rules or rule adjustments can be made
through simpler interfaces (like seq, frame, anchor_HSPs).

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
from collections import Counter, namedtuple, defaultdict
from string import Template
from operator import itemgetter, attrgetter
from math import floor

try:
    from Bio.Blast import NCBIXML
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError, e:
    sys.exit("Cannot import BioPython modules; please install it.")

from rules import get_closest_relative_anchor_HSP
from rules import get_anchor_HSPs
from rules import predict_ORF_frameshift, predict_ORF_missing_5prime
from rules import predict_ORF_vanilla
from rules import get_ORF_overlaps_5prime_HSP
from rules import annotate_ORF, START_CODONS, STOP_CODONS
from templates import anchor_hsps_repr, hsp_repr
from templates import orf_repr, contig_sequence_repr

GTF_FIELDS = ("seqname", "source", "feature", "start",
              "end", "score", "strand", "frame", "group")

class ORF():
    """
    ORF represents an ORF prediction.
    """

    def __init__(self, start, end, query_start, query_end,
                 frame, no_start=False, no_stop=False):
        self.start = start
        self.end = end
        self.query_start = query_start
        self.query_end = query_end
        self.frame = frame
        self.no_stop = no_stop
        self.no_start = no_start
        
    @property
    def length_bp(self):
        if None in (self.query_end, self.query_start):
            return None
        return abs(self.query_end - self.query_start)

    @property
    def length_aa(self):
        if None in (self.query_end, self.query_start):
            return None
        return int(floor(abs(self.query_end - self.query_start)/3))

    def __repr__(self):
        info = dict(start=self.start, end=self.end,
                    query_end=self.query_end,
                    query_start=self.query_start,
                    length_bp=self.length_bp,
                    length_aa=self.length_aa,
                    frame=self.frame,
                    missing_start=self.no_start,
                    missing_stop=self.no_stop)
        return Template(orf_repr).substitute(info)

    def get_orf(self, query_seq, include_stop=True):
        if self.frame < 0:
            query_seq = query_seq.reverse_complement()
        end = self.query_end
        if include_stop:
            end += 3

        # if we don't have a query start (is None), we make it the
        # absolute value of the frame minus 1, which is the 0-indexed
        # start.
        start = self.query_start
        if start is None:
            start = abs(self.frame) - 1
        return query_seq[start:end]

class HSP():
    """
    HSP represents a High-scoring Segment Pair from BLAST. BioPython
    has an HSP class, but for ORF prediction, a smaller subset of the
    attributes are needed.
    """
    def __init__(self, e, identities, length, percent_identity, title,
                 query_start, query_end, sbjct_start, sbjct_end, frame):
        self.e = e
        self.identities = identities
        self.length = length
        self.percent_identity = percent_identity
        self.title = title
        self.query_start = query_start
        self.query_end = query_end
        self.sbjct_start = sbjct_start
        self.sbjct_end = sbjct_end
        self.frame = frame

    def __repr__(self):
        info = dict(identities=self.identities,
                    percent_identity=round(self.percent_identity, 3),
                    length=self.length,
                    e=round(self.e, 12),
                    frame=self.frame,
                    query_start=self.query_start,
                    query_end=self.query_end,
                    sbjct_start=self.sbjct_start,
                    sbjct_end=self.sbjct_end)
        return Template(hsp_repr).substitute(info)


class AnchorHSPs():
    def __init__(self, most_5prime, most_3prime, strand):
        self.most_5prime = most_5prime
        self.most_3prime = most_3prime
        self.strand = strand

    def _indent_hsp(self, hsp, tab="  "):
        out = ""
        for line in repr(hsp).split("\n"):
            out += tab + line + "\n"
        return out

    def __repr__(self):
        info = dict(strand=self.strand,
                    most_3prime=self._indent_hsp(self.most_3prime),
                    most_5prime=self._indent_hsp(self.most_5prime))
        
        return Template(anchor_hsps_repr).substitute(info)

    def __iter__(self):
        for x in [self.most_5prime, self.most_3prime, self.strand]:
            yield x
         
        
class ContigSequence():
    """
    ContigSequence represents an assembled contig, that may be coding
    or non-coding. It contains all information about the its sequence
    and the blastx results to its relatives.
    """

    def __init__(self, query_id, sequence):
        """
        Initialize a ContigSequence with a contig ID and sequence. The
        contig ID must correspond to the same one used in the blastx
        results.
        """
        # core data attributes
        self.query_id = query_id
        self.seq = sequence
        self.len = len(sequence)

        # information added by blastx results
        self.all_relatives = defaultdict(list)

        # ORF and annotation attributes
        self.orf = None
        self.all_orfs = None
        self.annotation = dict()

        # These are run-specific (i.e. they are gathered from
        # arguments) parameters, and are hidden. The are stored
        # because __repr__() should output the same results as the
        # command line script run would, if python -i was called.
        self._e_value = None
        self._pi_range = None
        
    def __repr__(self):
        cr_anchor_hsp = get_closest_relative_anchor_HSP(self.get_anchor_HSPs())
        info = dict(id=self.query_id,
                    length=self.len,
                    num_relatives=self.num_relatives,
                    majority_frame=self.majority_frame,
                    majority_frameshift=self.majority_frameshift,
                    missing_stop=self.get_annotation('missing_stop'),
                    missing_start=self.get_annotation('missing_start'),
                    contains_stop=self.get_annotation('contains_stop'),
                    cr_frameshift=self.get_annotation('closest_relative_frameshift'),
                    orf_hsp_coverage=self.get_annotation('orf_hsp_coverage'),
                    missing_5prime=self.missing_5prime(self.get_anchor_HSPs()),
                    full_length_orf=self.get_annotation('full_length'),
                    anchor_hsps=repr(cr_anchor_hsp),
                    orf=repr(self.orf))
        if self.orf is not None:
            info['seq'] = self.orf.get_orf(self.seq)
        else:
            info['seq'] = None
        
        out = Template(contig_sequence_repr).substitute(info)
        return out
    
    def annotate_contig(self):
        """
        Annotate the contig, based on some attributes. We add these to
        a dictionary so cross-contig annotation counting is easier.
        """
        self.annotation['has_relatives'] = self.has_relatives
        self.annotation['num_relatives'] = self.num_relatives
        self.annotation['length'] = self.len
        self.annotation['majority_frameshift'] = self.majority_frameshift
        self.annotation['any_frameshift'] = self.any_frameshift
        
    def get_annotation(self, key):
        """
        Get annotation from key, or return None if it doesn't exist.
        """
        return self.annotation.get(key, None)

    def add_orf_prediction(self, orf):
        """
        Add a ORF.
        """
        self.orf = orf

    def add_annotation(self, annotation):
        """
        Annotation is just a dictinary of key-value pairs; this append
        those (or updates them if they exist)
        """
        self.annotation = dict(self.annotation.items() + annotation.items())
        

    def get_relatives(self):
        """
        Return relatives that pass thresholding filters.
        
        The `add_relative` method adds relatives' HSPs to a dictionary
        attribute, `all_relatives`. However, in most cases, we want to
        use a subset of these relatives that satisfy requirements
        based on phylogenetic requirements, i.e. requiring a relative
        HSP have a percent identity consistent with evolutionary
        distance. These constraints are run (via command line)
        specific, and are `_e_value` and `_pi_range`.

        If `_e_value` or `_pi_range` are None, they are not used for
        filtering `all_relatives`.
        """
        if self._e_value is None and self._pi_range is None:
            return self.all_relatives

        # little funcs for e-value filtering
        e_thresh = lambda x: x.e <= self._e_value

        filtered_relatives = defaultdict(list)
        for relative, hsps in self.all_relatives.items():

            filters = [(self._e_value, e_thresh)]
            # make a custom filter closure for this relative's range;
            # if a relative's range is None, we don't filter on it.
            
            if self._pi_range is not None:
                rng = self._pi_range[relative]
                in_range = (lambda x:
                            rng is None or rng[0] <= x.percent_identity <= rng[1])
                filters.append((self._pi_range, in_range))

            for h in hsps:
                if all([fun(h) for arg, fun in filters if arg is not None]):
                    filtered_relatives[relative].append(h)

        return filtered_relatives

    @property
    def num_relatives(self):
        """
        Return the number of relatives.
        """
        num_relatives = 0
        for relative, hsps in self.all_relatives.items():
            if len(hsps) > 0:
                num_relatives += 1
        return num_relatives

    @property
    def has_relatives(self):
        """
        Return True or False depending on whether the number of
        relatives is greater than 0.

        """
        return self.num_relatives > 0
        
        
    def gff_dict(self):
        """
        Return a dictionary of some key attribute's values,
        corresponding to a GFF file's columns.

        Note that GFFs are 1-indexed, so we add one to positions.
        """
        out = dict()
        out["seqname"] = self.query_id
        out["source"] = "findorf"
        out["feature"] = "predicted_orf"
        # we increment the start because GTF is 1-indexed, but not for
        # the end, since we want the ORF to (but not including) the
        # stop codon.
        out["start"] = self.orf.query_start + 1 if self.orf is not None else '.'
        out["end"] = self.orf.query_end if self.orf is not None else '.'
        out["score"] = "."

        if self.majority_frame is not None:
            out["strand"] = self.majority_frame/abs(self.majority_frame)
        else:
            out["strand"] = "."

        if self.majority_frame is not None:
            # GFF uses frames in [0, 2]
            out["frame"] = abs(self.majority_frame) - 1
        else:
             out["frame"] = "."
        out["group"] = "."
        return out

    def gtf_dict(self):
        """
        Return a dictionary corresponding to the columns of a GTF
        file.
        """
        self.annotate_contig()
        # a GTF's file's "group" column contains a merged set of
        # attributes, which in ContigSequence's case are those below
        group = "; ".join(["%s %s" % (k, v) for k, v in self.annotation.iteritems()])
        out = self.gff_dict()
        out["group"] = group
        return out

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

            self.all_relatives[relative].append(hsp)

    @property
    def frames(self):
        """
        Calculate and return the identity counts by relative and
        frame. This is used for both frameshift and any_frameshift.

        """

        frame_counts = defaultdict(Counter)
        for relative, hsps in self.all_relatives.iteritems():
            # count the number of identities per each relative's HSP
            for h in hsps:
                f = h.frame
                frame_counts[relative][f] += h.identities

        return frame_counts

    @property
    def majority_frame(self):
        """
        The `majority_frame` attribute indicates the majority, based
        on the number of *identities* that agree on a frame.
        
        This has the advantage that longer HSPs are weighted more
        heavily in the calculations. Furthermore, more distant
        relatives will likely be more divergent in terms of protein
        identity, so this provides a natural way of weighting by
        evolutionary distance.
        """
        if not self.has_relatives:
            return None

        frame = Counter()
        for relative, counts in self.frames.items():
            if len(counts) == 1:
                # no frameshifts in this relative
                frame[counts.keys()[0]] += sum(counts.values())

        if len(frame):
            majority_frame, count = frame.most_common(1)[0]
            return majority_frame

        return None
                
                
    @property
    def majority_frameshift(self):
        """
        Returns True of False if there's a frameshift in the majority
        of relatives, weighted by their identities.
        
        """

        if not self.has_relatives:
            return None

        frameshifts = Counter()
        for relative, counts in self.frames.items():
            frameshifts[len(counts.keys()) > 1] += sum(counts.values())
            
        return frameshifts[True] >= frameshifts[False]

    @property
    def any_frameshift(self):
        """
        Return if there are any relatives with frameshifts.
        
        """
        if not self.has_relatives:
            return None
    
        return any([len(c.keys()) > 1 for r, c in self.frames.iteritems()])

    @property
    def is_reversed(self):
        """
        Return True if the query is reversed.
        
        """
        if not self.has_relatives:
            return None

        return self.majority_frame < 0        

    def get_anchor_HSPs(self):
        """
        Get the 5'-most and 3'-most HSPs for each relative and put
        them in a tuple.

        We a generic `get_anchor_HSPs` method here, both because this
        method is useful outside of the `ContigSequence` class and
        because I wanted to unit test it outside of the
        `ContigSequence` class.

        Also note that `get_relatives()` will look at `_e_value` and
        `_pi_range`, run specific thresholding on which relatives to
        consider.
        """
        
        return get_anchor_HSPs(self.get_relatives(), self.is_reversed)

    def missing_5prime(self, anchor_hsps, qs_thresh=16, ss_thresh=40):
        """
        Return True if the anchor_hsps indicate a missing 5'-end of
        this contig.

        `qs_start` and `ss_thresh` are in amino acids.
        Each HSP has a query start and a subject start. A missing
        5'-end would look like this (in the case that the HSP spans
        the missing part):
        
                       query start
                      |   HSP
                 |------------------------------------------| contig
                      |||||||||||
              |.......|---------| subject
           subject
            start
        
        We infer missing 5'-end based on the query start position
        (compared to a threshold, `qs_thresh`) and the subject start
        position (`ss_thresh`). Starting late in the subject and early
        in the query probably means we're missing part of a protein.
        
        """

        if not len(anchor_hsps):
            return None
        
        missing_5prime = Counter()

        for relative, hsps in anchor_hsps.iteritems():
            most_5prime, most_3prime, strand = hsps
            # note that for reverse strand: query_start is really
            # query_end, but we compare it to the difference between
            # query_end and query_length.
            query_start = most_5prime.query_start
            sbjct_start = most_5prime.sbjct_start

            if strand > 0:
                m = query_start <= qs_thresh and sbjct_start >= ss_thresh
                missing_5prime[m] += 1
            else:
                # takte the query start and subtract it from length to
                # put everything on forward strand.
                qs = abs(query_start - self.len) + 1 # blast results are 1-indexed
                m = qs <= qs_thresh and sbjct_start >= ss_thresh
                missing_5prime[m] += 1
                  
        return missing_5prime[True] >= missing_5prime[False]

    def inconsistency(self):
        pass

    def generic_predict_ORF(self, e_value=None, pi_range=None):
        """
        The central dispatcher/logic behind ORF prediction.

        This method calls various functions in `rules`, which are
        side-effect free on this ContigSequence object.

        This is important because one could want to interactively run
        different predictions using different `e_value` and
        `pi_range`, without tangling up the `ContigSequence` object.
        """

        if self.has_relatives:
            # let's indicate these are a bit more private.
            self._e_value = e_value
            self._pi_range = pi_range

            anchor_hsps = self.get_anchor_HSPs()
            
            seq = self.seq

            ## The primary cases handled
            
            # we have relatives; we can predict the ORF
            if self.majority_frameshift:
                # we have a majority frameshift, so we let this
                # function use the anchor HSPs to deduce the 5'-most
                # frame
                missing_5prime = self.missing_5prime(anchor_hsps)
                orfs = predict_ORF_frameshift(seq, anchor_hsps, missing_5prime)

                # these could have separate missing 5'-ends, so
                # annotate those
                is_missing_5prime = self.missing_5prime(anchor_hsps) is True
                self.add_annotation({'missing_5prime':is_missing_5prime})
            elif self.missing_5prime(anchor_hsps):
                # no frameshift, we know the frame (at least in a majority of cases)
                self.add_annotation({'missing_5prime':True})
                frame = self.majority_frame
                orfs = predict_ORF_missing_5prime(seq, frame)
            else:
                self.add_annotation({'missing_5prime':False})
                frame = self.majority_frame
                orfs = predict_ORF_vanilla(seq, frame)

            ## Now, we try take each ORF list and find the 5'-most
            ## candidate that overlaps the 5'-most HSP
            if orfs is not None:
                self.all_orfs = orfs
                overlap_tuple = get_ORF_overlaps_5prime_HSP(orfs, anchor_hsps)
                if overlap_tuple is None:
                    self.add_annotation({'orf_hsp_coverage':False})
                    return

                relative, orf = overlap_tuple
                self.add_annotation({'relative_orf_hsp_overlap':relative})

                # we may not have a majority frameshift, but maybe the
                # closest relative does
                cr_ahsp = anchor_hsps[relative]
                cr_fs = cr_ahsp.most_5prime.frame != cr_ahsp.most_3prime.frame
                self.add_annotation({'closest_relative_frameshift':cr_fs})
                
                # with an ORF, we annotate it
                orf_annotation = annotate_ORF(anchor_hsps, orf)
                self.add_orf_prediction(orf)
                self.add_annotation(orf_annotation)
                self.add_annotation({'orf_hsp_coverage':True})

                    
