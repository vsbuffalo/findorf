"""
RangesFeatures.py

"""

from collections import defaultdict, Counter
from copy import deepcopy

def nested_attrgetter(attr_1, attr_2):
        """
        Like operator.attrgetter, but for nested objects like
        AnchorHSPs.

        This may look strange, but recall attrgetter returns a
        function, so this is essentially currying.
        """
        
        return lambda x: attrgetter(attr_1)(attrgetter(attr_2)(x))

class SeqRange(object):
    """
    A very simple range on a sequence.    
    """

    def __init__(self, start, end, strand, seqname=None, seqlength=None):
        """
        Make a new SeqRange.
        
        """
        try:
            if None not in (end, start):
                assert(end >= start)
        except AssertionError:
            raise TypeError("start must be <= end")

        self.start = start
        self.end = end
        self.strand = strand
        self.seqname = seqname
        self.seqlength = seqlength

    def __len__(self):
        """
        The length; end should always be greater than or equal to
        start (equal in case of SNP), so we assert this first.

        Forward:
        
        No end:
        start = 6
        ATTATAGGAGGAGTCTAGAATAG
              |----------------

        No start:
        end = 18
        ATTATAGGAGGAGTCTAGAATAG
        -----------------|
        
        Reverse:
        No end:
        start = 6
        ATTATAGGAGGAGTCTAGAATAG
        -----------------|      

        No start:
        end = 14
        ATTATAGGAGGAGTCTAGAATAG
                |--------------

        """
        if None not in (self.start, self.end):
            return self.end - self.start

        if self.start is None:
            # missing start, the length is 0 to the end.
            return self.end

        if self.end is None:
            # missing end, the length is the start to the end of the
            # sequence.
            if seqlength is None:
                return None
            return seqlength - self.start

    def overlaps(self, other, allow_unbounded=True):
        """
        Given another SeqRange, return true of false if it overlaps.

          |---------------| a
                  |------------| b
    
         |---------------| a
                         |----| b

        or              
                     |------------| a
               |----------| b

        """

        if other.__class__.__name__ != "SeqRange":
            raise TypeError("overlaps() require SeqRange object.")

        if other.strand != self.strand:
            raise TypeError("SeqRange objects must both have same strand")

        if not allow_unbounded:
            if None in (self.start, self.end, other.start, other.end):
                return None
            return other.start <= self.end and self.start <= other.end

        # unbounded case; for now, let's let only `other` contain
        # unbounded.
        if None in (self.start, self.end):
            raise ValueError("overlaps() only allows None "
                             "to be in the 'other' SeqRange.")
        if other.seqlength is None:
            raise ValueError("overlaps() requires 'other' to have seqlength "
                             "if it's missing a start or end position.")
        if other.start is None:
            start = 0
        if other.end is None:
            end = other.seqlength

        return start <= self.end and self.start <= end
        
        
class HSP(SeqRange):
    """
    HSPs are like BioPython's SeqFeatures: a range on a sequence,
    stranded, with some metadata/annotation.
    
    """
    
    def __init__(self, e, identities, length, percent_identity, title,
                 query_start, query_end, sbjct_start, sbjct_end, frame):

        strand = -1 if frame < 0 else 1

        # initialize the seq range using the query
        super(HSP, self).__init__(start=query_start, end=query_end,
                                       strand=strand)

        self.strand = strand
        self.e = e
        self.identities = identities
        self.length = length
        self.percent_identity = percent_identity
        self.title = title
        self.sbjct_start = sbjct_start
        self.sbjct_end = sbjct_end
        self.frame = frame

    def put_on_forward_strand(self, query_length):
        """
        Return a new copy of this HSP as if it were on the forward
        strand.
        
        No side effects; do not actual mutate this object to put it on
        forward strand.
        """
        if self.frame > 0:
            return self
        hsp = deepcopy(self)
        hsp.strand = 1
        hsp.frame = abs(self.frame)

        hsp.start = query_length - self.end + 1
        hsp.end = query_length - self.start + 1
        assert(hsp.end >= hsp.start)
        
        return hsp

class AnchorHSPs():
    """
    AnchorHSPs represent two HSPs, the 5'-most and 3'-most on a query,
    which could be the same if there's just HSP.
    """

    def __init__(self, hsps):
        """
        Create a new AnchorHSPs set, given a list of HSPs. It will
        find the 5'-most and 3'-most.

        """
        consistent_strand = len(set([h.frame < 0 for h in hsps])) == 1
        if not consistent_strand:
            raise TypeError("HSPs must have same strand")

        is_reversed = h[0].frame < 0

        hsp_1 = sorted(hsps, key=attrgetter('end'), reverse=True)[0]
        hsp_2 = sorted(hsps, key=attrgetter('start'))[0]

        if is_reversed:
            # if reversed, the HSP closest to the protein N-terminus
            # is the one with the latest end position
            self.most_5prime = hsp_1
            self.most_3prime = hsp_2
        else:
            # on the forward strand, the opposite is true.
            self.most_5prime = hsp_2
            self.most_3prime = hsp_1


    def __iter__(self):
        for x in [self.most_5prime, self.most_3prime, self.strand]:
            yield x

    def overlaps_5prime(self, other):
        """
        Check if another SeqRange overlaps the 5'-end anchor HSP.
        """
        
        return self.most_5prime.overlaps(other)

class RelativeHSPs():
    """
    RelativeAnchorHSPs' primary attribute is a default dict of
    relatives' HSPs.

    """

    def __init__(self):
        self.relatives = defaultdict(list)
        self.anchor_hsps = None
        self._e_value = None
        self._pi_range = None

    def add_relative_hsp(self, relative, hsp):
        self.relatives[relative].append(hsp)

    def get_relatives(self, e_value, pi_range):
        """
        The `add_relative` method adds relatives' HSPs to a dictionary
        attribute, `relatives`. However, in most cases, we want to
        use a subset of these relatives that satisfy requirements
        based on phylogenetic requirements, i.e. requiring a relative
        HSP have a percent identity consistent with evolutionary
        distance. These constraints are run (via command line)
        specific.
        """

        if e_value is None and pi_range is None:
            return self.relatives
        
        # little funcs for e-value filtering
        e_thresh = lambda x: x.e <= self._e_value

        filtered_relatives = defaultdict(list)
        for relative, hsps in self.relatives.items():

            filters = [(e_value, e_thresh)]
            # make a custom filter closure for this relative's range;
            # if a relative's range is None, we don't filter on it.
            
            if pi_range is not None:
                rng = pi_range[relative]
                in_range = (lambda x:
                            rng is None or rng[0] <= x.percent_identity <= rng[1])
                filters.append((pi_range, in_range))

            for h in hsps:
                if all([fun(h) for arg, fun in filters if arg is not None]):
                    filtered_relatives[relative].append(h)

        # returned filtered relatives
        new = RelativeHSPs()
        new._e_value = e_value
        new._pi_range = pi_range
        new.relatives = filtered_relatives
        return new

    def __len__(self):
        """
        The number of relatives.
        """
        return len(self.relatives)

    def has_relatives(self):
        return len(self.relative_hsps) > 0

    @property
    def _frames(self):
        """
        Calculate and return the identity counts by relative and
        frame. This is used for both frameshift and any_frameshift.

        """

        frame_counts = defaultdict(Counter)
        for relative, hsps in self.relatives.iteritems():
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
        for relative, counts in self._frames.items():
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
        for relative, counts in self._frames.items():
            frameshifts[len(counts.keys()) > 1] += sum(counts.values())
            
        return frameshifts[True] >= frameshifts[False]

    @property
    def any_frameshift(self):
        """
        Return if there are any relatives with frameshifts.
        
        """
        if not self.has_relatives:
            return None
    
        return any([len(c.keys()) > 1 for r, c in self._frames.iteritems()])

    @property
    def anchor_hsps(self):
        """
        Get (and add to the anchor_hsps attribute) the 5' and 3' HSPs.
        """
        for relative, hsps in self.relatives.items(): 
            self.anchor_hsps[relative] = AnchorHSPs(hsps)

        return self.anchor_hsps

    def closest_relative_anchor_hsps(self, key='e', which='most_5prime'):
        """
        Get (and add to the anchor_hsps attribute) the 5' and 3' HSPs
        of the closest_relative.
        """
        if not self.has_relatives:
            return None

        if key not in ('e', 'identity'):
            raise TypeError("key must be 'e' or 'identity'")

        reverse = False
        if key == 'identity':
            # more identities = closer; contrast more e-value =
            # fruther
            reverse = True

        tmp = sorted(self.anchor_hsps.iteritems(),
                     key=nested_attrgetter(key, which),
                     reverse=reverse)
        
        return tmp[0]

    
    def missing_5prime(self, qs_thresh=16, ss_thresh=40):
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

        if not self.has_relatives:
            return None
        
        missing_5prime = Counter()

        for relative, anchor_hsps in self.anchor_hsps.iteritems():
            most_5prime, most_3prime = anchor_hsps
            strand = anchor_hsps.strand
            
            # note that for reverse strand: the 5'-most is really
            # query_end, but we compare it to the difference between
            # query_end and query_length.
            query_start = most_5prime.query_start
            sbjct_start = most_5prime.sbjct_start

            if strand > 0:
                m = query_start <= qs_thresh and sbjct_start >= ss_thresh
                missing_5prime[m] += 1
            else:
                # take the query start and subtract it from length to
                # put everything on forward strand.
                qs = abs(query_start - self.len) + 1 # blast results are 1-indexed
                m = qs <= qs_thresh and sbjct_start >= ss_thresh
                missing_5prime[m] += 1
                  
        return missing_5prime[True] >= missing_5prime[False]


    def orf_overlaps_closest_relative_5prime_hsp(self, orf, query_length):
        """
        Does SeqRange `orf` (ALWAYS ON FORWARD STRAND), overlap the 5'
        anchor hsp? The 5'-anchor HSP is strand-aware: see AnchorHSPs.

        The forward strand is handled only because ORFs only make
        sense on the forward strand since transcription and
        translation are 5' to 3'.

        Currently, the 5' anchor HSP is based on the contig reference,
        which if the majority frame > 0, we can easily use existing
        overlap function with the `orf` coordinates.

        However, if the `orf` is on the reverse strand (majority frame
        < 0), we need to translate the 5' anchor HSP to the
        corresponding 5'-anchor HSP on the forwards strand to compute
        overlaps. Recall, like Bioconductor's GRanges, to calculate
        overlaps, we must be on the same strand.

        Also, we allow one-sided overlaps.
        """

        m5p_ahsp = self.get_closest_relative_anchor_HSP()

        if m5p_ahsp.frame < 0:
            m5p_ahsp = m5p_ahsp.put_on_forward_strand(query_length)

        return m5p_ahsp.overlaps(orf.range)


class ORF():
    """
    An ORF or ORF candidate.
    """
    
    def __init__(self, query_start, query_end, frame,
                 no_start=None, no_stop=None):
        """
        Note that query_start and query_send are on the forward
        strand.
        """
        self.start = query_start
        self.end = query_end
        self.frame = frame
        self.no_start = no_start
        self.no_stop = no_stop

        self.range = SeqRange(start, end, strand=1)

    def get_sequence(self, contig):
        if contig.__class__.__name__ != "Contig":
            raise TypeError("'contig' must be a Contig object")

        if frame < 0:
            return contig.seq.reverse_complement()[self.start:self.end]
        return contig.seq[self.start:self.end]

    def __len__(self):
        return len(self.range)

    
