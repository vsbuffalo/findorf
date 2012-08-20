"""
RangesFeatures.py contains *lightweight* classes for storing the bare
minimum amount of information and functional of ranges on sequences to
predict an ORF. Lightweight range structures are used (i.e. no
metainformation, no interval tree backend) because we don't have to do
too many overlap calculations.

"""
import pdb
from collections import defaultdict, Counter
from copy import deepcopy
from operator import attrgetter
from string import Template
        
def indent(string, level=2, char=' '):
    s = string.split("\n")
    return char*level + '\n'.join([char*level + line for line in s])

def nested_attrgetter(attr_1, attr_2):
        """
        Like operator.attrgetter, but for nested objects like
        AnchorHSPs.

        This may look strange, but recall attrgetter returns a
        function, so this is essentially currying.

        Note that we use x[1], as this is fed in tuples of dict's
        (key, values).
        """
        
        return lambda x: attrgetter(attr_1)(attrgetter(attr_2)(x[1]))

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

        if None not in (self.start, self.end, other.start, other.end):
            return other.start <= self.end and self.start <= other.end
        else:
            if not allow_unbounded:
                return None

        # unbounded case; for now, let's let only `other` contain
        # unbounded.
        if None in (self.start, self.end):
            raise ValueError("overlaps() only allows None "
                             "to be only in the 'other' SeqRange;"
                             "the SeqRange must have specified start/end positiosn.")

        if None in (other.start, other.end) and other.seqlength is None:
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

    def __repr__(self):
        info = (self.start, self.end,
                self.sbjct_start, self.sbjct_end,
                self.frame, round(self.percent_identity, 4),
                round(self.e, 4))
        return "HSP(qs:%s, qe:%s, ss:%s, se:%s, f:%s, pe:%s, e:%s)" % info

    def __str__(self):
        out = """
query start/end: $start/$end
subject start/end: $sbjct_start/$sbjct_end
frame: $frame
e-value: $e
identities: $identities
length: $length
percent identity: $percent_identity
"""
        return Template(out).substitute(self.__dict__)        
        

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

        is_reversed = hsps[0].frame < 0
        self.strand = -1 if is_reversed else 1
        
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

    def __repr__(self):
        return "AnchorHSPs(5': %s; 3': %s)" % (repr(self.most_5prime),
                                               repr(self.most_3prime))
    def __str__(self):
        hsps = (indent(str(self.most_5prime)), indent(str(self.most_3prime)))
        return "most 5':\n%s\nmost 3':\n%s\n" % hsps

    def __iter__(self):
        for x in [self.most_5prime, self.most_3prime, self.strand]:
            yield x

    def orf_overlaps_5prime(self, orf, query_length):
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
        if self.most_5prime.frame < 0:
            # pdb.set_trace()
            m5p_ahsp = self.most_5prime.put_on_forward_strand(query_length)
        else:
            m5p_ahsp = self.most_5prime
            
        return m5p_ahsp.overlaps(orf.range)

class RelativeHSPs():
    """
    RelativeAnchorHSPs' primary attribute is a default dict of
    relatives' HSPs.

    """

    def __init__(self):
        self.relatives = defaultdict(list)
        self.anchor_hsps = dict()
        self._e_value = None
        self._pi_range = None

    def _pretty_anchor_hsps(self, ahsps=None):
        """
        Return a formatted string of relatives and their anchor
        HSPs. Expects a list of (relative, AnchorHSPs) tuple, or if
        None, will print all anchor HSPs.
        """
        if ahsps is None:
            self.get_anchor_hsps()
            ahsps = self.anchor_hsps.items()
        out = ""
        for relative, hsps in ahsps:
            out += "\nrelative: %s" % relative
            out += indent(str(hsps))
        return out

    def __repr__(self):
        self.get_anchor_hsps()
        closest_rel, cr_ahsps = self.closest_relative_anchor_hsps()
        info = (len(self), self._pretty_anchor_hsps([(closest_rel, cr_ahsps)]))
        return "RelativeHSPs\n%s relatives\nAnchorHSPs:\n%s" % info

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
            # no filtering required
            return self
        
        # little funcs for e-value filtering
        e_thresh = lambda x: x.e <= e_value

        filtered_relatives = defaultdict(list)
        for relative, hsps in self.relatives.items():

            filters = [(e_value, e_thresh)]
            # make a custom filter closure for this relative's range;
            # if a relative's range is None, we don't filter on it.
            
            if pi_range is not None:
                rng = pi_range[relative]
                in_range = (lambda x:
                            rng is None or rng[0] <= 100*x.percent_identity <= rng[1])
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

    @property
    def has_relatives(self):
        return len(self.relatives) > 0

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
    def inconsistent_strand(self):
        """
        In some cases, we may have a majority frameshift, but also
        because the HSPs are on different strands. This is a very
        degenerate case, and should be annotated as such.
        """

        if not self.has_relatives:
            return None

        inconsistent_strand = Counter()
        for relative, counts in self._frames.items():
            # look for differing frame, differing strand
            strands = [f/abs(f) for f in counts.keys()]
            inconsistent_strand[len(set(strands)) > 1] += 1

        return inconsistent_strand[True] >= inconsistent_strand[False]

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

    def get_anchor_hsps(self):
        """
        Get (and add to the anchor_hsps attribute) the 5' and 3' HSPs.

        We ignore cases in which there is a relative with HSPs on
        different strands (but allow for differing frame)
        """
        # not majority inconsistent strand, individual relatives could
        # still have HSPs on differing strands.
        assert(not self.inconsistent_strand)

        for relative, hsps in self.relatives.items():
            strands = [h.frame/abs(h.frame) for h in hsps]
            if len(set(strands)) > 1:
                continue
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

        anchor_hsps = self.get_anchor_hsps()
        tmp = sorted(anchor_hsps.items(),
                     key=nested_attrgetter(key, which),
                     reverse=reverse)
        
        return tmp[0]

    
    def missing_5prime(self, query_length, qs_thresh=16, ss_thresh=40):
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
            most_5prime, most_3prime, strand = anchor_hsps
            
            # note that for reverse strand: the 5'-most is really
            # query_end, but we compare it to the difference between
            # query_end and query_length.
            query_start = most_5prime.start
            sbjct_start = most_5prime.sbjct_start

            if strand > 0:
                m = query_start <= qs_thresh and sbjct_start >= ss_thresh
                missing_5prime[m] += 1
            else:
                # take the query start and subtract it from length to
                # put everything on forward strand.
                qs = abs(query_start - query_length) + 1 # blast results are 1-indexed
                m = qs <= qs_thresh and sbjct_start >= ss_thresh
                missing_5prime[m] += 1
                  
        return missing_5prime[True] >= missing_5prime[False]


class ORF():
    """
    An ORF or ORF candidate.
    """
    
    def __init__(self, query_start, query_end, query_length, frame,
                 no_start=None, no_stop=None):
        """
        Note that query_start and query_send are on the forward
        strand.
        """
        self.start = start = query_start
        self.end = end = query_end
        self.frame = frame
        self.no_start = no_start
        self.no_stop = no_stop

        if no_start:
            start = 0
        if no_stop:
            end = query_length
        
        self.range = SeqRange(start, end, strand=1)

    def abs_start(self):
        """
        abs_start returns the start position without None.
        """
        if self.start is None:
            if not self.no_start:
                raise ValueError("inconsistent ORF: no_start and start don't agree")
            return 0
        else:
            return self.start

    def abs_end(self, query_length):
        if self.end is None:
            if not self.no_stop:
                raise ValueError("inconsistent ORF: no_stop and end don't agree")
            return query_length
        else:
            return self.end

        
    def __repr__(self):
        return "ORF(%s-%s, %s)" % (self.start, self.end, self.frame)

    def get_sequence(self, contig, include_stop=True):
        if contig.__class__.__name__ != "Contig":
            raise TypeError("'contig' must be a Contig object")
        end = self.end
        if include_stop:
            end = end + 3
        if self.frame < 0:
            return contig.seq.reverse_complement()[self.start:end]
        return contig.seq[self.start:end]

    def __len__(self):
        return len(self.range)

