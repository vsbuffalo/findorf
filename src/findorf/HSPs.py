"""
HSPs.py

"""

from collections import defaultdict

class SeqRange():
    """
    A very simple range on a sequence.
    
    """

    def __init__(self, start, stop, strand):
        """
        Make a new SeqRange.
        
        """
        try:
            assert(self.stop >= self.start)
        except AssertionError:
            raise TypeError("start must be <= stop")

        self.start = start
        self.stop = stop
        self.strand = strand

    def __len__(self):
        """
        The length; stop should always be greater than or equal to
        start (equal in case of SNP), so we assert this first.
        """
        return self.stop - self.start

    def overlaps(self, other):
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

        return other.start <= self.stop and self.start <= other.stop
        
class HSP(SeqRange):
    """
    HSPs are like BioPython's SeqFeatures: a range on a sequence,
    stranded, with some metadata/annotation.
    
    """

    
    def __init__(self, e, identities, length, percent_identity, title,
                 query_start, query_end, sbjct_start, sbjct_end, frame):

        strand = -1 if frame < 0 else 1

        # initialize the seq range using the query
        super(SeqRange, self).__init__(start=query_start, stop=query_end,
                                       strand=strand)
        
        self.e = e
        self.identities = identities
        self.length = length
        self.percent_identity = percent_identity
        self.title = title
        self.sbjct_start = sbjct_start
        self.sbjct_end = sbjct_end
        self.frame = frame

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
        if other.__class__.__name__ != "SeqRange":
            raise TypeError("overlaps() require SeqRange object.")

        if other.strand != self.strand:
            raise TypeError("SeqRange objects must both have same strand")

        return self.most_5prime.overlaps(other)

class RelativeHSPs():
    """
    RelativeAnchorHSPs' primary attribute is a default dict of
    relatives' HSPs.

    """

    def __init__(self):
        self.relatives = defaultdict(list)
        self._e_value = None
        self._pi_range = None
        
    def add_relative_hsp(self, relative, hsp):
        self.relatives[relative].append(hsp)

    def get_relatives(self, e_value, pi_range):
        """
        The `add_relative` method adds relatives' HSPs to a dictionary
        attribute, `all_relatives`. However, in most cases, we want to
        use a subset of these relatives that satisfy requirements
        based on phylogenetic requirements, i.e. requiring a relative
        HSP have a percent identity consistent with evolutionary
        distance. These constraints are run (via command line)
        specific.
        """

        if self.e_value is None and self.pi_range is None:
            return self.relatives
        
        # little funcs for e-value filtering
        e_thresh = lambda x: x.e <= self._e_value

        filtered_relatives = defaultdict(list)
        for relative, hsps in self.all_relatives.items():

            filters = [(self.e_value, e_thresh)]
            # make a custom filter closure for this relative's range;
            # if a relative's range is None, we don't filter on it.
            
            if self.pi_range is not None:
                rng = self.pi_range[relative]
                in_range = (lambda x:
                            rng is None or rng[0] <= x.percent_identity <= rng[1])
                filters.append((self.pi_range, in_range))

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

    
