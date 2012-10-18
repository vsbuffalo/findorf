"""
contig.py contains the Contig class. This new version is a depature
from the previous version in that (1) it relies upon BioRanges and (2)
it stores uses objects from BioPython directly.

"""

from BioRanges.lightweight import Range, SeqRange, SeqRanges
from collections import Counter, namedtuple
from operator import itemgetter
from itertools import groupby

# named tuples used to improve readability 
AnchorHSPs = namedtuple('AnchorHSPs', ['most_5prime', 'most_3prime'])

DEFAULT_MIN_EXPECT = 10

def _HSP_to_dict(hsp):
    """
    Convert an HSP object from BioPython's HSP class to a dictionary
    for use with SeqRanges. This is a lightweight format, since we
    don't need everything from the full class. We don't store query
    start and end positions because this will be stored in the
    SeqRange object.

    We also do some sanity checking here.
    """
    # the BioPython parser doesn't give us a non-zero second
    # frame (which is for use with non-blastx parsers).
    assert(hsp.frame[1] is 0)
    
    # blastx has protein subjects, so this should always be the case
    assert(hsp.sbjct_start < hsp.sbjct_end)

    percent_identity = hsp.identities/float(hsp.align_length)

    hsp = dict(expect=hsp.expect,
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
        self.has_relative = False
        
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
            strand = "-" if hsp.frame[0] < 0 else "+"
            assert(qstart < qend)

            data = _HSP_to_dict(hsp)
            data.update({"relative":relative, "title":best_alignment.title})
            seqrng = SeqRange(Range(qstart, qend),
                              seqname=self.record.id,
                              strand=strand,
                              data=data)
            self.hsps.append(seqrng)

        self.has_relative = True
    
    def get_anchor_HSPs(self, min_expect=DEFAULT_MIN_EXPECT):
        """
        Get the 5'-most and 3'-most HSPs, handling the possibility the
        contig is in the reverse orientation.

        Note that this is O(n). We could potentially add methods to
        SeqRanges to make this faster, but this would depend on using
        a tree of some sort during construction.
        """
        if not self.has_relative:
            return None

        filtered_hsps = filter(lambda x: x['expect'] <= min_expect, self.hsps)

        # We get the outermost HSPs
        i = sorted(range(len(self.hsps)), key=lambda k: self.hsps.end[k], reverse=True)[0]
        j = sorted(range(len(self.hsps)), key=lambda k: self.hsps.start[k])[0]

        if self.majority_frame(min_expect) < 0:
            # negative strand; 5'-most HSP is that with the largest
            # query end
            return AnchorHSPs(self.hsps[i], self.hsps[j])
        else:
            # positive strand; 5-most HSP is that with the smallest
            # query start
            return AnchorHSPs(self.hsps[j], self.hsps[i])

    def count_frames(self, min_expect=DEFAULT_MIN_EXPECT):
        """
        Count the frames (by identities in that frame) of all HSPs,
        for use with majority_frame() nad majority_frameshift()
        methods.
        """        
        # filter SeqRange objects by whether they meet min e-value
        # requirements
        filtered_hsps = filter(lambda x: x['expect'] <= min_expect, self.hsps)
        
        frames = Counter()
        for hsp in filtered_hsps:
            frames[hsp['frame']] += hsp['identities']
        return frames
        
    def majority_frame(self, min_expect=DEFAULT_MIN_EXPECT):
        """
        Get the majority frame by looking at relatives' HSPs on the
        contig. The frame with the most identities backing it up is
        the majority frame.
        """
        if not self.has_relative:
            return None
        frames = self.count_frames(min_expect)
        if len(frames):
            majority_frame, _ = frames.most_common(1)[0]
            return majority_frame
        return None

    def any_frameshift(self, min_expect=DEFAULT_MIN_EXPECT):
        """
        Return whether there's any frameshift by looking at
        relatives' HSPs on the contig.
        """
        if not self.has_relative:
            return None
        frames = self.count_frames(min_expect)
        if len(frames):
            return len(set(frames)) > 1
        return None

    def majority_frameshift(self, min_expect=DEFAULT_MIN_EXPECT):
        """
        Return whether there's a majority frameshift by looking at
        relatives' HSPs on the contig. Majority frameshift is defined
        as whether there's a frameshift in most relatives.

        There's always the possibility that we're hitting a paralog
        with a frameshift in distant realtives via BLASTX. This is why
        our 'score' is based on the identities in each HSP. If not, a
        low-identity hit from relative with a frameshift would have
        been as equally weighted as a relative with high identity.
        """
        if not self.has_relative:
            return None

        filtered_hsps = filter(lambda x: x['expect'] <= min_expect, self.hsps)
        frames = [(h['relative'], h['frame'], h['identities']) for h in filtered_hsps]
        frames = sorted(frames, key=itemgetter(0))
        frameshifts = Counter()
        
        # here, we group by relative to see if the majority of
        # relatives have a frameshift. Each of the HSPs of each
        # relative are grouped and the number of identities is kept as
        # a tally of support
        for relative, hsps_info in groupby(frames, itemgetter(0)):
            # look to see whether this relative has HSPs in different
            # frames (that is, the set of frames has cardinality > 1)
            hsps_info = list(hsps_info)
            has_multiple_frames = len(set(map(itemgetter(1), hsps_info))) > 1

            num_identities = map(itemgetter(2), hsps_info)
            frameshifts[has_multiple_frames] += sum(num_identities)

        return frameshifts[True] >= frameshifts[False]

    def missing_5prime(self, qs_thresh=16, ss_thresh=40,
                       min_expect=DEFAULT_MIN_EXPECT):
        """
        Return True if the anchor HSPsx indicate a missing 5'-end of
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

        if not self.has_relative:
            return None

        most_5prime, most_3prime = self.get_anchor_HSPs()

        if most_5prime.strand == "-":
            missing = (most_5prime.start <= qs_thresh and
                       most_5prime['sbjct_start'] >= ss_thresh)
        else:
            qs = abs(most_5prime.start - most_5prime.width) + 1 # blast results are 1-indexed
            missing = qs <= qs_thresh and most_5prime['sbjct_start'] >= ss_thresh

        return missing
    
        
