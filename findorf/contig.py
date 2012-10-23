"""
contig.py contains the Contig class. This new version is a depature
from the previous version in that (1) it relies upon BioRanges and (2)
it stores uses objects from BioPython directly.

Strange cases:
 - k61_contig_13995

"""

import pdb
from collections import Counter, namedtuple
from operator import itemgetter
from itertools import groupby

from BioRanges.lightweight import Range, SeqRange, SeqRanges
from orfprediction import get_all_orfs

# named tuples used to improve readability 
AnchorHSPs = namedtuple('AnchorHSPs', ['relative', 'most_5prime', 'most_3prime'])

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
        self.pfam_domains = SeqRanges()
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

            # Adjust BLAST's 1-based indexing to our 0-based indexing.
            qstart = hsp.query_start - 1
            qend = hsp.query_end - 1
            strand = "-" if hsp.frame[0] < 0 else "+"
            assert(qstart <= qend)

            data = _HSP_to_dict(hsp)
            data.update({"relative":relative, "title":best_alignment.title})
            seqrng = SeqRange(Range(qstart, qend),
                              seqname=self.record.id,
                              strand=strand,
                              seqlength=len(self.record.seq),
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

        if self.get_strand(min_expect) == "-":
            # negative strand; 5'-most HSP is that with the largest
            # query end
            return AnchorHSPs(self.hsps[i]['relative'], self.hsps[i], self.hsps[j])
        else:
            # positive strand; 5-most HSP is that with the smallest
            # query start
            return AnchorHSPs(self.hsps[j]['relative'], self.hsps[j], self.hsps[i])

    def get_strand(self, min_expect=DEFAULT_MIN_EXPECT):
        """
        Get a strand (+, -), a step we can do before guess frame.

        We need strand to infer 5'-anchor HSPs, which we need if we
        have a frameshift, so this must be found before frame.
        """
        if not self.has_relative or self.inconsistent_strand(min_expect):
            return None

        filtered_hsps = filter(lambda x: x['expect'] <= min_expect, self.hsps)

        strands = [h.strand for h in filtered_hsps]
        # assert cardinality of strand set is 1
        assert(len(set(strands)) == 1)
        return strands[0]

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

    def inconsistent_strand(self, min_expect=DEFAULT_MIN_EXPECT):
        """
        In some cases, we may have a majority frameshift, but also
        because the HSPs are on different strands. This is a very
        degenerate case, and should be annotated as such.
        """
        filtered_hsps = filter(lambda x: x['expect'] <= min_expect, self.hsps)
        frames = [(h['relative'], h['frame'], h['identities']) for h in filtered_hsps]
        frames = sorted(frames, key=itemgetter(0))
        inconsistent_strand = Counter()
        
        for relative, hsps_info in groupby(frames, itemgetter(0)):
            hsps_info = list(hsps_info)
            frames = map(itemgetter(1), hsps_info)
            strands = [f/abs(f) for f in frames]
            
            inconsistent_strand[len(set(strands)) > 1] += 1

        return inconsistent_strand[True] >= inconsistent_strand        

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

        If missing_5prime is True, then we also want to consider the
        case that there is not start codon 5' of the 5'-most HSP --
        this is a case where the start is 100% missing.
        
        """

        if not self.has_relative:
            return None

        most_5prime_relative, most_5prime, most_3prime = self.get_anchor_HSPs(min_expect)
        
        if most_5prime.strand == "-":
            missing = (most_5prime.start <= qs_thresh and
                       most_5prime['sbjct_start'] >= ss_thresh)
        else:
            qs = abs(most_5prime.start - most_5prime.width) + 1 # blast results are 1-indexed
            missing = qs <= qs_thresh and most_5prime['sbjct_start'] >= ss_thresh

        return missing

    def add_pfam(self, domain_hit_seqrange):
        """
        Add PFAM domain hit (from HMMER). Note that all of the
        coordinate conversion is done via add_pfam_domain_hits()
        function in the hmmer module.
        """
        self.pfam_domains.append(domain_hit_seqrange)

    def more_5prime_pfam_domain(self, most_5prime_hsp, frame, min_expect=DEFAULT_MIN_EXPECT):
        """
        Return PFAM domain more 5' prime of supplied SeqRange object
        (which should be the 5' anchor HSP), or None of if there is none.

        Note that all PFAM domains are on the positive strand, since
        PFAM domains found via HMMSCAN were in protein space.
        """
        if len(self.pfam_domains) == 0:
            return False # no PFAM domains, so nothing more 5'
        if most_5prime_hsp.strand == "-":
            most_5prime_hsp = most_5prime_hsp.forward_coordinate_transform()

        # subset PFAM domain is on same frame
        pfam_frames = self.pfam_domains.getdata("frame")
        pfam_same_frame = [seqrng for seqrng in self.pfam_domains if seqrng.frame == frame]

        # take 5'-most PFAM domain. Note these are all on forward strand
        most_5prime_pfam = sorted(pfam_same_frame, key=lambda x: x.start)[0]

        if most_5prime_pfam.start < most_5prime_hsp:
            return most_5prime_pfam
        return None
        
    def predict_orf(self, method='5prime-hsp', use_pfam=True,
                    qs_thresh=16, ss_thresh=40, min_expect=DEFAULT_MIN_EXPECT):
        """
        Predict ORF based on one of two methods:

        1. 5'-most beginning ORF that overlaps 5'-most HSP. This
        procedure errors on the side of too much protein sequence.

        2. ORF starting at the start codon 5' of the 5'-most HSP.

        These are the core two methods for choosing an ORF in the case
        when we:

        - don't suspect missing 5'-end
        - don't suspect a frameshift

        """        
        if not self.has_relative or self.inconsistent_strand(min_expect):
            return None
        
        ## 0. Get strand and anchor HSPs.
        strand = self.get_strand(min_expect)
        most_5prime_relative, most_5prime, most_3prime = self.get_anchor_HSPs(min_expect)

        ## 1. Try to infer frame
        ## 1.a Look for frameshift
        has_majority_frameshift = self.majority_frameshift(min_expect)
        if has_majority_frameshift:
            # Our frame is that of the 5'-most HSP
            frame = most_5prime['frame']
            assert({}[most_5prime.strand] == strand)
        else:
            ## 1.d Finally, infer frame in the vanilla case
            frame = self.majority_frame(min_expect)

        # assert our strand according to strand & frame are consistent
        numeric_strand = {"+":1, "-":-1}[strand]
        assert(numeric_strand == frame/abs(frame))
        
        ## If the frame is negative, we must do a
        ## coordinate transform of the anchor HSPs SeqRange objects so
        ## that they are on the forward orientation (as ORF candidates
        ## would be)
        if frame < 0:
            most_5prime, most_3prime = (most_5prime.forward_coordinate_transform(),
                                        most_3prime.forward_coordinate_transform())

        ## Check for PFAM frames, if necessary
        if use_pfam:
            more_5prime_pfam = self.more_5prime_pfam_domain(most_5prime, frame)
            # TODO annotate PFAM frameshifts
            if more_5prime_pfam is not None:
                most_5prime = more_5prime_pfam

        ## 3. Look for missing 5'-end
        missing_5prime = self.missing_5prime(qs_thresh, ss_thresh, min_expect)

        ## 4. Get all ORFs
        orf_candidates = get_all_orfs(self.record, frame)
        assert(len(orf_candidates) > 0) # we should always have the
                                        # open ended option

        ## 4.a If we don't have a missing 5'-end, remove candidates
        ## that are missing start codon.
        if missing_5prime:
            no_starts = orf_candidates.getdata('no_start')
            tmp = SeqRanges()
            for i, no_start in enumerate(no_starts):
                if not no_start:
                    tmp.append(orf_candidates[i])
            orf_candidates = tmp

        ## 6. ORF Prediction: subset ORFs by those that overlap the
        ## 5'-most HSP
        overlapping_candidates = orf_candidates.subsetByOverlaps(most_5prime)
        if len(overlapping_candidates):
            ## 6.a Method-dependent ORF selection. Method (a): 5'-most
            ## start codon.
            if method == '5prime-most':
                orf_i = range(len(overlapping_candidates))
                tmp = sorted(orf_i, key=lambda x: overlapping_candidates[x].start)
                if len(tmp) > 0:
                    orf_range_i = tmp[0]
                else:
                    return None          
            elif method == '5prime-hsp':
                # which of the overlapping candidates have a start
                # position 5' of the most 5' HSP?
                five_prime_of_hsp_i = filter(lambda i: overlapping_candidates[i].start <= most_5prime.start,
                                             range(len(overlapping_candidates)))
                # let's sort these by start position now
                if len(five_prime_of_hsp_i) > 0:
                    five_prime_of_hsp_i = sorted(five_prime_of_hsp_i,
                                                 key=lambda i: overlapping_candidates[i].start)
                    orf_range_i = five_prime_of_hsp_i[0]
                else:
                    # if no ORF candidates that overlap a 5' HSP have
                    # a start position 5' of the anchor HSP, we take
                    # the 5'-most ORF overlapping candidate and assert
                    # that it's start position is 3' of the 5' HSP
                    # start.
                    orf_i = range(len(overlapping_candidates))
                    tmp = sorted(orf_i, key=lambda x: overlapping_candidates[x].start)
                    orf_range_i = tmp[0]
                    assert(overlapping_candidates[orf_range_i].start > most_5prime.start)
            else:
                raise ValueError("method must be either '5prime-most' or '5prime-hsp'")
        else:
            # no candidates overlap the most 5prime HSP
            return None

        ## 6. Internal stop codon check TODO
        return overlapping_candidates[orf_range_i]
