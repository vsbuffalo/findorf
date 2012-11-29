"""
contig.py contains the Contig class. This new version is a depature
from the previous version in that (1) it relies upon BioRanges and (2)
it stores uses objects from BioPython directly.

Strange cases:
 - k61_contig_13995

TODO add stop codon number of BP
"""

import pdb
from collections import Counter, namedtuple
from operator import itemgetter
from itertools import groupby

try:
    from Bio.Alphabet import IUPAC
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError:
    sys.exit("Cannot import BioPython modules; please install it.")

from BioRanges.lightweight import Range, SeqRange, SeqRanges
from orfprediction import get_all_orfs, ORFTypes, count_5prime_ATG

# named tuples used to improve readability 
AnchorHSPs = namedtuple('AnchorHSPs', ['relative', 'most_5prime', 'most_3prime'])

NUCLEOTIDES = IUPAC.IUPACAmbiguousDNA.letters
MASK_CHAR = "X"
DEFAULT_MIN_EXPECT = 10
ANNOTATION_FIELDS = ["most_5prime_relative", # string
                     "pfam_extended_5prime", # boolean
                     "num_orf_candidates", # integer
                     "contig_len", # integer
                     "num_relatives", # integer
                     "majority_frameshift", # boolean
                     "internal_stop", # boolean
                     "distant_start", # boolean
                     "most_5prime_query_start", # integer
                     "most_5prime_sbjct_start", # integer
                     "diff_5prime_most_start_and_orf", # integer
                     "num_5prime_ATG"] # integer 


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
        self.annotation = dict().fromkeys(ANNOTATION_FIELDS)
        self.annotation["contig_len"] = len(self.record.seq)
        self.hsps = SeqRanges()
        self.pfam_domains = SeqRanges()
        self.has_relative = False
        self.orf = None
        self.orf_candidates = None
        self.orf_type = None
        
    @property
    def seq(self):
        """
        Return the sequence of the contig.
        """
        return self.record.seq
    
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

        if self.orf is not None:
            out["start"] = self.orf.start + 1
            out["end"] = self.orf.end + 1
        else:
            out["start"] = "."
            out["end"] = "."
            
        out["score"] = "."

        maj_frame = self.orf['frame'] if self.orf is not None else None
        if maj_frame is not None:
            out["strand"] = maj_frame/abs(maj_frame)
        else:
            out["strand"] = "."
            
        if maj_frame is not None:
            # GFF uses frames in [0, 2]
            out["frame"] = abs(maj_frame) - 1
        else:
            out["frame"] = "."
            out["group"] = "."
        return out
        
    def gtf_dict(self):
        """
        Return a dictionary corresponding to the columns of a GTF
        file.
        """
        orf_anno = {"orf_type":self.orf_type.type,
                    "no_prediction_reason":self.orf_type.reason}
        anno = dict(self.annotation.items() + orf_anno.items())
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
            seq = self.orf_seq
            # we get this from ORF so we don't have to re-look at min expect
            frame = self.orf["frame"]
            desc = self.description + " translated from frame %s" % str(frame)
            return SeqRecord(seq=seq.seq.translate(), id=self.id,
                             description=self.annotated_description)
        return None

    @property
    def orf_seq(self):
        """
        Return the nucleotide sequence record
        """
        if self.orf is not None:
            seq = self.seq
            if self.orf['frame'] < 0:
                seq = seq.reverse_complement()
            seq = self.orf.sliceseq(seq)
            return SeqRecord(seq=seq, id=self.id, description=self.annotated_description)
        return None

    @property
    def orf_masked(self):
        """
        Return the original contig sequence (as BioPython SeqRecord) with the predicted ORF
        masked.
        """
        if self.orf is None:
            return self.record        
        if self.orf['frame'] < 0:
            seq = self.orf.maskseq(self.seq.reverse_complement(), MASK_CHAR)
        else:
            seq = self.orf.maskseq(self.seq)

        # let's put in an assertion here that we're not losing any
        # sequence. Unit tests (tests/test_contig.py) cover some this
        # too.
        assert(seq.count(MASK_CHAR) == self.orf.width)
        return SeqRecord(seq=seq, id=self.id, description=self.description)
        
    @property
    def id(self):
        """
        Return the sequence header ID.
        """
        return self.record.id

    @property
    def annotated_description(self):
        """
        Return a longer, annotated description of any ORFs found. 
        """
        if self.orf is None:
            return self.description
        hsp_id = self.orf["most_5prime_hsp"]["title"].split(" ")[0]
        pfam_extension = self.annotation["pfam_extended_5prime"]
        internal_stop = self.annotation["internal_stop"]
        majority_frameshift = self.annotation["majority_frameshift"]
        if not internal_stop and not majority_frameshift:
            msg = "predicted ORF (type '%s') based on protein '%s' from relative '%s'"
        else:
            if internal_stop and not majority_frameshift:
                pg_type = "contains premature stop codon"
            elif internal_stop and majority_frameshift:
                pg_type = "majority frameshift and contains premature stop codon"
            elif majority_frameshift:
                pg_type = "majority frameshift"
            else:
                raise ValueError
            msg = "predicted ORF (type '%s', likely " + pg_type + ") based on protein '%s' from relative '%s'"
        msg = msg % (self.orf_type.type, hsp_id, self.orf["most_5prime_hsp"]["relative"])
        if pfam_extension:
            msg += " with PFAM domain extension"
        return msg + "; " + " ".join(self.description.split(" ")[1:])


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

        Note that there is no requirement these have the same
        relative. Any homologous sequence from any relative can set
        the 5' or 3' HSP. This is important to note because when we do
        internal stop codon checks, we could end up with the situation
        that an HSP is in a masked region (i.e. due to transposable
        elements), and this is *not* evidence of a internal stop codon.
        """
        if not self.has_relative:
            return None

        filtered_hsps = filter(lambda x: x['expect'] <= min_expect, self.hsps)
        if len(filtered_hsps) < 1:
            return None

        # life is easier if we turn these back into a
        # SeqRange. Eventually, BioRanges.lightweight.SeqRanges should
        # have a method for this.
        tmp = SeqRanges()
        for fhsp in filtered_hsps:
            tmp.append(fhsp)
        filtered_hsps = tmp

        # We get the outermost HSPs
        i = sorted(range(len(filtered_hsps)), key=lambda k: filtered_hsps.end[k], reverse=True)[0]
        j = sorted(range(len(filtered_hsps)), key=lambda k: filtered_hsps.start[k])[0]

        if self.get_strand(min_expect) == "-":
            # negative strand; 5'-most HSP is that with the largest
            # query end
            return AnchorHSPs(filtered_hsps[i]['relative'], filtered_hsps[i], filtered_hsps[j])
        else:
            # positive strand; 5-most HSP is that with the smallest
            # query start
            return AnchorHSPs(filtered_hsps[j]['relative'], filtered_hsps[j], filtered_hsps[i])

    def get_strand(self, min_expect=DEFAULT_MIN_EXPECT):
        """
        Get a strand (+, -), a step we can do before guess frame.

        We need strand to infer 5'-anchor HSPs, which we need if we
        have a frameshift, so this must be found before frame.
        """
        if not self.has_relative or self.inconsistent_strand(min_expect):
            return None

        filtered_hsps = filter(lambda x: x['expect'] <= min_expect, self.hsps)
        if len(filtered_hsps) < 1:
            return None
        
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
        if len(filtered_hsps) < 1:
            return None
        
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
        if len(filtered_hsps) < 1:
            return None

        return len(set([seqrng["frame"]/abs(seqrng["frame"]) for seqrng in filtered_hsps])) > 1

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
        if len(filtered_hsps) < 1:
            return None
        
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

    def add_pfam(self, domain_hit_seqrange):
        """
        Add PFAM domain hit (from HMMER). Note that all of the
        coordinate conversion is done via add_pfam_domain_hits()
        function in the hmmer module.
        """
        self.pfam_domains.append(domain_hit_seqrange)

    def internal_stop_codon(self, orf_end, orf_frame):
        """
        Check if there are any _non-masked_ HSPs more 3' than the ORF
        end position (everything on forward strand)

        DEPRECATED
        """
        if not self.has_relative or self.orf is None:
            return None
        masked_letters = NUCLEOTIDES.lower()
        for hsp in self.hsps:
            if orf_frame < 0:
                hsp = hsp.forward_coordinate_transform()
                seq = str(self.seq.reverse_complement())
            else:
                seq = str(self.seq)
            contains_masked = any(letter in masked_letters for letter in hsp.sliceseq(seq))
            if hsp.start > orf_end and not contains_masked:
                return True
        return False

    def majority_internal_stop(self, buffer_bp=60, min_expect=DEFAULT_MIN_EXPECT):
        """
        A more conservative internal stop codon detection approach. In
        this case, much like we do when looking at majority
        frameshifts, we just consider HSPs grouped by protein.

        buffer_bp is the threshold by which the end must pass the end
        of the ORF:

                      ORF end
        ------------------|   buffer_bp
        ---------------------------|--------| HSP end
        
        This function looks at cases where there is an HSP that
        overlaps the ORF, but has an ending that satifies the criteria
        above.

        We use strictly greater than (>) here rather than >= (as we do
        in majority_frameshift) because we're comparing number of
        relatives rather than number of identities.
        """
        if not self.has_relative or self.orf is None or self.inconsistent_strand(min_expect):
            return None

        filtered_hsps = filter(lambda x: x['expect'] <= min_expect, self.hsps)
        if len(filtered_hsps) < 1:
            return None

        # we join title with relatives. title should be unique by
        # relative, but why make assumptions?
        name_join = lambda x, y: "%s-%s" % (x, y)

        # combine HSPs by alignment/relative. Currently findorf uses
        # only the top alignment
        filtered_hsps = [(name_join(h['relative'], h['title']), h) for h in filtered_hsps]
        filtered_hsps = sorted(filtered_hsps, key=itemgetter(0))        

        internal_stop = Counter()
        for alignment, hsps_grouped in groupby(filtered_hsps, itemgetter(0)):
            hsps_grouped = list(hsps_grouped)
            hsps = map(itemgetter(1), hsps_grouped)

            # at this point, we're guaranted consistent strand, so if
            # the first has negative frame, they all do and need to be
            # transformed.
            if hsps[0]['frame'] < 0:
                hsps = [h.forward_coordinate_transform() for h in hsps]

            # if any HSP of a protein overlaps the ORF, we consider
            # its end position.
            if any(self.orf.overlaps(h) for h in hsps):
                max_end = max(h.end for h in hsps)
                if max_end > self.orf.end + buffer_bp:
                    internal_stop[True] += 1
                    continue
            internal_stop[False] += 1

        return internal_stop[True] > internal_stop[False]

    def more_5prime_pfam_domain(self, most_5prime_hsp, frame, min_expect=DEFAULT_MIN_EXPECT):
        """
        Return PFAM domain more 5' prime of supplied SeqRange object
        (which should be the 5' anchor HSP), or None of if there is none.

        Note that all PFAM domains are on the positive strand, since
        PFAM domains found via CAN were in protein space.
        """
        if len(self.pfam_domains) == 0:
            return None # no PFAM domains, so nothing more 5'
        if most_5prime_hsp.strand == "-":
            most_5prime_hsp = most_5prime_hsp.forward_coordinate_transform()

        # subset PFAM domain is on same frame
        pfam_frames = self.pfam_domains.getdata("frame")
        pfam_same_frame = [seqrng for seqrng in self.pfam_domains if seqrng['frame'] == frame]

        # take 5'-most PFAM domain. Note these are all on forward strand
        most_5prime_pfams = sorted(pfam_same_frame, key=lambda x: x.start)

        if len(most_5prime_pfams) > 0 and most_5prime_pfams[0].start < most_5prime_hsp.start:
            return most_5prime_pfams[0]
        return None
        
    def predict_orf(self, method='5prime-hsp', use_pfam=True, min_expect=DEFAULT_MIN_EXPECT):
        """
        Predict ORF based on one of two methods:

        1. 5'-most beginning ORF that overlaps 5'-most HSP. This
        procedure errors on the side of too much protein sequence.

        2. ORF starting at the start codon 5' of the 5'-most HSP.

        These are the core two methods for choosing an ORF in the case
        when we:

        - don't suspect missing 5'-end
        - don't suspect a frameshift

        TODO: this is a huge method; in the future this should be
        refactored and maybe put in a new module.
        """        
        if not self.has_relative:
            self.orf_type = ORFTypes(None, "no_relative")
            return None
        if self.inconsistent_strand(min_expect):
            self.orf_type = ORFTypes(None, "inconsistent_strand")
            return None

        # even though every function does this, we do it here to
        # return None if none pass thresholds.
        filtered_hsps = filter(lambda x: x['expect'] <= min_expect, self.hsps)
        if len(filtered_hsps) < 1:
            self.orf_type = ORFTypes(None, "none_passed_expect_thresh")
            return None
        self.annotation["num_relatives"] = len(set([s['relative'] for s in filtered_hsps]))

        ## 0. Get strand and anchor HSPs.
        strand = self.get_strand(min_expect)
        most_5prime_relative, most_5prime, most_3prime = self.get_anchor_HSPs(min_expect)
        self.annotation['most_5prime_relative'] = most_5prime_relative

        ## 1. Try to infer frame
        ## 1.a Look for frameshift
        has_majority_frameshift = self.majority_frameshift(min_expect)
        self.annotation["majority_frameshift"] = has_majority_frameshift
        if has_majority_frameshift:
            # Our frame is that of the 5'-most HSP
            frame = most_5prime['frame']
        else:
            ## 1.d Finally, infer frame in the vanilla case
            frame = self.majority_frame(min_expect)

        # assert our strand according to strand & frame are consistent
        numeric_strand = {"+":1, "-":-1}[strand]
        assert(int(numeric_strand) == int(frame/abs(frame)))
        
        ## If the frame is negative, we must do a
        ## coordinate transform of the anchor HSPs SeqRange objects so
        ## that they are on the forward orientation (as ORF candidates
        ## would be)
        if frame < 0:
            most_5prime, most_3prime = (most_5prime.forward_coordinate_transform(),
                                        most_3prime.forward_coordinate_transform())

        most_5prime_hsp = most_5prime # reference for annotation, in case of PFAM extension

        ## Check for PFAM frames, if necessary
        if use_pfam:
            more_5prime_pfam = self.more_5prime_pfam_domain(most_5prime, frame)
            if more_5prime_pfam is not None:
                most_5prime = more_5prime_pfam
            self.annotation["pfam_extended_5prime"] = more_5prime_pfam is not None
            
        ## 4. Get all ORFs
        orf_candidates = get_all_orfs(self.record, frame)
        self.orf_candidates = orf_candidates
        self.annotation["num_orf_candidates"] = len(orf_candidates)
        if len(orf_candidates) == 0:
            # why would we have no ORF candidate at all? usually there
            # should be the open-ended case. However, if a sequence's
            # first codon is a stop codon and no start codons are
            # found, there can be no ORF.
            self.orf_type = ORFTypes(None, "no_orf_candidates")
            return None

        ## 6. ORF Prediction: subset ORFs by those that overlap the
        ## 5'-most HSP
        overlapping_candidates = orf_candidates.subsetByOverlaps(most_5prime)
        if len(overlapping_candidates):
            ## 6.a Method-dependent ORF selection. Method (a): 5'-most
            ## start codon. If there is none, we take the open-ended
            ## case.
            if method == '5prime-most':
                orf_i = range(len(overlapping_candidates))
                tmp = sorted(orf_i, key=lambda x: overlapping_candidates[x].start)
                assert(len(tmp) > 0)
                orf_range_i = tmp[0]
                # assert(not overlapping_candidates[orf_range_i]["no_start"])
            elif method == '5prime-hsp':
                # which of the overlapping candidates have a start
                # position 5' of the most 5' HSP?
                five_prime_of_hsp_i = filter(lambda i: overlapping_candidates[i].start <= most_5prime.start,
                                             range(len(overlapping_candidates)))
                # let's sort these by start position now, reversing so
                # that the latest ORF candidate that overlaps is chosen
                if len(five_prime_of_hsp_i) > 0:
                    five_prime_of_hsp_i = sorted(five_prime_of_hsp_i,
                                                 key=lambda i: overlapping_candidates[i].start,
                                                 reverse=True)
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
            self.orf_type = ORFTypes(None, "no_overlap")
            return None
        orf = overlapping_candidates[orf_range_i]
        self.orf = orf
        self.orf["frame"] = frame
        self.orf["most_5prime_hsp"] = most_5prime_hsp
        
        # check for ORF type, and annotate
        self.orf_type = ORFTypes(self.orf)

        ## 6. Internal stop codon check
        self.annotation["internal_stop"] = self.majority_internal_stop()

        ## 7. Annotate other 5' start sites.
        if orf is not None:
            self.annotation["num_5prime_ATG"] = count_5prime_ATG(self.seq, frame, orf.start)

        ## 8. Annotate the furthest 5 ORF candidate start position's
        ## difference with current orf start position (ignoring open
        ## ended cases)
        if self.orf is not None:
            starts = [x.start for x in orf_candidates if not x["no_start"] and orf.start > x.start]
            if len(starts) == 0:
                self.annotation["diff_5prime_most_start_and_orf"] = 0
            else:
                tmp = orf.start - max(starts)
                assert(tmp > 0)
                self.annotation["diff_5prime_most_start_and_orf"] = tmp

        ## Annotate the data used in the 5'-most HSP, specifically
        ## subject and query start
        self.annotation["most_5prime_query_start"] = most_5prime_hsp.start
        self.annotation["most_5prime_sbjct_start"] = most_5prime_hsp['sbjct_start']
        return orf
    
if __name__ == "__main__":
    # code for testing this function interactively
    import cPickle
    a = cPickle.load(open("joined_blastx_dbs.pkl"))

    # for key, contig in a.iteritems():
    #     anch = contig.get_anchor_HSPs()
    #     frame = contig.majority_frame()
    #     fs = contig.majority_frameshift()
    #     if frame is None or fs:
    #         continue
    #     # get most 5' PFAM
    #     tmp = anch.most_5prime
    #     if frame < 0:
    #         tmp = anch.most_5prime.forward_coordinate_transform()
        
    #     more_5prime_pfam = contig.more_5prime_pfam_domain(tmp, frame)
    #     if more_5prime_pfam is None:
    #         continue
    #     if not anch.most_5prime.overlaps(more_5prime_pfam):
    #         print key
