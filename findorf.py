# findorf.py -- orf prediction and annotation
info = """
findorf.py: ORF prediction and annotation via blastx results of close
  relatives.

findorf.py works by reading in the XML blastx results of mulitple
queries (mRNA contigs) against several databases (run
separately). Each XML blastx results file should correspond to a
particular organism's protein database. With each of these relatives,
the consensus ORF is found, and other attributes of the contig are
added.

There are two primary operations of findorf.py:

1. Join all the XML blastx results with the contig FASTA file.

2. Predict ORFs and annotate contigs based on the data from the join
operation.

These are done seperately, since one may wish to change the parameters
and output from the predict command without having to re-run the join
operation.



TODO

1. Is there a stop codon between the 5'-most HSP and the end of the
sequence?

2. Use the first frame of a frameshifted contig to predict ORF.

3. If there's a stop codon in in from the contig start to the query start, we ignore it.

Strange cases:

k36_contig_9886 starts with a start codon, but has a missing 5'-end.

k21_contig_36350: has a stop codon in first codon and missing 5'-end


k36_contig_9886
k51_contig_10673

k61_contig_20415 - problem detecting frameshift, not in BLAST results.

k26_contig_24653

k26_contig_22146 - frameshift, but orf start/end
"""

import sys
import pdb
import cPickle
import csv
from collections import Counter, namedtuple
from string import Template
from operator import itemgetter, attrgetter
try:
    from Bio.Blast import NCBIXML
    from Bio.Alphabet import IUPAC, generic_dna, DNAAlphabet
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError, e:
    sys.exit("Cannot import BioPython modules; please install it.")
import argparse
import os

from sequtils import STOP_CODONS, START_CODONS, GTF_FIELDS
from sequtils import get_codons

OrfSet = namedtuple('OrfSet', ['start', 'stop', 'length', 'rank'])

def mean(x):
    """
    The arithematic mean.
    """
    return float(sum(x))/len(x) if len(x) > 0 else float('nan')

def get_orfs(seq_in_frame, missing_5prime=False):
    """
    Return the (start, stop) positions for all ORFs. If missing_5prime
    is True, we start from the beginning of the sequence, not the
    first start codon.

    Only ORFs falling under the following conditions are appended to
    the list:

     - full start and stop positions (i.e. both integers)

     - 0 start position, if missing_5prime is True.

     - start position, None stop position (i.e. started open reading
       frame, hit end of sequence before stop codon)

     - None start position, stop position (i.e. a possible ORF with
       missing start).

    Output tuples contain:

     - start position
     - stop position
     - length
     - rank (ordered position in the sequence)
    """
    # Make a list of the codonsb
    seq = str(seq_in_frame).upper()
    codons = get_codons(seq)

    # Initialize values
    start_pos = None if not missing_5prime else 0
    stop_pos = None
    orfs = list()
    rank = 0
    reading = missing_5prime # if we have a missing 5'-end, we assume
                             # we are reading.
    found_any_stop_codon = False

    # for the most part, this loop should function like a ribosome
    # would, except for the possibility that it will add two special
    # cases: missing start codon but hit stop codon, and missing stop
    # codon after hitting a stard codon.
    # Todo: add non-missing 5'-end and 
    for codon, pos in codons:
        if not reading and codon in START_CODONS:
            start_pos = pos
            reading = True
        if not reading and codon in STOP_CODONS and not found_any_stop_codon:
            # we've found our first stop codon, which could indicate
            # an ORF with a stop codon outside of the assembled contig
            orfs.append(OrfSet(None, pos, pos, rank))
            rank += 1
            found_any_stop_codon = True
        if reading and codon in STOP_CODONS and not found_any_stop_codon:
            stop_pos = pos
            orf_length = stop_pos - start_pos
            assert(orf_length >= 0)
            orfs.append(OrfSet(start_pos, stop_pos, orf_length, rank))
            rank += 1
            reading = False
            start_pos = None
            stop_pos = None
            found_any_stop_codon = True
    if reading:
        orf_length = len(seq) - start_pos
        orfs.append(OrfSet(start_pos, stop_pos, orf_length, rank))
    
    return orfs
        
def put_seq_in_frame(seq, frame):
    """
    Take a sequence and transform it to into the correct frame.
    """
    if frame < 0:
        seq = seq.reverse_complement()
        frame = -1*frame
    if not frame in range(1, 4):
        raise Exception, "improper frame: frame must be in [1, 3]"
    return seq[(frame-1):]
    
class ContigSequence():
    """
    A ContigSequence is an assembled contig, that may be coding or
    non-coding.
    """

    def __init__(self, query_id, sequence):
        """
        Initialize a ContigSequence via a BioPython SeqRecord.
        """
        # configurations
        self.e_value_thresh = None
        
        # data attributes
        self.query_id = query_id
        self.all_relatives = dict()
        self.num_hsps = Counter()
        self.seq = sequence
        self.query_length = len(sequence)
        
        # some annotation attributes
        self.missing_start = None
        self.missing_stop = None
        self.full_length_orf = None
        self.orf_start = None
        self.orf_stop = None
        self.orf=None
        self.start_tuples = dict()
        
    def __repr__(self):
        """
        Pretty print all info.
        """

        info = dict(id=self.query_id, length=self.query_length, num_relatives=self.num_relatives,
                    consensus_frame=self.consensus_frame, majority_frame=self.majority_frame,
                    any_frameshift=self.any_frameshift, majority_frameshift=self.majority_frameshift,
                    missing_start=self.missing_start, missing_stop=self.missing_stop,
                    missing_5prime=self.missing_5prime, full_length_orf=self.full_length_orf,
                    orf_start = self.orf_start, orf_stop=self.orf_stop, seq=self.orf)
        
        out = Template("""
ContigSequence element for ID: $id
Length: $length
Number of relatives: $num_relatives

# Frames - these values are in [-3, -2, -1, 1, 2, 3]. GTF/GFF uses [0, 1, 2]
Consensus frame: $consensus_frame
Majority frame: $majority_frame

# Frameshift
Any frameshift: $any_frameshift
Majority frameshift: $majority_frameshift

# ORF Integrity
Missing start codon: $missing_start
Missing stop codon: $missing_stop
5'-end likely missing: $missing_5prime

# Predicted ORF - these values are 0-indexed
ORF is full length: $full_length_orf
ORF start: $orf_start
ORF stop: $orf_stop
ORF seq: $seq

# Relatives Start Sites
""").substitute(info)

        for relative, start_tuple in self.start_tuples.iteritems():
            query_start, sbjct_start, strand = start_tuple
            rel_info = (relative, sbjct_start, query_start, {1:"+", -1:"-"}[strand])
            out += ("%s\n    subject start: %s\n    query start/end"
            " (if strand forward/reverse): %s\n    strand: %s\n" % rel_info)

        # in later versions, we could use a templating engine...
        if self.has_relatives:
            out += "\n# Relative Identities in Frames\n"
            for relative, count_frames in self.frames_identities.iteritems():
                if len(count_frames):
                    out += "%s\n" % relative
                for frame, identities in count_frames.iteritems():
                    out += "  frame: %s\n  identities:  %s\n\n" % (frame, identities)
                    
        return out

    @property
    def relatives(self):
        """
        self.all_relatives is propagated by add_relatives with all
        relatives. findorf predict allows a user to specify an
        e-value, so this property retuns only those items in
        all_relatives with an e-value less than or equal to the
        self.e_value_thresh attribute. If this attribute is None,
        all_relatives is returned.
        """
        
        if self.e_value_thresh is None:
            return self.all_relatives

        passed_thresh = dict()
        for relative, hsps in self.all_relatives.items():
            passed_thresh[relative] = [h for h in hsps if h['e'] <= self.e_value_thresh]

        return passed_thresh
    

    @property
    def num_relatives(self):
        """
        The number of relatives with HSPs.
        """
        num_relatives = len([_ for _, hsps in self.relatives.items() if len(hsps) > 0])
        return num_relatives

    @property
    def max_identities(self):
        """
        Return the relatives, ordered by the number of identities.
        """
        if not self.has_relatives:
            return None

        identities = Counter()
        for relative, hsps in self.relatives.iteritems():
            for h in hsps:
                identities[relative] += h["identities"]

        return identities.most_common()

    def gff_dict(self):
        """
        Return a dictionary of the some key attribute's values, for
        export to a file via the csv module.

        Note that GFFs are 1-indexed, so we add one to positions.
        """
        out = dict()
        out["seqname"] = self.query_id
        out["source"] = "findorf"
        out["feature"] = "predicted_orf"
        out["start"] = self.orf_start + 1 if self.orf_start is not None else "."
        out["end"] = self.orf_stop + 1 if self.orf_stop is not None else "."
        out["score"] = "."
        out["strand"] = self.majority_frame/abs(self.majority_frame) if self.majority_frame is not None else "."
        out["frame"] = abs(self.majority_frame) - 1 if self.majority_frame is not None else "." # GFF uses frames in [0, 2]
        out["group"] = "."
        return out

    def gtf_dict(self):
        """
        Output a GTF file, which carries attributes in the group
        column.
        """

        attributes = dict(full_length_orf=self.full_length_orf,
                          majority_frameshift=self.majority_frameshift,
                          any_frameshift=self.any_frameshift,
                          missing_5prime=self.missing_5prime,
                          number_relatives=len(self.relatives))

        group = "; ".join(["%s %s" % (k, v) for k, v in attributes.iteritems()])
        out = self.gff_dict()
        out["group"] = group
        return out

    def get_hsp_frames(self):
        """
        For each relative, count how the number of occurences of a
        certain frame in HSPs. This creates the attribute
        `frames`. Also, use a counter to keep track of the frames
        accross HSPs and relatives, with the attribute `all_frames`.

        This method sets up the stage for methods like
        consensus_frame, majority_frame, consensus_frameshift,
        majority_frameshift by building up `frames` and `all_frames`.
        """
        self.frames = dict()
        self.all_frames = Counter()

        # for each hsp in a relative, record its frame in both the
        # across-relative set and create a set within each relative.
        for relative, hsps in self.relatives.iteritems():
            self.frames[relative] = Counter()

            # count the number of frames per each relative's HSP
            for h in hsps:
                # TODO check for other tuple elements?
                f = h['frame'][0]
                self.frames[relative][f] += 1
                
                self.all_frames[f] += 1

    @property
    def consensus_frame(self):
        """
        The `consensus_frame` attribute is based on if `all_frames`
        has 1 key (i.e. every relative's HSPs are the same). This is
        the strongest evidence a frame could have (although this is
        dependent on number of HSPs to relatives!). The lack of a
        consensus (due to either differing HSP frames in a relative or
        differing relative hits) leads this to take the value of None.
        """
        if not self.has_relatives:
            return None
        return self.all_frames.keys()[0] if len(self.all_frames) == 1 else None

    @property
    def relative_majority_frame(self):
        """
        The `majority_frame` attribute indicates the majority, based
        on the number of relatives that agree on a frame. We don't
        Note that relatives HSPs in differing frames are *not* counted
        in the majority voting.

        Note that if there is a consensus_frame, majority_frame =
        consensus_frame. We need to do this because a consensus frame
        could have majority ties, which would force the majority to be
        None when there is a consensus. We want majority=None to
        indicate a non-clear majority.
        """
        if not self.has_relatives:
            return None

        if self.consensus_frame is not None:
            return self.consensus_frame
        
        relative_agree_frames = Counter()
        for relative, count_frames in self.frames.iteritems():
            if len(count_frames) == 1: # all HSPs agree (or there's just one)
                only_frame = count_frames.keys()[0]
                relative_agree_frames[only_frame] += 1

        if len(relative_agree_frames):
            # check if there is there a tie, set None if so
            tie = len(set([c for frame, c in relative_agree_frames.most_common(2)])) == 1
            return relative_agree_frames.most_common(1)[0][0] if not tie else None
        else:
            return None

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

        As with `relative_majority_frame`, if `majority_consensus`
        frame is set, this will equal that.
        """
        if not self.has_relatives:
            return None

        if self.consensus_frame is not None:
            return self.consensus_frame
        
        frames = Counter()
        for relative, count_frames in self.frames_identities.iteritems():
            if len(count_frames) == 1: # all HSPs agree (or there's just one)
                only_frame = count_frames.keys()[0]
                frames[only_frame] += count_frames.values()[0]

        if len(frames):
            # check if there is there a tie, set None if so
            tie = len(set([c for frame, c in frames.most_common(2)])) == 1
            return frames.most_common(1)[0][0] if not tie else None
        else:
            return None

    @property
    def relative_majority_frameshift(self):
        """
        Returns True of False indicating if there's a frameshift in
        the majority of relatives. We do not count relatives with one
        HSP as there isn't sufficient information.

        True is returned if there's a tie.
        """
        if not self.has_relatives:
            return None

        frameshifts = Counter()

        for relative, frames in self.frames.iteritems():
            # values corresponds to num HSPs per frame - less than 2,
            # we don't have enough information to infer frameshift.
            if sum(frames.values()) < 2:
                continue
            frameshifts[len(frames.keys()) > 1] += 1

        # handle corner case where there's a tie of zeros (not a
        # frameshift).
        if sum(frameshifts.values()) == 0:
            return False

        return frameshifts[True] >= frameshifts[False]

    @property
    def frames_identities(self):
        """
        Count the number of identities (not HSPs, as in self.frames)
        that are in a relative's HSPS.
        """
        
        if not self.has_relatives:
            return None

        frames = dict()
        for relative, hsps in self.relatives.iteritems():
            frames[relative] = Counter()

            # count the number of identities per each relative's HSP
            for h in hsps:
                # TODO check for other tuple elements?
                f = h['frame'][0]
                frames[relative][f] += h['identities']

        return frames

    @property
    def majority_frameshift(self):
        """
        A refinement of what is now
        relative_majority_frameshift. relative_majority_frameshift has
        problem that hits in a relative that lead to one HSPs are not
        counted in the majority decision. This unfairly penalized HSPs
        that take up the entire query sequence.

        Here, number of identities of each HSP are incorporated, such
        that if the majority of total idenities are in relative's
        alignment with an frameshift, we say it's a majority
        frameshift. This is a finer-version than just looking at the
        majority of relatives with frameshifts because it (1) takes it
        account evolutionary distance implicitly and (2) means we
        don't have to ignore relatives' alignents with only one HSP.

        Note too that we don't have to use proporitions here because
        we're only looking at the best alignment, which does not
        differ in length across the relatives' blast results.
        """
        if not self.has_relatives:
            return None

        frameshifts = Counter()
        for relative, frames in self.frames_identities.iteritems():
            frameshifts[len(frames.keys()) > 1] += sum(frames.values())

        self.frameshift_identities = frameshifts
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
    def consenus_frameshift(self):
        """
        Returns a sorted tuple if all relatives agree on a
        frameshift. Not implemeted now intentionally.
        """
        pass
        

    def get_hsp_start_tuples(self):
        """
        We want to find the 5'-most HSPs. Sorting by query_start is
        not sufficient because reverse-strand cDNA contigs would have
        a minimum query_start value corresponding to the 3'-most
        subject. Note that we're aligning to proteins with blastx, so
        proteins subjects will always be read left to right: start
        translation to stop.

        One could also look minimize subject_start, as the 5'-most HSP
        will have the smallest subject start. However, Ksenia had a
        really nice exception when this logic breaksdown: if there are
        overlapping HSPs, both starting at subject one (due to a
        repeated domain), and one HSP starts at query position 1400
        (assume reverse strand) and the other starts at query position
        1100, then the better ORF would anchored by the 5'-most, so
        minimizing subject start is insufficient.

        So now, we must have a consensus strand before getting the
        5'-most HSP, as we use the minimize query start position. If
        there's a consensus frame, we use its sign; if not, we use the
        majority. If neither exist, we don't calculate this.
        """
        if not self.has_relatives:
            return None

        if self.majority_frame is not None:
            strand = self.majority_frame/abs(self.majority_frame)
            reverse = strand < 0 # higher query_start is the 5' most
                                 # when reverse strand
                                 
            for relative, hsps in self.relatives.iteritems():
                if len(hsps) == 0:
                    continue
                if reverse:
                    most_5prime = sorted(hsps, key=lambda x: x['query_end'], reverse=True)[0]
                    most_5prime_tuple = (most_5prime['query_end'], most_5prime['sbjct_start'], strand)
                else:
                    most_5prime = sorted(hsps, key=lambda x: x['query_start'], reverse=False)[0]
                    most_5prime_tuple = (most_5prime['query_start'], most_5prime['sbjct_start'], strand)
                    
                self.start_tuples[relative] = most_5prime_tuple

    @property
    def missing_5prime(self, qs_thresh=16, ss_thresh=40):
        """
        This attribute indicates where it's predicted that the 5'-end
        of the contig is missing, inferred by whether most relatives
        have an HSP in this region.

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
        

        Defaults are chosen based on some test cases.
        """
        if not self.has_relatives:
            return None

        self.get_hsp_start_tuples()

        missing_5prime = Counter()
        
        for relative, start_tuple in self.start_tuples.iteritems():
            # note that for reverse strand: query_start is really
            # query_end, but we compare it to the difference between
            # query_end and query_length.
            query_start, sbjct_start, strand = start_tuple
            if strand > 0:
                missing_5prime[query_start <= qs_thresh and sbjct_start >= ss_thresh] += 1
            else:
                missing_5prime[abs(query_start - self.query_length) <= qs_thresh and sbjct_start >= ss_thresh] += 1

        if missing_5prime[True] >= missing_5prime[False]:
            return True
        return False

    def find_5prime_most_orf(self, seq_in_frame):
        """
        Find the 5'-most ORF that contains an HSPs.
        """
        orfs = get_orfs(seq_in_frame, missing_5prime=self.missing_5prime)
        self.all_orfs = orfs
        
        if not len(orfs):
            return (None, None, None, None)

        full_orfs = [o for o in orfs if None not in o]
        # Remove cases where there is no HSP in this ORF.
        
        if not len(full_orfs):
            # we return the 5'-most orf
            return sorted(orfs, key=attrgetter('start'))[0]

        return sorted(full_orfs, key=attrgetter('start'))[0]


    def predict_orf(self):
        """
        The `predict_orf` gathers all information to make an ORF
        prediction and annotate this object with information about the
        prediction.

        This involves:

        1. Frame prediction
        
        2. Predicting whether the 5'-end of the sequence is missing
        from presence of HSPs.

        3. Finding the 5'-most ORF that contains HSPs.
        """
        # No relatives? Can't predict ORF
        if not self.has_relatives:
            return None

        frame = self.majority_frame

        # No frame? Can't predict an ORF. TODO handle frameshift?
        if frame is None:
            return None

        seq_in_frame = put_seq_in_frame(self.seq, frame)

        # General note: Be cautious when reading this section: there
        # are two sets of positions: one relative to the sequence once
        # put in frame, the other relative to the raw sequence.

        # relative to sequence in frame; the "_if" refers to in frame
        tmp = self.find_5prime_most_orf(seq_in_frame)
        start_codon_pos_if, stop_codon_pos_if, orf_length, rank = tmp

        # these are booleans indicating whether a valid start of stop
        # found; useful if the orf_start orf_stop positions are set
        # because they're the end of the sequence.
        contains_start = start_codon_pos_if is not None
        contains_stop = stop_codon_pos_if is not None

        # find_longest_orf() works on the sequence already in the
       # frame.  We need to put it back in the coordinates for the
        # original sequence by adding abs(frame-1) (- 1 accounts for
        # difference between 1-indexed BLAST hits and 0-index Python
        # strings). The coordinates are relative to the original
        # sequence. Also, find_longest_orf returns None if it cannot
        # find a start or stop codon in the sequence in frame.
        start_codon_pos = start_codon_pos_if + abs(frame) - 1 if start_codon_pos_if is not None else abs(frame)
        stop_codon_pos = stop_codon_pos_if + abs(frame) - 1 if stop_codon_pos_if is not None else self.query_length

        # these are for slicing the sequence-in-frame
        start_pos_if = start_codon_pos_if if start_codon_pos_if is not None else abs(frame)
        stop_pos_if = stop_codon_pos_if if stop_codon_pos_if is not None else self.query_length
        
        # If start or stop positions found by find_longest_orf are
        # None, this means that we should use the original sequence's
        # start and stop position, adjusted for frame. Since these are
        # sequence-relative, not sequence-in-frame-relative, they make
        # up the attributes.
        self.orf_start = start_codon_pos
        self.orf_stop = stop_codon_pos
        self.orf_pos = (self.orf_start, self.orf_stop)

        if contains_start and contains_stop and not self.missing_5prime:
            self.full_length_orf = True
            self.orf = seq_in_frame[start_pos_if:stop_pos_if]
        elif contains_stop and not contains_start:
            self.orf = seq_in_frame[:stop_pos_if]
        elif not contains_stop and contains_start:
            self.orf = seq_in_frame[start_pos_if:]
        else:
            self.orf = None

        self.missing_start = not contains_start
        self.missing_stop = not contains_stop

        
    def report_first_hsp_starts(self):
        """
        A report method: print out how consistent the start subject
        sites are.
        """
        if not self.has_relatives:
            return None

        info_values = dict()
        msg = Template("""
## ORF Start Prediction
Number of relatives where query's earliest HSP starts at 1: $num_start_1/$num_rels
Number of 
""").substitute()

    @property
    def has_relatives(self):
        """
        Number of relatives > 0?
        """
        return self.num_relatives > 0
        
    def add_relative_alignment(self, relative, blast_record):
        """
        Given a relative and a BioPython BLAST alignment objects,
        extract and store the relevant parts of the _best_ alignment
        only.
        """
        relative_exists = self.all_relatives.get(relative, False)

        if not relative_exists:
            self.all_relatives[relative] = dict()
        else:
            raise Exception, "relative '%s' already exists for this ContigSequence" % relative

        if len(blast_record.alignments) == 0:
            # no alignments, so we dont have any info to add for this
            # relative.
            return 

        self.all_relatives[relative] = list()

        # TODO check: are these gauranteed in best first order?
        best_alignment = blast_record.alignments[0]
        for hsp in best_alignment.hsps:
            hsp_dict = {'align_length':best_alignment.length,
                        'align_title':best_alignment.title,
                        'e':hsp.expect,
                        'identities':hsp.identities,
                        'frame':hsp.frame,
                        'query_start':hsp.query_start,
                        'query_end':hsp.query_end,
                        'sbjct_start':hsp.sbjct_start,
                        'sbjct_end':hsp.sbjct_end}
            
            self.num_hsps[relative] += 1
            self.relatives[relative].append(hsp_dict)
                        

def parse_blastx_args(args):
    """
    Take a list of args and return a dictionary of names and files. If
    any arg is of the key:value sort, the name will be the
    user-defined key; otherwise they will be extracted from the file
    name.
    """
    blastx_files = dict()
    for arg in args:
        tmp = arg.split(":")
        if len(tmp) == 2:
            name, value = tmp
            if blastx_files.get(name, False):
                msg = "key '%s' dy exists in blastx args" % name
                raise argparse.ArgumentTypeError(msg)
            blastx_files[name] = value
        elif len(tmp) == 1:
            blastx_files[os.path.splitext(os.path.basename(arg))[0]] = arg
        else:
            msg = "malformed key-value pair in argument: '%s'" % arg
            raise argparse.ArgumentTypeError(msg)

    # now open all file handles
    try:
        handles = dict([(k, open(f, 'r')) for k, f in blastx_files.items()])
    except IOError, e:
        sys.exit("could not open BLASTX XML result '%s' - no such file.\n" % e.filename)
    return handles


def join_blastx_results(args):
    """
    Each BLAST XML result file contains different alignments, each
    with possibly multiple HSPs. The foreign key is the contig ID,
    which is also the BLAST query ID. This of course is also the key
    between the BLAST results and the reference contig FASTA file.

    This function builds a dictionary which each key being the contig
    ID and each value being another dictionary in which each
    'relative' is the BLAST result file and the values are the BLAST
    results.
    """
    # First we populate the results dictionary with all keys and None
    # value, primarily to keep track of cases where we don't have any
    # BLAST hits.

    results = dict()
    for record in SeqIO.parse(args.ref, "fasta"):
        results[record.id] = ContigSequence(record.id, record.seq)

    blast_files = parse_blastx_args(args.blastx)

    for relative, blast_file in blast_files.items():
        sys.stderr.write("processing BLAST file '%s'...\t" % relative)

        for record in NCBIXML.parse(blast_file):
            query_id = record.query.strip().split()[0]

            # add the relative's alignment information
            results[query_id].add_relative_alignment(relative, record)

        sys.stderr.write("done\n")

    # dump the joined results
    cPickle.dump(results, file=args.output)

def predict_orf(args):
    """
    
    """
    contig_seqs = cPickle.load(args.input)
    num_no_relatives = 0
    
    num_consensus_frames = 0
    num_majority_frames = 0
    num_hsps_different_frames = 0
    num_relatives_majority_frameshift = 0
    
    num_full_length_orf = 0
    num_missing_start = 0
    num_missing_stop = 0

    num_missing_5prime = 0
    total = 0

    # let's sum all the identities to see who comes out on top
    total_identities = Counter()

    for id, cs in contig_seqs.iteritems():
        # We set the e-value threshold based on the argument.
        cs.e_value_thresh = args.e_value
    
        cs.get_hsp_frames()
        cs.get_hsp_start_tuples()
        cs.predict_orf()

        # collect counts of various cases; note we use "is True", as
        # values can be None (default) or False.

        num_no_relatives += int(not cs.has_relatives)
        
        num_consensus_frames += int(cs.consensus_frame is not None)
        num_majority_frames += int(cs.consensus_frame is None and cs.majority_frame is not None)

        num_hsps_different_frames += int(cs.any_frameshift is True)
        num_relatives_majority_frameshift += int(cs.majority_frameshift is True)

        num_full_length_orf += int(cs.full_length_orf is True)
        num_missing_start += int(cs.missing_start is True)
        num_missing_stop += int(cs.missing_stop is True)
        num_missing_5prime += int(cs.missing_5prime is True)

        if cs.max_identities is not None:
            for relative, identities in cs.max_identities:
                total_identities[relative] += identities
        
        total += 1

    print "number *without* relatives (no BLAST hit):", num_no_relatives
    print "number with relatives (has at least one BLAST hit):", total - num_no_relatives
    print "total:", total
    print
    print "number of contigs where all relatives' frames agree:", num_consensus_frames
    print "number of contigs where an identity-weighted majority of relatives' frames agree (but no consensus frame):", num_majority_frames
    print "number of contigs with a possible frameshift (any relative's HSPs are in different frames):", num_hsps_different_frames
    print "number of contigs with a very likely frameshift (majority of identities are in relatives with HSPs in different frames) :", num_relatives_majority_frameshift
    print
    print "number of full length ORFs:", num_full_length_orf
    print "number of missing start codons:", num_missing_start
    print "number of missing stop codons:", num_missing_stop
    print "number of likely missing 5'-ends:", num_missing_5prime
    print "total identities:", ' '.join(["%s: %s" % (rel, count) for rel, count in total_identities.most_common()])
    print

    # GTF and GFF output
    if args.gtf is not None:
        with args.gtf as f:
            dw = csv.DictWriter(f, GTF_FIELDS, delimiter="\t")
            for cs in contig_seqs.values():
                dw.writerow(cs.gtf_dict())

    if args.gff is not None:
        with args.gff as f:
            # same fields as GFF (since group just becomes lumped
            # attributes with GTF)
            dw = csv.DictWriter(f, GTF_FIELDS, delimiter="\t")
            for cs in contig_seqs.values():
                dw.writerow(cs.gff_dict())

    if args.dense is not None:
        with args.dense as f:
            for cs in contig_seqs.values():
                if cs.has_relatives:
                    f.write("----------------%s" % str(cs))


    return contig_seqs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=info)
    subparsers = parser.add_subparsers(help="sub-commands")

    ## join arguments
    parser_join = subparsers.add_parser('join', help="join each BLASTX XML results file")
    parser_join.add_argument('blastx', type=str, 
                             nargs="+", help="the output XML of BLASTx")
    parser_join.add_argument('--output', type=argparse.FileType('wb'),
                             default="joined_blastx_dbs.pkl",
                             help="joined results of all the BLASTX results files (Python pickle file)")
    parser_join.add_argument('--ref', type=str, required=False,
                             help="the FASTA reference that corresponds to BLASTX queries.")

    parser_join.set_defaults(func=join_blastx_results)

    ## predict arguments
    parser_predict = subparsers.add_parser('predict', help="predict ORFs")
    parser_predict.add_argument('--input', type=argparse.FileType('rb'),
                                default="joined_blastx_dbs.pkl",
                                help="the joined results of all BLASTX files (Python pickle file; from sub-command joined)")
    parser_predict.add_argument('--gff', type=argparse.FileType('w'),
                                default=None,
                                help="filename of the GFF of ORF predictions")
    parser_predict.add_argument('--gtf', type=argparse.FileType('w'),
                                default=None,
                                help="filename of the GTF of ORF predictions")
    parser_predict.add_argument('--dense', type=argparse.FileType('w'),
                                default=None,
                                help="filename of the dense output file")
    parser_predict.add_argument('--full-length', action="store_true",
                                help="the FASTA reference that corresponds to BLASTX queries.")
    parser_predict.add_argument('-e', '--e-value', type=float,
                                default=10e-3,
                                help="e-value threshold (relative hit only include if less than this)")

    parser_predict.set_defaults(func=predict_orf)

    args = parser.parse_args()

    ## Argument checking
    results = args.func(args)

    
    
