## findorf.py -- orf prediction and annotation
info = """
findorf.py: ORF prediction and annotation via blastx results of close
  relatives.

findorf.py works by reading in the XML blastx results of mulitple
queries (mRNA contigs) against several databases (run
separately). Each XML blastx results file should correspond to a
particular organism's protein database. With each of these relatives,
the consensus ORF is found, and other attributes of the contig are
added.

Specifying blastx Results


TODO: is there a stop codon betwene the 5'-most HSP and the end of the
sequence?

k51_contig_10673

k61_contig_20415 - problem detecting frameshift, not in BLAST results.

k26_contig_24653

k26_contig_22146 - frameshift, but orf start/end
"""

STOP_CODONS = ("TAG", "TGA", "TAA")
START_CODONS = ("ATG")
GTF_FIELDS = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "group"]

import sys
import pdb
import cPickle
import csv
from collections import Counter
from string import Template
from operator import itemgetter
try:
    from Bio.Blast import NCBIXML
    from Bio.Alphabet import IUPAC, generic_dna, DNAAlphabet
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError, e:
    sys.exit("Cannot import BioPython modules; please install it.")

import argparse
import os

def mean(x):
    """
    The arithematic mean.
    """
    return float(sum(x))/len(x) if len(x) > 0 else float('nan')

def find_longest_orf(seq_in_frame, missing_start=False):
    """
    Return the (start, stop) positions for the longest ORF, or None if
    one is not found.

    Note that longest ORF will *not* (of course) be from the earliest
    M to the latest stop codon, since this is not a proper ORF.

    TODO: push each ORF on to a stack, take the longest.
    """
    seq = str(seq_in_frame)

    # these are in position order - TODO use namedtuple?
    codons = [(seq[pos:pos+3], pos) for pos in range(0, len(seq), 3)]
    start_pos = None if missing_start is False else 0
    stop_pos = None

    orfs = list()
    for codon in codons:
        if start_pos is None and codon[0] in START_CODONS:
            start_pos = codon[1]
        # note that we require start_pos is not None... we've must
        # have already found a start codon
        if stop_pos is None and codon[0] in STOP_CODONS and start_pos is not None:
            # note that we subract 1 because we don't want beginning
            # to first position of stop codon, but rather last
            # nucleotide in the ORF
            stop_pos = codon[1] - 1
        if start_pos is not None and stop_pos is not None:
            orfs.append((start_pos, stop_pos, abs(stop_pos - start_pos)))
            start_pos = None
            stop_pos = None
    #pdb.set_trace()
    if not len(orfs):
        return (None, None, None)
    return sorted(orfs, key=itemgetter(2), reverse=True)[0]

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
        # data attributes
        self.query_id = query_id
        self.relatives = dict()
        self.has_relatives = False
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

        info = dict(id=self.query_id, length=self.query_length, num_relatives=len(self.relatives),
                    consensus_frame=self.consensus_frame, majority_frame=self.majority_frame,
                    any_frameshift=self.any_frameshift, majority_frameshift=self.majority_frameshift,
                    missing_start=self.missing_start, missing_stop=self.missing_stop,
                    likely_missing_5prime=self.likely_missing_5prime, full_length_orf=self.full_length_orf,
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
5'-end likely missing: $likely_missing_5prime

# Predicted ORF - these values are 0-indexed
ORF is full length: $full_length_orf
ORF start: $orf_start
ORF stop: $orf_stop
ORF seq: $seq

# Relatives Start Sites
""").substitute(info)

        for relative, start_tuple in self.start_tuples.iteritems():
            sbjct_start, query_start, strand = start_tuple
            rel_info = (relative, sbjct_start, query_start, {1:"+", -1:"-"}[strand])
            out += ("%s\n    subject start: %s\n    query start/end"
            " (if strand forward/reverse): %s\n    strand: %s\n" % rel_info)
    
        return out


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
                          likely_missing_5prime=self.likely_missing_5prime,
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
    def majority_frame(self):
        """
        The `majority_frame` attribute indicates the majority, based
        on the number of relatives that agree on a frame. We don't
        take the number of HSPs, as one long, long HSP shouldn't be
        counted less than many tiny HSPs (in the future, we may look
        at number of identities in a frame). Majority across relatives
        seems like the best approach. Note that relatives HSPs in
        differing frames are *not* counted in the majority voting.

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
    def old_majority_frameshift(self):
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
        old_majority_frameshift. old_majority_frameshift has problem
        that hits in a relative that lead to one HSPs are not counted
        in the majority decision. This unfairly penalized HSPs that
        take up the entire query sequence.

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
                if reverse:
                    most_5prime = sorted(hsps, key=lambda x: x['query_end'], reverse=True)[0]
                    most_5prime_tuple = (most_5prime['query_end'], most_5prime['sbjct_start'], strand)
                else:
                    most_5prime = sorted(hsps, key=lambda x: x['query_start'], reverse=False)[0]
                    most_5prime_tuple = (most_5prime['query_start'], most_5prime['sbjct_start'], strand)
                    
                self.start_tuples[relative] = most_5prime_tuple

    @property
    def likely_missing_5prime(self, qs_thresh=16, ss_thresh=40):
        """
        This is an important function: we look at the query/subject
        start positions of the 5'most HSP. If the subject start is
        late (and the query start is early), it probably means that
        this contig is missing a 5'-end.

        This is based on a majority count procedure. To be
        conservative, ties go to "yes" - that the 5 prime is missing.

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


    def predict_orf(self):
        """
        Predict the ORF from the consensus or majority frame's largest
        ORF.

        This also checks for a stop codon (valid 3'-end) and start
        codon.
        """
        if not self.has_relatives:
            return None

        frame = self.consensus_frame if self.consensus_frame is not None else self.majority_frame

        # if we can't predict the frame, we can't predict an ORF.
        if frame is None:
            return None

        seq = put_seq_in_frame(self.seq, frame)


        # General note: Be cautious when reading this section: there
        # are two sets of positions: one relative to the sequence once
        # put in frame, the other relative to the raw sequence.

        # relative to sequence in frame; the "_if" refers to in frame
        start_codon_pos_if, stop_codon_pos_if, orf_length = find_longest_orf(seq, missing_start=self.likely_missing_5prime)

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
        start_codon_pos = start_codon_pos_if + abs(frame) - 1 if start_codon_pos_if is not None else frame
        stop_codon_pos = stop_codon_pos_if + abs(frame) - 1 if stop_codon_pos_if is not None else self.query_length

        # these are for slicing the sequence-in-frame
        start_pos_if = start_codon_pos_if if start_codon_pos_if is not None else frame
        stop_pos_if = stop_codon_pos_if if stop_codon_pos_if is not None else self.query_length
        
        # If start or stop positions found by find_longest_orf are
        # None, this means that we should use the original sequence's
        # start and stop position, adjusted for frame. Since these are
        # sequence-relative, not sequence-in-frame-relative, they make
        # up the attributes.
        self.orf_start = start_codon_pos
        self.orf_stop = stop_codon_pos
        self.orf_pos = (self.orf_start, self.orf_stop)

        if contains_start and contains_stop and not self.likely_missing_5prime:
            self.full_length_orf = True
            self.orf = seq[start_pos_if:stop_pos_if]
            self.missing_start = False
            self.missing_stop = False
        elif contains_stop and not contains_start:
            self.missing_start = True
            self.missing_stop = False
            self.orf = seq[:stop_pos_if]
        elif not contains_stop and contains_start:
            self.missing_stop = True
            self.missing_start = Fale
            self.orf = seq[start_pos_if:]
        else:
            self.orf = None
            self.missing_start = True
            self.missing_stop = True

        
    def report_first_hsp_starts(self):
        """
        A report method: print out how consistent the start subject
        sites are.
        """
        if not self.has_relatives:
            return None

        info_values = {}
        msg = Template("""
## ORF Start Prediction
Number of relatives where query's earliest HSP starts at 1: $num_start_1/$num_rels
Number of 
""").substitute()
        
    def add_relative_alignment(self, relative, blast_record):
        """
        Given a relative and a BioPython BLAST alignment objects,
        extract and store the relevant parts of the _best_ alignment
        only.
        """
        relative_exists = self.relatives.get(relative, False)

        if not relative_exists:
            self.relatives[relative] = dict()
        else:
            raise Exception, "relative '%s' already exists for this ContigSequence" % relative

        self.has_relatives = True
        

        # TODO check: are these gauranteed in best first order?
        self.relatives[relative] = list()
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
                msg = "key '%s' already exists in blastx args" % name
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

    num_likely_missing_5prime = 0
    total = 0
    
    for relative, cs in contig_seqs.iteritems():
        #pdb.set_trace()
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
        num_likely_missing_5prime += int(cs.likely_missing_5prime is True)
        
        total += 1

    print "number of contigs with all relatives' frames agree:", num_consensus_frames
    print "number of contigs where a majority of relatives' frames agree (but no consensus frame):", num_majority_frames
    print "number of contigs with a relative with HSPs in different frames:", num_hsps_different_frames
    print "number of contigs with a frameshift in the majority of relatives:", num_relatives_majority_frameshift
    print
    print "number of full length ORFs:", num_full_length_orf
    print "number of missing start codons:", num_missing_start
    print "number of missing stop codons:", num_missing_stop
    print "number of likely missing 5'-ends:", num_likely_missing_5prime
    print
    print "number *without* relatives (no BLAST hit):", num_no_relatives
    print "number with relatives (has at least BLAST hit):", total - num_no_relatives
    print "total:", total

    if args.gtf is not None:
        with args.gtf as f:
            dw = csv.DictWriter(f, GTF_FIELDS, delimiter="\t")
            for cs in contig_seqs.values():
                dw.writerow(cs.gtf_dict())

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
    parser_predict.add_argument('--full-length', action="store_true",
                                help="the FASTA reference that corresponds to BLASTX queries.")

    parser_predict.set_defaults(func=predict_orf)

    args = parser.parse_args()

    ## Argument checking
    results = args.func(args)

    
    
