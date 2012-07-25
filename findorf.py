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

"""

import sys
import pdb
import cPickle
from collections import Counter
from string import Template
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
    
def pp_dict_of_counters(x):
    """
    Pretty print a dict of counters.
    """
    out = ""
    for key, counter in x.iteritems():
        for value, count in counter.iteritems():
            out += "%s\t%s\t%s\n" % (key, value, count)
    return out

class ContigSequence():
    """
    A ContigSequence is an assembled contig, that may be coding or
    non-coding.
    """

    def __init__(self):
        """
        Initialize a ContigSequence via a BioPython SeqRecord.
        """
        self.relatives = dict()
        self.num_hsps = Counter()

        ## annotation attributes
        self.annotation = dict(full_length_orf=None, missing5prime=None,
                               missing3prime=None, likely_psuedogene=None)
        
    def get_hsp_frames(self):
        """
        For each relative, count how the number of occurences of a
        certain frame in HSPs. This creates the attribute
        `frames`. Also, use a counter to keep track of the frames
        accross HSPs and relatives, with the attribute `all_frames`.
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

    def all_relatives_agree_frame(self):
        """
        Do all relatives agree on the same frame?
        """
        return len(self.all_frames) == 1

    def all_relatives_hsps_agree_frame(self):
        """
        Does every HSP's frame agree within a relative?

        This has the side-effect that it sets the attribute
        `relatives_with_diff_hsp_frames`.

        Note that if all_relatives_agree_frame() is False, this
        necessarily must be so too.
        """
        relatives_with_diff_hsp_frames = [r for r, c in self.frames.items() if len(c.keys()) > 1]
        relatives_consistent = True if len(relatives_with_diff_hsp_frames) == 0 else False

        self.relatives_with_diff_hsp_frames = relatives_with_diff_hsp_frames

        return relatives_consistent


#     def repr_hsp_annotation(self):

#         """
        
#         """

#         info_values = {'num_hsps':repr(dict(self.num_hsps))[1:-1],
#                        'all_frames':', '.join(list(self.all_frames)),
#                        'frame_summary': pp_dict_of_counters(self.frames)}
        
#         msg = Template("""
# ## Frameshifts
# Number of HSPs per relative: $num_hsps
# Frames across all relatives: $all_frames
# Frames by relative:
# \tframe\tnumber of HSPs
# $frame_summary
# """).substitute(info_values) # TODO consensus
#         print msg
        

    def get_hsp_start_tuples(self):
        """
        For each relative, store the tuples (earliest HSP start,
        subject start). Ideally, this would be (1, 1). Cases of (1, x)
        where x > 1 indicates alignments of the earleist HSP where:

                 |----------------------------- query
        |-------------------------------------- subject

        Indicating that the contig is probably missing a 5'end. But
        this annotation is left for another method.
        """
        self.start_tuples = dict()
        for relative, hsps in self.relatives.iteritems():
            # first HSP by *query* position (TODO: mind strandedness)
            first_hsp = sorted(hsps, key=lambda x: x['query_start'])[0]
            start_tuple = (first_hsp['query_start'], first_hsp['sbjct_start'])
            self.start_tuples[relative] = start_tuples


    def all_relative_same_query_start(self):
        """
        Do all relatives have the exact same query start? This
        indicates consistency of the query.
        """
        return len(set([qs for qs, ss in self.start_tuples])) == 1

    def all_relaties_fuzzy_query_start(self, fuzzy=10):
        """
        Do all relatives have roughly the same query start? This is
        calculated as their abs(mean - x) < fuzzy.

        Note this does not discount a different frame... this may be a
        TOADD feature. For example, if all the relatives agree that
        the query is in the 3rd frame, this leaves less room for the
        fuzzy to kick in.
        """
        qs_mean = mean([qs for qs, ss in self.start_tuples])
        return all([abs(qs - qs_mean) < fuzzy for qs, ss, in self.start_tuples])

    def all_relatives_have_earliest_possible_start(self):
        """
        Do all query starts agree (same position), and is this
        position the earliest possible start (given the frame).

        Note this does not exclude the case that the subject start is
        much after the query start, possibly incidating a missing
        5'-end of a contig
        """
        if not self.all_relatives_agree_frame():
            return False
        frame = list(self.all_frames)[0]

        query_start = list(set([qs for qs, ss in self.start_tuples]))[0]
        

    
    def all_relatives_late_sbjct_start(self, fuzzy=10):
        """
        Do all relatives indicate a missing start, i.e. the majority
        of relatives have tuples (query_start < fuzzy, sbjct_start >
        fuzzy), indicating subject start is later than query start.

        TOADD: factor in evolutionary distance.
        """
        
        

    
    def report_first_hsp_starts(self):
        """
        A report method: print out how consistent the start subject
        sites are.
        """
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

        # TODO check: are these gauranteed in best first order?
        self.relatives[relative] = list()
        best_alignment = blast_record.alignments[0]
        for hsp in best_alignment.hsps:
            hsp_dict = {'align_length':best_alignment.length,
                        'align_title': best_alignment.title,
                        'e':hsp.expect,
                        'identities':hsp.identities,
                        'frame':hsp.frame,
                        'query_start':hsp.query_start,
                        'sbjct_start':hsp.sbjct_start,
                        'sbjct':hsp.sbjct}
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

    blast_files = parse_blastx_args(args.blastx)

    results = dict()

    for relative, blast_file in blast_files.items():
        
        sys.stderr.write("processing BLAST file '%s'...\t" % relative)

        for record in NCBIXML.parse(blast_file):
            query_id = record.query.strip().split()[0]

            # initialize a ContigSequence object for this query/contig
            if not results.get(query_id, False):
                results[query_id] = ContigSequence()

            # add the relative's alignment information
            results[query_id].add_relative_alignment(relative, record)
                    
        sys.stderr.write("done\n")
        
    cPickle.dump(results, file=args.output)



def predict_orf(args):
    """
    
    """
    # first, we index the FASTA reference
    ref_index = SeqIO.index(args.ref, "fasta")
    ref_ids = ref_index.keys()

    data = cPickle.load(args.input)
    for relative, hsps in data.items():
        hsps.get_hsp_frames()
        hsps.get_hsp_start_tuples()
    
    return {'ref':ref_index, 'data':data}

def put_seq_in_frame(seq, frame, alphabet=IUPAC.unambiguous_dna):
    """
    Take a sequence and transform it to into the correct frame.
    """
    # seq = Seq(seq_str, alphabet) # TODO: handle?
    if frame < 0:
        seq = seq.reverse_complement()
        frame = -1*frame
    if not frame in range(1, 4):
        pdb.set_trace()
        raise Exception, "improper frame: frame must be in [1, 3]"
    return seq[(frame-1):]

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
    parser_join.set_defaults(func=join_blastx_results)

    ## predict arguments
    parser_predict = subparsers.add_parser('predict', help="predict ORFs")
    parser_predict.add_argument('--input', type=argparse.FileType('rb'),
                                default="joined_blastx_dbs.pkl",
                                help="the joined results of all BLASTX files (Python pickle file; from sub-command joined)")
    parser_predict.add_argument('--ref', type=str, required=False,
                                help="the FASTA reference that corresponds to BLASTX queries.")
    parser_predict.add_argument('--full-length', action="store_true",
                                help="the FASTA reference that corresponds to BLASTX queries.")

    parser_predict.set_defaults(func=predict_orf)

    args = parser.parse_args()

    ## Argument checking
    results = args.func(args)

    
    
