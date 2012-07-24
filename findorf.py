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
try:
    from Bio.Blast import NCBIXML
    from Bio.Alphabet import IUPAC, generic_dna, DNAAlphabet
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError, e:
    sys.exit("Cannot import BioPython modules; please install it.")

import argparse
import os

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

        ## annotation attributes
        self.annotation = dict(full_length_orf=None, missing5prime=None,
                               missing3prime=None, likely_psuedogene=None)
        
    def get_hsp_frames(self):
        """
        Per each relative, check if all HSPs of the first alignment
        are in the same frame.  Return a dictionary of relatives and
        whether they are in the same frame.

        We care about a few things here:

        1. All the HSPs of the top hit all in the same frame, for each
        relative?

        2. Do all relatives agree on the frame?
        
        TODO we also want to compare the top hit *across* relatives
        """
        frames = dict()
        all_frames = set() # for across-relative frame comparison
        for relative, hsps in self.relatives.iteritems():
            f = [h['frame'][0] for h in hsps] # TODO check for other tuple elements?
            frames[relative] = list(set(f))
            all_frames.update(f)
        self.frames = frames
        self.all_relatives_same_frame = True if len(all_frames) == 1 else False

        # note whether there's a frameshift in the relatives TODO:
        # note whether it's one HSP or many, and also how good they
        # are.
        self.frame_shifts_by_relative = dict([(r, len(f) > 1) for r, v in frames.items()])
        self.num_hsps_by_relative = dict([(r, len(v)) for k, v in self.relatives.items()])

        #frame_shifts_in_all_relatives = [k for k, x in a.items() if all(x.frame_shifts_by_relative.values())]

    def get_relative_first_hsp_start(self):
        """
        Annotate whether cases where:

        1. The query starts with 1.

        2. The subject starts at a position > 1.

        3. The subject has a methionine.
        """
        # TODO earliest HSP: in all frames?
        self.relative_subject_starts = dict()
        for relative, hsps in self.relatives.iteritems():
            # TODO should we make this fuzzy?
            hsp_with_query_at_1 = [h for h in hsps if h['query_start'] == 1]
            if len(hsp_with_query_at_1):
                # get the hsp with query start at 1 with lowest e
                tmp = sorted(hsp_with_query_at_1, key=lambda x: x['e'])[0]
                self.relative_subject_starts[relative] = tmp['sbjct_start']


        
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

    return {'ref':ref_index, 'blast':cPickle.load(args.input)}

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
    parser_predict.set_defaults(func=predict_orf)

    args = parser.parse_args()

    ## Argument checking
    results = args.func(args)

    
    
