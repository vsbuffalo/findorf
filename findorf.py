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

    def __init__(self, record):
        """
        Initialize a ContigSequence via a BioPython SeqRecord.
        """
        self.record = record
        self.relatives = dict()

        ## annotation attributes
        self.annotation = dict(full_length_orf=None, missing5prime=None,
                               missing3prime=None, likely_psuedogene=None)
        self.all_hsps_same_frame = None
        
    def get_hsp_frames(self):
        """
        Per each relative, check if all HSPs are in the same frame.
        Return a dictionary of relatives and whether they are in the
        same frame.
        """
        
        

    def add_relative_alignment(self, relative, blast_record):
        """
        Given a relative and a BioPython BLAST alignment object,
        extract and store the relevant parts.
        """
        relative_exists = self.relatives.get(relative, False)

        if not relative_exists:
            self.relatives[relative] = dict()
        else:
            raise Exception, "relative '%s' already exists for this ContigSequence" % relative

        for alignment in blast_record:
            self.relatives[relative][alignment.title] = list()
            for hsp in alignment.hsps:
                hsp_dict = {'align_length':alignment.length,
                            'e':hsp.expect,
                            'identities':hsp.identities,
                            'frame':hsp.frame,
                            'query_start':hsp.query_start,
                            'sbjct_start':hsp.sbjct_start,
                            'sbjct':hsp.sbjct}
                
                self.relatives[relative][alignment.title].append(hsp_dict)

                


                
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
            if not results.get(query_id, False):
                results[query_id] = dict()
                
            results[query_id][relative] = dict()
            for align in record.alignments:
                results[query_id][relative][align.title] = list()
                for hsp in align.hsps:
                    results[query_id][relative][align.title].append({'align_length':align.length,
                                                                'e':hsp.expect,
                                                                'identities':hsp.identities,
                                                                'frame':hsp.frame,
                                                                'query_start':hsp.query_start,
                                                                'sbjct_start':hsp.sbjct_start,
                                                                'sbjct':hsp.sbjct
                                                                })

        sys.stderr.write("done\n")
        
    cPickle.dump(results, file=args.output)



def predict_orf(args):
    """
    
    """
    return cPickle.load(args.input)

def put_seq_in_frame(seq, frame, alphabet=IUPAC.unambiguous_dna):
    """
    Take a sequence and transform it to into the correct frame.
    """
    # seq = Seq(seq_str, alphabet) # TODO: handle?
    if frame < 0:
        seq.reverse_complement()
        frame = -1*frame
    if not frame in range(1, 4):
        pdb.set_trace()
        raise Exception, "improper frame: frame must be in [1, 3]"
    return seq[(frame-1):]

def annotate_orf(seq, table=1):
    """
    Annotate an ORF-candidate: given a sequence, add annotation attributes such
    as after tranlation:
      - full length
      - missing 5'
      - missing 3'

    etc.
    """
    pseq = seq.tranlate(table)
    

def hsps_different_frame(joined_results):
    """
    Iterate through all results and check that HSPs are in the same
    frame.
    """

    num_same_frame = 0
    num_diff_frame = 0
    full_length_orf = 0
    flo_list = list()
    start = "ATG"
    stop_codons = ["TAG", "TGA", "TAA"]
    for query_id, relative_alignments in joined_results.items():
        for relative, alignments in relative_alignments.items():
            for alignment, hsps in alignments.items():
                frames = [h['frame'][0] for h in hsps]
                if len(set(frames)) == 1:
                    num_same_frame += 1

                    # check if full length ORF
                    seq_in_frame = put_seq_in_frame(ref_index[query_id], frames[0]) # TODO check frame tuple

                    if start in seq_in_frame.seq and any([seq_in_frame.seq.find(c) > -1 for c in stop_codons]):
                        full_length_orf += 1
                        flo_list.append(seq_in_frame.seq)
                else:
                    num_diff_frame += 1
                
    return {'diferent_frame': num_diff_frame, 'same_frame':num_same_frame, 'full_orf':full_length_orf, "flos": flo_list}


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

    # # first, we index the FASTA reference
    ref_index = SeqIO.index(args.ref, "fasta")
    ref_ids = ref_index.keys()

    ## Argument checking
    results = args.func(args)

    a = hsps_different_frame(results)
    
    
