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
import cPickle
try:
    from Bio.Blast import NCBIXML
    from Bio.Alphabet import generic_dna, DNAAlphabet
    from Bio import SeqIO
except ImportError, e:
    sys.exit("Cannot import BioPython modules; please install it.")

import argparse
import os

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
            blastx_files[os.path.splitext(os.path.basename(arg))] = arg
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
    ID and each value being another dictionary in which each key is
    the BLAST result file and the values are the BLAST results.
    """

    blast_files = parse_blastx_args(args.blastx)

    results = dict()

    for key, blast_file in blast_files.items():
        
        sys.stderr.write("processing BLAST file '%s'...\t" % key)

        for record in NCBIXML.parse(blast_file):
            if not results.get(record.query, False):
                results[record.query] = dict()
                
            results[record.query][key] = dict()
            for align in record.alignments:
                results[record.query][key][align.title] = list()
                for hsp in align.hsps:
                    results[record.query][key][align.title].append(hsp)

        sys.stderr.write("done\n")
        
    return results

def predict_orf():
    """
    
    """
    pass


    
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
    parser_predict.add_argument('--input', help="the joined results of all BLASTX files (Python pickle file; from sub-command joined)")
    parser_predict.add_argument('--ref', type=str, required=False,
                                help="the FASTA reference that corresponds to BLASTX queries.")
    parser_predict.set_defaults(func=predict_orf)

    ## Argument checking
    args = parser.parse_args()
    args.func(args)
    
    print args
    # a = join_blastx_results(parse_blastx_args(args.blastx))
    # cPickle.dump(a, file=open("joined_blastx_dbs.pkl", "wb"))

    # # first, we index the FASTA reference
    # ref_index = SeqIO.index(args.ref, "fasta")
    # ref_ids = ref_index.keys()

