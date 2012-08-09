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

A ContigSequence is an object that summarizes the contig sequence and
relatives' information (from blastx). The `rules` module is a set of
functions applied to these objects, and some rules will require that
specific attributes of the ContigSequence be propagated.
"""

import sys
import pdb
import cPickle
import csv
from collections import Counter, namedtuple, defaultdict
from string import Template
from operator import itemgetter, attrgetter
try:
    from Bio.Blast import NCBIXML
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError, e:
    sys.exit("Cannot import BioPython modules; please install it.")
import argparse
import os

import templates
import rules
from ContigSequence import ContigSequence

def mean(x):
    """
    The arithematic mean.
    """
    return float(sum(x))/len(x) if len(x) > 0 else float('nan')

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
        msg = "could not open BLASTX XML result '%s' - no such file.\n"
        sys.exit(msg % e.filename)
    return handles

def predict_orf(args):
    """
    TODO
    """
    all_contig_seqs = cPickle.load(args.input)

    for query_id, contig_seq in all_contig_seqs:
        if contig_seq.has_relatives:
            rules.predict_ORF_vanilla(contig_seq)
    
    return contig_seqs
        

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
    # Load all sequences into new ContigSequence objects
    sys.stderr.write("loading all contig sequences...\t")
    results = dict()
    for record in SeqIO.parse(args.ref, "fasta"):
        results[record.id] = ContigSequence(record.id, record.seq)
    sys.stderr.write("done\n")

    # For each relative alignment, add the HSPs via ContigSequence's
    # add_relative method
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=info)
    subparsers = parser.add_subparsers(help="sub-commands")

    ## join arguments
    parser_join = subparsers.add_parser('join',
                                        help=("merge contigs FASTA file with "
                                              "XML BLASTX results"))
    parser_join.add_argument('blastx', type=str, 
                             nargs="+", help="the output XML of BLASTx")
    parser_join.add_argument('--output', type=argparse.FileType('wb'),
                             default="joined_blastx_dbs.pkl",
                             help=("joined results of all the BLASTX "
                              "results files (Python pickle file)"))
    parser_join.add_argument('--ref', type=str, required=False,
                             help=("the FASTA reference that corresponds "
                                   "to BLASTX queries."))

    parser_join.set_defaults(func=join_blastx_results)

    ## predict arguments
    parser_predict = subparsers.add_parser('predict', help="predict ORFs")
    parser_predict.add_argument('--input', type=argparse.FileType('rb'),
                                default="joined_blastx_dbs.pkl",
                                help=("the joined results of all BLASTX "
                                      "files (Python pickle file; from "
                                      "sub-command joined)"))
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
                                help=("the FASTA reference that corresponds "
                                "to BLASTX queries."))
    parser_predict.add_argument('-e', '--e-value', type=float,
                                default=10e-3,
                                help=("e-value threshold (relative hit "
                                      "only include if less than this)"))
    parser_predict.set_defaults(func=predict_orf)

    args = parser.parse_args()

    ## Run the appropriate step
    results = args.func(args)
