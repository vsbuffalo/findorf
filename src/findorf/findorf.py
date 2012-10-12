#!/usr/bin/env python
__version__ = 1.02

info = """
findorf - ORF prediction and RNA contig annotation using comparative genomics
version: %s

Vince Buffalo (vsbuffaloAAAA@gmail.com, sans poly-A tail)
Dubcovsky Lab, Plant Sciences, UC Davis

Basic usage: run `findorf join` to join XML BLASTX results against
separate relatives' databases, then run `findorf predict` to output
ORF predictions.

When specifying XML BLASTX result databases of each relative, use the
format: short_id:filename-of-blastx-results.xml, i.e.:

    at:athaliana.xml bd:brachy.xml

This is necssary for later specifying relative-specific percent
identity constraints (as integers out of 100), i.e.:

    at:78,90 bd:80,95
    
""" % __version__

import sys
import pdb
import csv
import argparse
import os
import code
from string import Template
from collections import Counter, namedtuple, defaultdict
import csv
import cPickle
from operator import itemgetter, attrgetter
try:
    from Bio.Blast import NCBIXML
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError, e:
    sys.exit("Cannot import BioPython modules; please install it.")

import templates
from contig import Contig, GTF_FIELDS
from utilities import pretty_summary, join_blastx
import blast
import findall

# Which annotation keys to include in counting.
SUMMARY_KEYS = set(['majority_frameshift', 'orf',
                    'any_frameshift', 'missing_5prime',
                    'has_relatives', 'closest_relative_anchor_hsps_diff_frame',
                    'missing_start', # see note at rules.predict_ORF_vanilla
                    'missing_stop', 'no_hsps_coverages',
                    'full_length', 'contains_stop'])

def go_interactive(contigs, summary):
    try:
        import readline
        import utilities
    except ImportError:
        pass
    else:
        import rlcompleter
        readline.parse_and_bind("tab: complete")                    
        
    code.InteractiveConsole(locals=dict(contigs=contigs, summary=summary,
                                        u=utilities)).interact()

def predict_orf(args):
    """
    First, parse the relative percenty identity arguments.
    """
    # initiate summary counting
    counter = Counter()
    counter['total'] = 0

    sys.stderr.write("[predict] loading contig objects...")
    all_contigs = cPickle.load(open(args.input, 'rb'))
    sys.stderr.write("\tdone.\n")
    
    pi_range_args = parse_percent_identity_args(args)

    sys.stderr.write("[predict] predicting and annotating ORFs.\n")
    for query_id, contig in all_contigs.items():
        # Predict ORF and update contig annotation
        if args.verbose:
            if counter['total'] % 1000 == 0:
               sys.stderr.write("\t%d out of %d contigs processed\r"
                                % (counter['total'], len(all_contigs)))
               sys.stderr.flush()
        orf = contig.predict_orf(args.e_value, pi_range_args)
        counter['total'] += 1

    if args.protein is not None:
        sys.stderr.write("[predict] writing protein sequences...")
        proteins = [x.protein for x in all_contigs.values() if x.protein is not None]
        SeqIO.write(proteins, args.protein, "fasta")
        args.protein.close()
        sys.stderr.write("\tdone.\n")

    if args.fasta is not None:
        sys.stderr.write("[predict] writing nucleotide sequences...")
        seqs = [x.orf_seq for x in all_contigs.values() if x.orf_seq is not None]
        SeqIO.write(seqs, args.fasta, "fasta")
        args.fasta.close()
        sys.stderr.write("\tdone.\n")
    
    if args.gtf is not None:
        sys.stderr.write("[predict] writing GTF...")
        gtffmt = '\t'.join(["$" + g for g in GTF_FIELDS]) + '\n'

        for c in all_contigs.values():
            args.gtf.write(Template(gtffmt).substitute(c.gtf_dict()))
        args.gtf.write("\n")
        args.gtf.close()
        sys.stderr.write("\tdone.\n")

    if args.dense is not None:
        sys.stderr.write("[predict] writing dense...")
        for c in all_contigs.values():
            args.dense.write(repr(c))
        args.dense.close()
        sys.stderr.write("\tdone.\n")

    if args.frameshift is not None:
        sys.stderr.write("[predict] writing frameshifts...")
        seqs = [SeqRecord(seq=c.seq, id=c.id) for c in all_contigs.values() if c.get_annotation('majority_frameshift')]
        SeqIO.write(seqs, args.frameshift, "fasta")
        args.frameshift.close()
        sys.stderr.write("\tdone.\n")

    if args.stop is not None:
        sys.stderr.write("[predict] writing contigs with internal stop codons...")
        seqs = [SeqRecord(seq=c.seq, id=c.id) for c in all_contigs.values() if
                not c.get_annotation('majority_frameshift') and
                c.internal_stop]
        SeqIO.write(seqs, args.stop, "fasta")
        args.stop.close()
        sys.stderr.write("\tdone.\n")

    if args.no_relatives is not None:
        sys.stderr.write("[predict] writing contigs with no relatives...")
        seqs = [SeqRecord(seq=c.seq, id=c.id) for c in all_contigs.values() if
                not c.get_annotation('has_relatives')]
        SeqIO.write(seqs, args.no_relatives, "fasta")
        args.no_relatives.close()
        sys.stderr.write("\tdone.\n")

    if args.five_prime_utrs is not None:
        sys.stderr.write("[predict] writing 5'-UTRs...")
        seqs = [SeqRecord(seq=c.three_prime_utr(), id=c.id + " 5'-UTR sequence") for c in all_contigs.values()]
        seqs = filter(lambda x: x.seq is not None, seqs)
        SeqIO.write(seqs, args.five_prime_utrs, "fasta")
        args.five_prime_utrs.close()
        sys.stderr.write("\tdone.\n")

    if args.three_prime_utrs is not None:
        sys.stderr.write("[predict] writing 3'-UTRs...")
        seqs = [SeqRecord(seq=c.five_prime_utr(), id=c.id + " 3'-UTR sequence") for c in all_contigs.values()]
        seqs = filter(lambda x: x.seq is not None, seqs)
        SeqIO.write(seqs, args.three_prime_utrs, "fasta")
        args.three_prime_utrs.close()
        sys.stderr.write("\tdone.\n")

    sys.stdout.write(pretty_summary(all_contigs.values()))
        
    if args.interactive:
        go_interactive(all_contigs, None)
    # # for debugging with python -i
    # return all_contigs

def join_blastx_results(args):
    """
    Bridge args input to join_blastx_results so that the latter has a
    cleaner interface when used with module import.
    """
    join_blastx(args.ref, args.blastx, args.output)
     
def parse_percent_identity_args(args):
    """
    Parse the percent identity threshold arguments, making them None
    if they are not specified.
    """
    try:
        if args.i is None:
            return None
        tmp = [tuple(re.split(r':|,', i)) for i in args.i if i is not None]
        pi_ranges = dict([(k, (l, u)) for k, l, u in tmp])
    except Exception, e:
        msg = ("error parsing percent identity thresholds; "
               "must be in format rel:x,y where rel is the relative "
               "identifier from the join operation and 0 < x, y <= 100.\n")
        raise argparse.ArgumentTypeError(msg)
    return pi_range_args

def run_blast(args):
    """
    Run blastx on the input file. All arguments are extracted from
    sys.argv and passed directly to blast+'s blastx.
    """
    blastx_args = blast.make_blast_args(args.blast_args) # last entry should be input file
    if len(blastx_args) > 0:
        sys.stderr.write("[blast] passing args '%s' to blast.\n" % blastx_args.items())

    sys.stderr.write("[blast] making blast calls...")
    databases = blast.extract_databases(args.databases)
    sys.stderr.write("\tdone.\n")
    
    blast.blast_all_relatives(args.input, databases, **blastx_args)

def findall_orfs(args):
    findall.findall(args.contigs, min_length=args.min_length, translate=args.translate)

def main():
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
    parser_join.add_argument('--ref', type=str, required=True,
                             help=("the FASTA reference that corresponds "
                                   "to BLASTX queries."))

    parser_join.set_defaults(func=join_blastx_results)

    ## predict arguments
    parser_predict = subparsers.add_parser('predict', help="predict ORFs")
    
    parser_predict.add_argument('--input', type=str,
                                default="joined_blastx_dbs.pkl",
                                help=("the joined results of all BLASTX "
                                      "files (Python pickle file; from "
                                      "sub-command joined)"))
    parser_predict.add_argument('--gtf', type=argparse.FileType('w'),
                                default=None,
                                help="filename of the GTF of ORF predictions")
    parser_predict.add_argument('--dense', type=argparse.FileType('w'),
                                default=None,
                                help="filename of the dense output file")
    # parser_predict.add_argument('--full-length', action="store_true",
    #                             help="output file of full length ORFs")
    #                             "to BLASTX queries"))
    parser_predict.add_argument('-e', '--e-value', type=float,
                                default=10e-3,
                                help=("e-value threshold (relative hit "
                                      "only include if less than this)"))

    parser_predict.add_argument('-v', '--verbose', action="store_true",
                                default=False,
                                help="Output a period every 1,000 contigs predicted")
    parser_predict.add_argument('-I', '--interactive', action="store_true",
                                default=False,
                                help="drop into interactive mode after prediction")
    
    parser_predict.add_argument('-i', type=str, nargs="+",
                                default=None,
                                help=("A relative-specific percent identity "
                                      "range in teh format at:99,100, where"
                                      "at is the same identifier used in join"))
    parser_predict.add_argument('-f', '--fasta', type=argparse.FileType('w'), 
                                default=None,
                                help="FASTA file to output full length ORFs")
    parser_predict.add_argument('-p', '--protein', type=argparse.FileType('w'), 
                                default=None,
                                help="FASTA file to output translated proteins")
    parser_predict.add_argument('-F', '--frameshift', type=argparse.FileType('w'), 
                                default=None,
                                help="FASTA file to output frameshifted ORFs")
    parser_predict.add_argument('-s', '--stop', type=argparse.FileType('w'), 
                                default=None,
                                help="FASTA file to output internal stop codon ORFs (but not frameshift)")
    parser_predict.add_argument('-n', '--no-relatives', type=argparse.FileType('w'), 
                                default=None,
                                help="FASTA file to output contigs without relatives")
    parser_predict.add_argument('-3', '--three-prime-utrs', type=argparse.FileType('w'), 
                                default=None,
                                help="FASTA file to output the 3'-UTRs, for disjoining")
    parser_predict.add_argument('-5', '--five-prime-utrs', type=argparse.FileType('w'), 
                                default=None,
                                help="FASTA file to output the 5'-UTRs, for disjoining")

    parser_predict.set_defaults(func=predict_orf)

    # ## blast arguments
    # parser_blast = subparsers.add_parser('blast', help="run blastx")
    # parser_blast.add_argument('-p', '--processes', type=int,
    #                           default=1,
    #                           help="the number of processes to distribute blast calls across")
    # parser_blast.add_argument('-a', '--blast-args', type=str,
    #                           help="a quoted set of blastx arguments")
    # parser_blast.add_argument('input', type=str,
    #                           help="the FASTA input file")
    # parser_blast.add_argument('databases', type=str, nargs="+",
    #                           help="blast relative databases")
    # parser_blast.set_defaults(func=run_blast)

    ## findall arguments
    parser_findall = subparsers.add_parser('findall', help="find all ORFs by brute force")
    # parser_findall.add_argument('-k', '--kullback-leibler', action="store_true",
    #                             default=True,
    #                             help="Pick the best ORF based on K-L divergence to known ORFs")
    parser_findall.add_argument('-f', '--full-length', type=argparse.FileType('r'),
                                default=None,
                                help="FASTA file of full length ORFs")
    parser_findall.add_argument('-t', '--translate', help="translate sequences", action="store_true")
    parser_findall.add_argument('-m', '--min_length', help="minimum length to output", type=int, default=30)
    parser_findall.add_argument('contigs', type=str, 
                                help="FASTA file of contigs without relatives")
    parser_findall.set_defaults(func=findall_orfs)
    args = parser.parse_args()

    ## Run the appropriate step
    return args.func(args)

if __name__ == "__main__":
    out = main()
