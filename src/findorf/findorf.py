#!/usr/bin/env python
__version__ = 0.9

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
# import pdb
import csv
import argparse
import os

from string import Template
from collections import Counter, namedtuple, defaultdict
import csv
import cPickle
from operator import itemgetter, attrgetter
try:
    from Bio.Blast import NCBIXML
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError, e:
    sys.exit("Cannot import BioPython modules; please install it.")

import templates
from Contig import Contig, GTF_FIELDS


# Which annotation keys to include in counting.
SUMMARY_KEYS = set(['majority_frameshift', 'orf',
                    'any_frameshift', 'missing_5prime',
                    'has_relatives', 'closest_relative_anchor_hsps_diff_frame',
                    'missing_start', # see note at rules.predict_ORF_vanilla
                    'missing_stop', 'no_hsps_coverages',
                    'full_length', 'contains_stop'])

def predict_orf(args):
    """
    First, parse the relative percenty identity arguments.
    """
    # initiate summary counting
    counter = Counter()
    counter['total'] = 0
    
    all_contig_seqs = cPickle.load(open(args.input, 'rb'))
    pi_range_args = parse_percent_identity_args(args)
    total = 0
    for query_id, contig_seq in all_contig_seqs.items():
        if args.verbose: # FIXME
            if counter['total'] % 1000 == 0:
                sys.stderr.write('.')

        # Predict ORF and update contig annotation
        contig_seq.generic_predict_ORF(args.e_value, pi_range_args)
        contig_seq.annotate_contig()

        # Increment counters for this contig's annotations.
        for attribute, value in contig_seq.annotation.iteritems():
            if attribute in SUMMARY_KEYS:
                counter[attribute] += contig_seq.get_annotation(attribute) is True

        counter['total'] += 1

    ## Output various formats, we can make this a single loop later.
    if args.dense is not None:
        with args.dense as f:
            for cs in all_contig_seqs.values():
                if cs.has_relatives:
                    f.write("----------------%s" % str(cs))
    if args.gtf is not None:
        with args.gtf as f:
            dw = csv.DictWriter(f, GTF_FIELDS, delimiter="\t")
            for cs in all_contig_seqs.values():
                dw.writerow(cs.gtf_dict())
    if args.fasta is not None:
        with args.fasta as f:
            for cs in all_contig_seqs.values():
                if None not in (cs.orf):
                    f.write(">%s\n%s\n" % (cs.query_id, cs.orf.get_orf(cs.seq)))
    if args.protein is not None:
        with args.protein as f:
            for cs in all_contig_seqs.values():
                if None not in (cs.orf):
                    f.write(">%s\n%s\n" % (cs.query_id, cs.orf.get_orf(cs.seq).translate()))

        
    sys.stderr.write(Template(templates.out).substitute(counter))

    if args.interactive:
        contigs = all_contig_seqs
        summary = counter
        import code
        import rules
        try:
            import readline
        except ImportError:
            pass
        else:
            import rlcompleter
            readline.parse_and_bind("tab: complete")                    
        code.InteractiveConsole(locals=dict(contigs=contigs, summary=summary, rules=rules)).interact()


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
    # Load all sequences into new Contig objects
    contigs = dict()
    sys.stderr.write("[join] loading all contig sequences...\t")
    for record in SeqIO.parse(args.ref, "fasta"):
        contigs[record.id] = Contig(record)
    sys.stderr.write("done.\n")
    
    # For each relative alignment, add the HSPs via Contig's
    # add_relative method
    blast_files = parse_blastx_args(args.blastx)

    for relative, blast_file in blast_files.items():
        sys.stderr.write("[join] processing BLAST file '%s'...\t" % relative)

        for record in NCBIXML.parse(blast_file):
            query_id = record.query.strip().split()[0]

            # add the relative's alignment information, which actually
            # adds HSPs to relative HSPs attribute
            contigs[query_id].add_relative_alignment(relative, record)

        sys.stderr.write("done.\n")

    # dump the joined contigs
    cPickle.dump(contigs, file=args.output)

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
    parser_join.add_argument('--ref', type=str, required=False,
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
    parser_predict.add_argument('--full-length', action="store_true",
                                help=("the FASTA reference that corresponds "
                                "to BLASTX queries"))
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
    parser_predict.add_argument('-o', '--fasta', type=argparse.FileType('w'), 
                                default=None,
                                help="FASTA file to output full length ORFs")
    parser_predict.add_argument('-p', '--protein', type=argparse.FileType('w'), 
                                default=None,
                                help="FASTA file to output translated proteins")

    parser_predict.set_defaults(func=predict_orf)

    args = parser.parse_args()

    ## Run the appropriate step
    results = args.func(args)

if __name__ == "__main__":
    main()
