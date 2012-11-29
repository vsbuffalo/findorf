__version__ = 1.02

INFO = """
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
import argparse
import cPickle

from hmmer import add_pfam_domain_hits
from blast import add_blastx_results
from orfprediction import predictall


def _join_relative_results(args):
    """
    Private function that bridges args input to join_blastx_results so
    that the latter has a cleaner interface when used with module
    import. Also, add any PFAM/HMMER input.
    """
    contigs = add_blastx_results(args.ref, args.blastx)
    if args.domain_hits is not None:
        add_pfam_domain_hits(contigs, args.domain_hits)
    # dump the joined contigs
    cPickle.dump(contigs, file=args.output)
    return contigs

def _predict_all_orfs(args):
    """
    Private function that bridges args input to
    orfprediction.predictall()
    """
    sys.stderr.write("[predict] loading contig objects...")
    contig_objects = cPickle.load(open(args.input, 'rb'))
    sys.stderr.write("\tdone.\n")

    possible_output = {'protein':args.protein, 'orf':args.orf,
                       'gtf':args.gtf, 'frameshift':args.frameshift,
                       'stop':args.stop, 'no_relatives':args.no_relatives,
                       'masked':args.masked}
    output_keys = filter(lambda k: possible_output[k] is not None, possible_output)
    to_output = dict([(k, possible_output[k]) for k in output_keys])
    
    method = '5prime-most' if args.most_5prime else '5prime-hsp'
    predictall(contig_objects, args.evalue, method, args.use_pfam,
               to_output, args.verbose)
    return contig_objects

def main():
    """
    main() is the entry point to findorf's functionality.
    """
    
    parser = argparse.ArgumentParser(description=INFO)
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
    parser_join.add_argument('--domain-hits', type=argparse.FileType('r'),
                             default=None,
                             help="PFAM domain hits via HMMER")
    parser_join.add_argument('--ref', type=str, required=True,
                             help=("the FASTA reference that corresponds "
                                   "to BLASTX queries."))

    parser_join.set_defaults(func=_join_relative_results)

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
    parser_predict.add_argument('-e', '--evalue', type=float,
                                default=1e-5,
                                help=("e-value threshold (relative hit "
                                      "only include if less than this)"))
    parser_predict.add_argument('-v', '--verbose', action="store_true",
                                default=False,
                                help="Output a period every 1,000 contigs predicted")
    parser_predict.add_argument('-o', '--orf', type=argparse.FileType('w'), 
                                default=None,
                                help="FASTA file to output ORFs")
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
    parser_predict.add_argument('-M', '--masked', type=argparse.FileType('w'), 
                                default=None,
                                help="Contigs with masked ORF region, for further iterations")
    parser_predict.add_argument('-u', '--use-pfam', action="store_true", default=False,
                                help="Use PFAM domains if they are available")
    parser_predict.add_argument('-m', '--most-5prime', default=False, action="store_true",
                                help="always use most 5' start codon")
    parser_predict.set_defaults(func=_predict_all_orfs)


    # Parse arguments and run the appropriate step
    args = parser.parse_args()
    args.func(args)
