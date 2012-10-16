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

from blast import join_blastx_results
import argparse

def _join_blastx_results(args):
    """
    Private function that bridges args input to join_blastx_results so
    that the latter has a cleaner interface when used with module
    import.
    """
    join_blastx_results(args.ref, args.blastx, args.output)

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
    parser_join.add_argument('--ref', type=str, required=True,
                             help=("the FASTA reference that corresponds "
                                   "to BLASTX queries."))

    parser_join.set_defaults(func=_join_blastx_results)

    # Parse arguments and run the appropriate step
    args = parser.parse_args()
    return args.func(args)

if __name__ == "__main__":
    main()
