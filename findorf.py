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
try:
    from Bio.Blast import NCBIXML
    from Bio.Alphabet import generic_dna, DNAAlphabet
    from Bio import SeqIO
except ImportError, e:
    sys.exit("Cannot import BioPython modules; please install it.")

import argparse
import os

parser = argparse.ArgumentParser(description=info)
parser.add_argument('blastx', type=str, 
                    nargs="+", help="the output XML of BLASTx")
parser.add_argument('--ref', type=str, required=False,
                    help="the FASTA reference that corresponds to BLAST queries.")

def parse_blastx_args(bargs):
    """
    Take a list of args and return a dictionary of names and files. If
    any arg is of the key:value sort, the name will be the
    user-defined key; otherwise they will be extracted from the file
    name.
    """
    blastx_files = dict()
    for arg in bargs:
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
    return blastx_files



    
if __name__ == "__main__":

    ## Argument checking
    args = parser.parse_args()

    print parse_blastx_args(args.blastx)


    # # first, we index the FASTA reference
    # ref_index = SeqIO.index(args.ref, "fasta")
    # ref_ids = ref_index.keys()

    
