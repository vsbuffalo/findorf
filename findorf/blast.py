"""
Functions for:

 - processing arguments related to BLAST databases.

 - joining seperate relative's BLAST results, making new Contig
   objects for each contig and taking BLASTX alignments and converting
   HSPs of the top alignment into a SeqRanges object.
"""
import sys
import cPickle

try:
    from Bio import SeqIO
    from Bio.Blast import NCBIXML
    from Bio.SeqRecord import SeqRecord
except ImportError:
    sys.exit("Cannot import BioPython modules; please install it.")

import contig

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

def make_blast_args(args_str):
    """
    From a series of arguments passed, form arg:param dictionary.
    """
    args = args_str.split()
    keys = args[::2]
    values = args[1::2]

    if len(keys) != len(values):
        raise ValueError("error parsing blastx arguments; no positional"
                         " arguments should be passed outside of input' and "
                         "'databases'")
    return dict(zip(keys, values))

def add_blastx_results(ref, blastx):
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
    for record in SeqIO.parse(ref, "fasta"):
        contigs[record.id] = contig.Contig(record)
    sys.stderr.write("done.\n")
    
    # For each relative alignment, add the HSPs via Contig's
    # add_relative method
    blast_files = parse_blastx_args(blastx)

    for relative, blast_file in blast_files.items():
        sys.stderr.write("[join] processing BLAST file '%s'...\t" % relative)

        for record in NCBIXML.parse(blast_file):
            query_id = record.query.strip().split()[0]

            # add the relative's alignment information, which actually
            # adds HSPs to relative HSPs attribute
            contigs[query_id].add_alignment(relative, record)

        sys.stderr.write("done.\n")

    return contigs
