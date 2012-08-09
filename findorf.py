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
    from Bio.Alphabet import IUPAC, generic_dna, DNAAlphabet
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError, e:
    sys.exit("Cannot import BioPython modules; please install it.")
import argparse
import os

import templates

## Biological constants
STOP_CODONS = set(("TAG", "TGA", "TAA"))
START_CODONS = set(("ATG"))
GTF_FIELDS = ("seqname", "source", "feature", "start",
              "end", "score", "strand", "frame", "group")

## Named Tuples for lightweight object storage
OrfSet = namedtuple('OrfSet', ['start', 'stop', 'length', 'rank'])
HSP = namedtuple('HSP', ['e', 'identities', 'length',
                         'percent_identity', 'title',
                         'query_start', 'query_end',
                         'sbjct_start', 'sbjct_end',
                         'frame'])
def mean(x):
    """
    The arithematic mean.
    """
    return float(sum(x))/len(x) if len(x) > 0 else float('nan')

class ContigSequence():
    """
    A ContigSequence is an assembled contig, that may be coding or
    non-coding.
    """

    def __init__(self, query_id, sequence):
        """
        Initialize a ContigSequence with a contig ID and sequence. The
        contig ID must correspond to the same one used in the blastx
        results.
        """
        # core data attributes
        self.query_id = query_id
        self.seq = sequence
        self.len = len(sequence)

        # information added by blastx results
        self.all_relatives = dict()
        
    def __repr__(self):
        """
        A representation of the object for dense output and
        interactive debugging.
        """

        info = dict(id=self.query_id, length=self.query_length,
                    num_relatives=self.num_relatives,
                    consensus_frame=self.consensus_frame,
                    majority_frame=self.majority_frame,
                    any_frameshift=self.any_frameshift,
                    majority_frameshift=self.majority_frameshift,
                    missing_start=self.missing_start,
                    missing_stop=self.missing_stop,
                    missing_5prime=self.missing_5prime,
                    full_length_orf=self.full_length_orf,
                    orf_start = self.orf_start,
                    orf_stop=self.orf_stop, seq=self.orf)
        
        out = Template(templates.contig_seq_repr).substitute(info)

        for relative, start_tuple in self.start_tuples.iteritems():
            query_start, sbjct_start, strand = start_tuple
            rel_info = (relative, sbjct_start, query_start, {1:"+", -1:"-"}[strand])
            out += ("%s\n    subject start: %s\n    query start/end"
            " (if strand forward/reverse): %s\n    strand: %s\n" % rel_info)

        # in later versions, we could use a templating engine...
        if self.has_relatives:
            out += "\n# Relative Identities in Frames\n"
            for relative, count_frames in self.frames_identities.iteritems():
                if len(count_frames):
                    out += "%s\n" % relative
                for frame, identities in count_frames.iteritems():
                    out += "  frame: %s\n  identities:  %s\n\n" % (frame, identities)
                    
        return out

    def get_relatives(self, e_value=None, identity=None):
        """
        Return relatives that pass thresholding filters.
        
        The `add_relative` method adds relatives' HSPs to a dictionary
        attribute, `all_relatives`. However, in most cases, we want to
        use a subset of these relatives that satisfy requirements
        based on phylogenetic needs, i.e. requiring a relative HSP
        have a percent identity consistent with evolutionary distance.

        If `e_value` or `identity` are None, they are not used for
        filtering `all_relatives`.
        """
        if e_value is None and identity is None:
            return self.all_relatives

        filtered_relatives = dict()
        for relative, hsps in self.all_relatives.items():
            filtered_relatives

        return passed_thresh
    

    def gff_dict(self):
        """
        Return a dictionary of some key attribute's values,
        corresponding to a GFF file's columns.

        Note that GFFs are 1-indexed, so we add one to positions.
        """
        out = dict()
        out["seqname"] = self.query_id
        out["source"] = "findorf"
        out["feature"] = "predicted_orf"
        out["start"] = self.orf_start + 1 if self.orf_start is not None else "."
        out["end"] = self.orf_stop + 1 if self.orf_stop is not None else "."
        out["score"] = "."

        if self.majority_frameshift is not None:
            out["strand"] = self.majority_frame/abs(self.majority_frame)
        else:
            out["strand"] = "."

        if self.majority_frame is not None:
            # GFF uses frames in [0, 2]
            out["frame"] = abs(self.majority_frame) - 1
        else:
             out["frame"] = "."
        out["group"] = "."
        return out

    def gtf_dict(self):
        """
        Return a dictionary corresponding to the columns of a GTF
        file.
        """

        # a GTF's file's "group" column contains a merged set of
        # attributes, which in ContigSequence's case are those below
        attributes = dict(full_length_orf=self.full_length_orf,
                          majority_frameshift=self.majority_frameshift,
                          any_frameshift=self.any_frameshift,
                          missing_5prime=self.missing_5prime,
                          number_relatives=len(self.relatives))

        group = "; ".join(["%s %s" % (k, v) for k, v in attributes.iteritems()])
        out = self.gff_dict()
        out["group"] = group
        return out

    def add_relative_alignment(self, relative, blast_record):
        """
        Given a relative and a BioPython BLAST alignment objects,
        extract and store the relevant parts of the _best_ alignment
        only.
        """
        relative_exists = self.all_relatives.get(relative, False)

        if not relative_exists:
            self.all_relatives[relative] = dict()
        else:
            msg = "relative '%s' already exists for this ContigSequence"
            raise Exception, msg % relative

        if len(blast_record.alignments) == 0:
            # no alignments, so we dont have any info to add for this
            # relative.
            return 

        self.all_relatives = defaultdict(list)

        # TODO check: are these guaranteed in best first order?
        best_alignment = blast_record.alignments[0]
        for hsp in best_alignment.hsps:
            percent_identity = hsp.identities/float(hsp.align_length)

            hsp = HSP(e=hsp.expect,
                      identities=hsp.identities,
                      length=hsp.align_length,
                      percent_identity=percent_identity,
                      title=best_alignment.title,
                      query_start=hsp.query_start,
                      query_end=hsp.query_end,
                      sbjct_start=hsp.sbjct_start,
                      sbjct_end=hsp.sbjct_end,
                      frame=hsp.frame)

            self.all_relatives[relative].append(hsp)


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
    pass

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
    
