"""
utilities.py contains utilities that handle creating codons, puting
sequences in frame, etc.

"""
from collections import Counter
import sys
import contig
import cPickle

from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from RangedFeatures import ORF

## Biological constants TODO get from BioPython
STOP_CODONS = set(["TAG", "TGA", "TAA"])
START_CODONS = set(["ATG"])

def get_codons(seq, frame):
    """
    Return a list of (codon, position_in_orf, pos_in_forward_query)
    tuples, on the foward strand.
    """
    # for user friendliness with interactive shell, if seq is a
    # string, convert it to Bio.Seq.Seq
    if type(seq) is str:
        seq = Seq(seq)
        
    if frame < 0:
        seq = seq.reverse_complement()
        
    frame = abs(frame)
    tmp = [(str(seq[pos:pos+3]), pos-(frame-1), pos) for
           pos in range(frame-1, len(seq), 3)]
    
    # remove the last string if not a full codon.
    return [(c, po, pfq) for c, po, pfq in tmp if len(c) == 3]


def get_all_orfs(seq, frame, query_length,
                 in_reading_frame=False, add_partial=False):
    """
    Generic ORF finder; it returns a list of all ORFs as they are
    found, given codons (a list if tuples in the form (codon,
    position)) from `get_codons`.
    """
    
    all_orfs = list()
    orf_start_pos = None
    query_start_pos = None

    codons = get_codons(seq, frame)
    
    # Note that query_pos is forward strand.
    for codon, orf_pos, query_pos in codons:
        codon = codon.upper()
        if codon in START_CODONS and not in_reading_frame:
            in_reading_frame = True
            orf_start_pos = orf_pos
            query_start_pos = query_pos
            continue
        if in_reading_frame and codon in STOP_CODONS:
            # a full reading frame, unless we haven't hit any stop
            # yet, then we're in the possible ORF from the start of
            # the query to the end.
            all_orfs.append(ORF(query_start_pos, query_pos, query_length,
                                frame, no_start=query_start_pos is None,
                                no_stop=False))
            
            in_reading_frame = False
            orf_start_pos = None
            query_start_pos = None
            
    # add any partial ORFs, and mark as having no stop codon
    if in_reading_frame:
        all_orfs.append(ORF(query_start_pos, query_pos, query_length,
                            frame, no_start=query_start_pos is None,
                            no_stop=True))
    if add_partial:
        for codon, orf_pos, query_pos in codons:
            if codon in STOP_CODONS:
                all_orfs.insert(0, ORF(None, query_pos, query_length,
                                       frame, no_start=True,
                                       no_stop=True))
                break

    return all_orfs

def pretty_summary(contigs):
    summary = summarize_contigs(contigs)
    
    inorder_terms = ['has_orf',                      
                     'full_length', 'missing_start', 'missing_stop', 'internal_stop',
                     'majority_frameshift', 'missing_5prime',
                     'inconsistent_strand', 'has_relatives']

    out = ""
    for term in inorder_terms:
        out += term + ":\n"
        for key, count in summary[term].items():
            out += "  %s: %d\n" % (key, count)

    return out

def summarize_contigs(contigs):
    """
    Summarize annotation.
    """

    terms = ['majority_frameshift', 'missing_5prime', 'hsp_orf_overlap',
             'inconsistent_strand', 'has_orf', 'has_relatives',
             'full_length', 'contig_length',
             'missing_start', 'missing_stop', 'internal_stop',
             'num_relatives', 'num_orf_candidates', 'closest_relative']

    annotation_summary = dict([(t, Counter()) for t in terms])

    if type(contigs) is dict:
        contigs = contigs.values()
    
    for contig in contigs:
        all_anno = contig.get_annotation()
        for key, value in all_anno.items():
            annotation_summary[key][value] += 1

    return annotation_summary


def put_seq_in_frame(seq, frame):
    """
    Take a sequence (of Bio.Seq class) and transform it to into the
    correct frame.
    """
    if seq.__class__.__name__ != "Seq":
        seq = Seq(seq)
    if frame < 0:
        seq = seq.reverse_complement()
        frame = -1*frame

    if not frame in range(1, 4):
        raise Exception, "improper frame: frame must be in [1, 3]"
    return seq[(frame-1):]


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


def join_blastx(ref, blastx, output):
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
            contigs[query_id].add_relative_alignment(relative, record)

        sys.stderr.write("done.\n")

    # dump the joined contigs
    cPickle.dump(contigs, file=output)
