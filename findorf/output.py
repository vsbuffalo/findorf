"""
Functions for outputting different files. These return closures bound
to the open file handle, and these closures have one argument: the
contig object to write the appropriate data from.
"""
import sys
from string import Template

try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
except ImportError:
    sys.exit("Cannot import BioPython modules; please install it.")

GTF_FIELDS = ("seqname", "source", "feature", "start",
              "end", "score", "strand", "frame", "group")        

def protein_writer(contigs, file):
    sys.stderr.write("[predict] writing protein sequences...")
    proteins = filter(lambda x: x.protein is not None, contigs.values())
    SeqIO.write(proteins, file, "fasta")
    file.close()
    sys.stderr.write("\tdone.\n")


def orf_writer(contigs, file):
    sys.stderr.write("[predict] writing nucleotide sequences...")
    seqs = [x.orf_seq for x in contigs.values() if x.orf_seq is not None]
    SeqIO.write(seqs, file, "fasta")
    file.close()
    sys.stderr.write("\tdone.\n")

def gtf_writer(contigs, file):
    sys.stderr.write("[predict] writing GTF...")
    gtffmt = '\t'.join(["$" + g for g in GTF_FIELDS]) + '\n'
    for c in contigs.values():
        file.write(Template(gtffmt).substitute(c.gtf_dict()))
    file.write("\n")
    file.close()
    sys.stderr.write("\tdone.\n")

def frameshift_writer(contigs, file):
    sys.stderr.write("[predict] writing frameshifts...")
    seqs = [SeqRecord(seq=c.seq, id=c.id) for c in contigs.values()
            if c.annotation['majority_frameshift']]
    SeqIO.write(seqs, file, "fasta")
    file.close()
    sys.stderr.write("\tdone.\n")

def stop_writer(contigs, file):
    sys.stderr.write("[predict] writing contigs with internal stop codons...")
    seqs = [SeqRecord(seq=c.seq, id=c.id) for c in contigs.values() if
            not c.annotation['majority_frameshift'] and
            c.annotation['internal_stop']]
    SeqIO.write(seqs, file, "fasta")
    file.close()
    sys.stderr.write("\tdone.\n")

def no_relatives_writer(contigs, file):
    sys.stderr.write("[predict] writing contigs with no relatives...")
    seqs = [SeqRecord(seq=c.seq, id=c.id) for c in contigs.values() if
            c.annotation['num_relatives'] == 0]
    SeqIO.write(seqs, file, "fasta")
    file.close()
    sys.stderr.write("\tdone.\n")

def five_prime_utrs_writer(contigs, file):
    sys.stderr.write("[predict] writing 5'-UTRs...")
    seqs = [SeqRecord(seq=c.three_prime_utr(), id=c.id + " 5'-UTR sequence") for c in contigs.values()]
    seqs = filter(lambda x: x.seq is not None, seqs)
    SeqIO.write(seqs, file, "fasta")
    file.close()
    sys.stderr.write("\tdone.\n")

def three_prime_utrs_writer(contigs, file):
    sys.stderr.write("[predict] writing 3'-UTRs...")
    seqs = [SeqRecord(seq=c.five_prime_utr(), id=c.id + " 3'-UTR sequence") for c in contigs.values()]
    seqs = filter(lambda x: x.seq is not None, seqs)
    SeqIO.write(seqs, file, "fasta")
    file.close()
    sys.stderr.write("\tdone.\n")

WRITERS = {'protein':protein_writer, 'orf':orf_writer,
           'gtf':gtf_writer, 'frameshift':frameshift_writer,
           'stop':stop_writer, 'no_relatives':no_relatives_writer,
           'five_prime_utrs':five_prime_utrs_writer,
           'three_prime_utrs':three_prime_utrs_writer}
