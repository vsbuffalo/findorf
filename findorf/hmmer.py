"""
HMMER tabular output parser.

This is needed because HMMER (somewhat stupidly) outputs human
readable tabular output with variable number of spaces as
delimiters. This wouldn't be a problem, but the last column definitely
contains spaces. Furthermore, is that a guarantee that domain doesn't
have a space? If it does, each numeric column could be offset by one
(silently!).

Even more ridiculous is that if one tries to make a parser based on
fixed widths, this still isn't sufficient to parse HMMER output. Why?
Because the column headers are split across two rows. This is why we
have to specify columns.
"""
import pdb
import re

HMMER_COLS = [
"target_name",
"accession",
"tlen",
"query_name",
"accession",
"qlen",
"evalue_fullseq",
"score_fullseq",
"bias_fullseq",
"num_domain",
"total_domain",
"cEvalue_domain",
"iEvalue_domain",
"score_domain",
"bias_domain",
"hmm_from",
"hmm_to",
"ali_from",
"ali_to",
"env_from",
"env_to",
"acc",
"description"]

def make_hmmer_parser(hmmer_file):
    """
    Return a function that parses HMMER files, based on looking at the
    third row's column widths.
    """
    col_widths_values = list()
    with open(hmmer_file) as f:
        for i, line in enumerate(f):
            if i == 1:
                header = line.strip()
            if i == 2:
                assert line.startswith("#-")
                # get max column width based on these delimiters
                splits = re.findall(r"[ \#]+-+", line)
                col_widths_values = map(len, splits)
                break
    pos = 0
    widths = list()
    columns = dict()
    num_columns = len(col_widths_values)
    for i, width in enumerate(col_widths_values):
        colname = HMMER_COLS[i]
        if i+1 == num_columns:
            end = None
        else:
            end = pos+width
        columns[colname] = (pos, end)
        pos += width
        
    def parser(file):
        lines = list()
        with open(file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                row = dict()
                for field, pos in columns.iteritems():
                    start, end = pos
                    row[field] = line[start:end].strip()
                lines.append(row)
        return lines

    return parser

