"""
HMMER tabular output parser.

This is needed because HMMER (somewhat stupidly) outputs human
readable tabular output with variable number of spaces as
delimiters. This wouldn't be a problem, but the last column definitely
contains spaces. Furthermore, is that a guarantee that domain doesn't
have a space? If it does, each numeric column could be offset by one
(silently!).
"""
import pdb
import re

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
        colname = header[pos:pos+width]
        if i+1 == num_columns:
            end = None
        else:
            end = pos+width
        assert(len(colname) > 0)
        columns[colname.strip()] = (pos, end)
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
    
if __name__ == "__main__":
    import hmmer
    test_file = "All6frames.PfamA.E1e-3.domE1.domtbl"
    parser = hmmer.make_hmmer_parser(test_file)
    parsed = parser(test_file)
