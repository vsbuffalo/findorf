# hmmerfix.py - fix terrible output format (fixed width) of HMMER

# We make a simple assumption about the data: any non-delimiter spaces
# are in the last column. Under this assumption, we build a regular
# expression programmatically with strict types. Strict typing and the
# number of grouped elements ensure some safety.

import re
import sys
from collections import OrderedDict
from operator import itemgetter

# aliased parsers: here we need a match a different occurence then
# string, but we want the same type
def to_end(x):
    return str(x)

COMMENT_CHAR = "#"
NEW_DELIM = "\t"
DOMTBLOUT_FIELDS = (("target_name", str),
                    ("target_accession", str),
                    ("tlen", int),
                    ("query_name", str),
                    ("query_accession", str),
                    ("qlen", int),
                    ("seq_evalue", float),
                    ("seq_score", float),
                    ("seq_bias", float),
                    ("domain_num", int),
                    ("total_domains", int),
                    ("domain_cevalue", float),
                    ("domain_ievalue", float),
                    ("domain_score", float),
                    ("domain_bias", float),
                    ("hmm_from", int),
                    ("hmm_to", int),
                    ("ali_from", int), 
                    ("ali_to", int),
                    ("env_from", int),
                    ("env_to", int),
                    ("acc", float),
                    ("description", to_end))

matchers = dict(float=r"(?P<%s>[\d\.e+-]+)",
                int=r"(?P<%s>[\d]+)",
                str=r"(?P<%s>[^\s]+)",
                to_end=r"(?P<%s>.*$)")

def build_matcher(fields, matchers, delim="\s+"):
    """
    Build a function that takes a tuple of (field, function) tuples
    and build a regular expression from it. `function` must be a
    function to convert the type, and have a __name__ attribute in
    matchers.
    """
    re_parts = list()
    processors = OrderedDict(fields)
    for name, func in fields:
        matcher = matchers[func.__name__] % name
        re_parts.append(matcher)
    re_string = delim.join(re_parts)
    sys.stderr.write("Using regular expression:\n%s\n\n" % re_string)
    FIELD_PARSER = re.compile(re_string)
    def parser(line):
        parsed_line = FIELD_PARSER.match(line)
        if parsed_line is None:
            raise ValueError("error parsing line; check fields and matchers.\n%s" % line)
        parsed_line_dict = parsed_line.groupdict()
        line_dict = OrderedDict((k, processors[k](parsed_line_dict[k])) for k in map(itemgetter(0), fields))
        if len(line_dict) != len(fields):
            raise ValueError("differing number of values extracted than specified, check fields.\n%s" % line)
        return line_dict
    return parser

if __name__ == "__main__":
    domtblout_file = open(sys.argv[1])
    domtblout_parser = build_matcher(DOMTBLOUT_FIELDS, matchers)

    header = True
    for line in domtblout_file:
        if line.startswith("#"):
            continue
        entry_dict = domtblout_parser(line)
        if header:
            sys.stdout.write(NEW_DELIM.join(str(s) for s in entry_dict.keys()) + "\n")
            header = False
        sys.stdout.write(NEW_DELIM.join(str(s) for s in entry_dict.values()) + "\n")
    
    
