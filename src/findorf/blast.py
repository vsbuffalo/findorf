"""
blast.py contains functions to launch standalone blast runs and join
the output.

"""

from Bio.Blast.Applications import NcbiblastxCommandline
from os.path import join

REQUIRED_BLASTX_PARAMS = dict(num_alignments=1,
                              num_descriptions=1, outfmt=5)

DEFAULT_BLASTX_PARAMS = dict(evalue=0.001)

def make_all_relative_blastx_calls(seqfile, relatives, outdir="blasts",
                         processes=2, **blast_params):
    """
    Given a sequences file, and a dictionary of relatives; blastx all
    sequences against each, in parallel.
    """

    if any([True for param in blast_params.keys() if param
            in REQUIRED_BLASTX_PARAMS.keys()]):
        raise ValueError("BLASTX parameters cannot change REQUIRED_BLASTX_PARAMS")

    params = dict(REQUIRED_BLASTX_PARAMS.items() + DEFAULT_BLASTX_PARAMS.items())
    params = dict(params.items() + blast_params.items())
    
    blasters = list()
    for relative, database in relatives.items():
        try:
            cmd = NcbiblastxCommandline(query=seqfile, db=database, 
                                        out=join(outdir, "%s.xml" % relative),
                                        **params)
        except ValueError, error:
            msg = "a BLASTX parameter was specified that does not exist:\n"
            msg += "\t" + error.message
            raise ValueError(msg)
        
        blasters.append(cmd)

    return blasters

