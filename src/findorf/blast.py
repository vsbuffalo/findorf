"""
blast.py contains functions to launch standalone blast runs and join
the output.

"""
import tempfile
import sys
from Bio.Blast.Applications import NcbiblastxCommandline
from os.path import join, basename, expanduser, exists
from multiprocessing import Pool
import subprocess

REQUIRED_BLASTX_PARAMS = dict(num_alignments=1,
                              num_descriptions=1, outfmt=5)

DEFAULT_BLASTX_PARAMS = dict(evalue=0.001)


def blast_all_relatives(input_file, databases, num_processes=4, outdir=None, **blastx_args):
    """
    The core blasting function; called internally and via command like
    subcommand 'blast' (which is first run through findorf.run_blast
    to parse args). This is 

    `input_file` is a string filename, `daabases` is a dict of
    relative:db-file.fasta key/values, and `blastx_args` are all
    blastx arguments.
    """
    sys.stderr.write("[blast] making blast calls...")
    blast_calls = make_all_relative_blastx_calls(input_file, databases, outdir=outdir,
                                                 **blastx_args)
    sys.stderr.write("\tdone.\n")

    # formats = dict(5='tab', 6='xml')

    if num_processes == 1:
        for call in blast_calls:
            sys.stderr.write("[blast] calling blast (no MP)...")
            call()
            sys.stderr.write("\tdone.\n")
    else:
        sys.stderr.write("[blast] calling blast across worker pool...\n")
        p = Pool(processes=num_processes)

        # force the command via the BioPython standalone blast objects
        blast_cmds = map(str, blast_calls)
        results = p.map(blast_cmd_caller, blast_cmds)
        if 1 in set(results):
            # one blast call ended in exit status 1; find out which
            # and error out.
            failed_cmds = [cmd for i, cmd in enumerate(blast_cmds) if results[i] == 1]
            raise ValueError("[blast] error: the following blastx commands had "
                             "non-zero exit status:\n" + 
                             '\n'.join(failed_cmds))
        sys.stderr.write("[blast] blast worker pool complete.\n")

    # sys.stderr.write("[blast] outputting a string to be passed to join process "
    #                  "via command line.\n")
    # sys.stdout.write(' '.join(["%s:%s" % (k, v) for k, v in databases.items()]))


def blast_cmd_caller(x):
    sys.stderr.write("[blast] running blast command '%s'" % str(x))
    return subprocess.call(x, shell=True)

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

def extract_databases(db_list):
    """
    Make a dictionary out of rel:database.fasta files.
    """
    # extract database information
    databases = dict()
    for db in db_list:
        try:
            name, database = db.split(":")
        except ValueError:
            raise ValueError("incorrectly formatted database argument: "
                             "must be of form 'relative:database_file.fasta'")
        if not exists(expanduser(database)):
            raise ValueError("database '%s' does not exist" % database)
        databases[name] = expanduser(database)

    return databases

def call_blast(blast_commands):
    pass


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
                                        out=join(outdir, "%s.%s" % (relative, format)),
                                        **params)
        except ValueError, error:
            msg = "a BLASTX parameter was specified that does not exist:\n"
            msg += "\t" + error.message
            raise ValueError(msg)
        
        blasters.append(cmd)

    return blasters

