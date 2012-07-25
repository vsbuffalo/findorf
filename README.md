# findorf: ORF prediction and annotation via blastx results of close relatives

## Computational Approach

The BLASTx XML format returns `iteration` XML nodes that correspond to
each query. These *should* be in the same order as the contig
reference FASTA file. With this, we could reduce our memory footprint
by opening up file handles for all BLASTx files and the reference
FASTA file and looping through them together.

Currently each XML result file is being looped through to propagate
the dictionary. This is a tad safer, a tad easier, and a bit slower.

## TODO

 - right now, we're looking at the earliest HSP by query_start
   position, but this ignores negative strand cases. FIX