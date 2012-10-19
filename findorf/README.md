# findorf

## Design

 - separate BLASTX runs against relatives, effect on K-A statistics.
 
 - identity-based scoring of frames

## Join

BLASTX and PFAM domain results for each relative *separately*
essentially gives us a ragged array of ranges (HSPs or PFAM domains)
over a contig. These need to be joined by contig ID, which is what
`findorf join` does. All relevant information provided in this step is
stored in a `Contig` object, which

