# findorf

## Design

 - separate BLASTX runs against relatives, effect on K-A
   statistics. Why? because we also want a close relative with a short
   hit (= larger e-value) to be considered alongside a longer hit from
   a more distant relative. Sure this may not happen, but what about
   truncated proteins in databases?
 
 - identity-based scoring of frames

## Join

BLASTX and PFAM domain results for each relative *separately*
essentially gives us a ragged array of ranges (HSPs or PFAM domains)
over a contig. These need to be joined by contig ID, which is what
`findorf join` does. All relevant information provided in this step is
stored in a `Contig` object, which

