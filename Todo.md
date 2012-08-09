# Todo

## Things to Do

1. Is there a stop codon between the 5'-most HSP and the end of the
sequence?

2. Use the first frame of a frameshifted contig to predict ORF.

3. If there's a stop codon in in from the contig start to the query start, we ignore it.

## Cases to Check


 - k36_contig_9886 starts with a start codon, but has a missing 5'-end.

 - k21_contig_36350: has a stop codon in first codon and missing 5'-end

 - k36_contig_9886

 - k51_contig_10673

 - k61_contig_20415 - problem detecting frameshift, not in BLAST results.

 - k26_contig_24653

 - k26_contig_22146 - frameshift, but orf start/end

