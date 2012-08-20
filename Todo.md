# Todo

## Things to Do

- Study codon bias in real vs fake ORFs

- Check that orf_overlaps_5prime is working correctly.

- If there's a stop codon in in from the contig start to the query start, we ignore it

- Add qs_start, etc missing 5prime params as args

## Cases to Check
 - k36_contig_11881 - premature stop
 - k41_contig_27125 - tons of HSPs

 - k26_contig_4655 misassembled, finding ORF in 5'

 - lots of different frames: k31_contig_16146

 - k31_contig_37537 missing start and frameshift

 - k61_contig_59200 no frameshift, one in ORFpredictor

 - k36_contig_9886 starts with a start codon, but has a missing 5'-end.

 - k21_contig_36350: has a stop codon in first codon and missing 5'-end

 - k36_contig_9886

 - k51_contig_10673

 - k61_contig_20415 - problem detecting frameshift, not in BLAST results.

 - k26_contig_24653

 - k26_contig_22146 - frameshift, but orf start/end


## Neat contigs

 - k41_contig_61398 r-gene