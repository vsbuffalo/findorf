# findorf: ORF prediction of de novo transcriptome assemblies

`findorf` is an ORF prediction and transcriptome contig annotation
tool that uses the results of separate BLASTX results against close
relatives. It also annotates cases where a contig appears to be
missing a 5'-end, a predicted ORF is missing a 3'-end, internal stop
codons, and frameshifts.

## Requirements

`findorf` requires Python (version >= 2.7) and
[BioPython](http://biopython.org).

## Installation

Installation is easy: just clone or download this repository and enter
in the directory:

    python setup.py install

## `findorf` Join

`findorf` first joins the contig sequence FASTA file with the results
of each separate BLASTX against relatives using the `join`
subcommand. This is to ensure that if prediction is run with different
parameters, this step is not unnecessarily repeated.

    findorf join --ref contigs.fasta at:blast-a_thaliana_alt.xml bd:blast-b_distachyon_alt.xml \
      zm:blast-z_mays_alt.xml os:blast--o_sativa_alt.xml

Note that it is *highly* recommended organism abbreviation names are
provided (otherwise they'll be extracted from the basename). These are
then used throughout the predict stage as identifiers.

## `findorf` Predict

The `predict` subcommand predicts ORFs and annotates contigs and
ORFs. Annotation refers not to biological or functional annotation
(there are many great pieces of software for this), but rather
annotation about the state of the contig as what its relative hits
have to say about it. The following attributes are gathered for each
contig, before making an ORF prediction.

`predict` takes many options, for varying types of output.

    findorf predict --gtf orfs.gtf --protein proteins.fasta --fasta orfs.fasta \
      --dense dense.out -I -v

In this case, `findorf` would predict ORFs (`predict`), output
translated proteins (`--protein`), nucleotide ORFs (`--fasta`), GTF
(`--gtf`), and a dense output (`--dense`), be verbose about it (`-v`),
and then go interactive (`-I`) to allow Python-speaking users to look
at the data more closely.

Entering `findorf predict --help` will list all options.

### Majority Frame

`findorf` determines the **majority frame** across all relatives. The
majority is calculated by the number of identities per each frame:
more identities in a certain frame, more evidence that this is the
correct frame. In practice, there is very little disagreementacross
relatives about the frame in our test data.

### Majority and Any Frameshift

As with finding the majority frame, finding the **majority
frameshift** is based on the number of identities. If there are more
identities from HSPs disagreeing on the frame (per relative) than not,
we say there's a majority frameshift.

There's also **any frameshift**, which indicates if any relative says
there's a frameshift.

### Inconsistent Strand

`findorf` also checks that HSPs are on the same strand (allowing for
different frames due to frameshift). This could be due do a local
translocation, mis-assembly (conjoined contigs), or overlap.

### Anchor HSPs

For each relative, we calculate the **anchor HSPs**, which are the 5'-
and 3'-most HSPs on a contig (respecting strand). These are done for
each relative a contig has, subsetting by e-value and per-relative
percent identity thresholds (via the method `get_relatives()`). If a
relative's HSP falls outside of these thresholds, it is not
considered. This has two primary goals:

1. Let phylogenetic a priori beliefs about how close a relative should
guide the ORF prediction and annotation process.

2. Guard against potentially incorrect databases.

Note that if a contig has a single long HSP, this becomes the 5'- and
3'-most HSP.

### Missing 5'-end

`findof` annotates cases of a **missing 5'-end** based on whether our
5'-most anchor HSP corresponds to a subject protein that appears to be
cut. For example, suppose our contig's 5'-most HSP has a query start
position of 1 (the first nucleotide; assume this aligned to frame
1). If our subject protein starts 45 amino acids in, it's very likely
our contig is missing the 5'-end.

However, to increase sensitivity, `findorf` has fuzzy defaults: the
query start of the 5'-most HSP must be 16 nucleotides or less into the
contig sequence, and the subject start must be greater than 40 amino
acids. See the figure below for an illustration of these parameters.

          
             query start
              |   HSP
        |------------------------------------------| contig
              |||||||||||
      |.......|---------| subject
    subject
     start
        

## ORF Prediction

Finally, with these necessary requisite attributes, `findorf` can
proceed and try to predict the ORF. There are three separate ORF
prediction procedures, depending on the attributes of the contig
gathered from the relatives: (1) predict in the case of a frameshift,
(2) predict in teh case of a missing 5'-end, and (3) predict in the
case of a plain (vanilla) case. Note that it's possible that the data
reveal the case of a frameshift and missing 5'-end, which is handleded
by the frameshift prediction function.

All three cases lead to **candidate** ORFs. The chosen candidate is
the one that overlaps the 5'-most anchor HSP of the closest relative.

In some cases, a candidate HSP could be missing a start or stop codon;
in these cases, essentially `findorf` looks for one-sided overlap with
the 5'-most anchor HSP.

### Majority Frameshift Prediction

If there's a majority frameshift, `finforf` uses the
`predict_ORF_frameshift` procedure. First, `findorf` checks if there's
also a missing 5'-end. If so, it sets `in_reading_frame` to True, to
allow an ORF with a missing start codon to be added to the candidate
list. If not, `in_reading_frame` is False, meaning the procedure must
find a full length ORF in the 5'-end.

Because a frameshift could still lead to a functional protein, this
prediction procedure operators as a ribosome would: it reads in the
frame of the 5'-most anchor HSP of the closest relative.

### Mising 5'-end Prediction

If a majority of relatives indicate that the 5'-end is missing from a
contig, `findorf` sets `in_reading_frame` to True, assuming that this
contig is missing the start codon. This is the first ORF of the
candidate ORFs.

### Vanilla Prediction

If there's no frameshift or suspected missing 5'-end, `findorf` simply
finds all ORFs in a sequence and adds them to the candidate list.

### Candidates to Prediction, and ORF Annotation

With these candidates, `findorf` then predicts the ORF that overlaps
the 5'-most anchor HSP of the closest relative.

With an ORF prediction, `findorf` then annotates this particular ORF,
indicating whether it missing a start or stop codon. Also, `findorf`
then goes through each relative's 3'-most anchor HSPs and notes if
there's any cases where these lie outside this ORF. If so, this ORF is
annotated possibly containing an internal stop codon, as there 3'
regions of the contig that contain alignment HSPs. All relatives are
used for this step for maximum sensitivity.

## Annotation Inconsistencies 

In some cases, annotation `findorf` produces may appear
inconsistent. `findorf` does not try to push every possible case into
a category. Furthermore, some steps work with the closest relative
when it makes more sense (finding the best ORF out of a list), while
others, for maximum sensitivity, use all relatives.

This can lead to tricky cases where the ORF prediction or annotation
terms disagree. For example, missing 5'-end could be true, while
missing start codon could be false. Why could this occur?

Recall missing 5'-end is a *majority* rule; the majority of relatives
have to agree that this is the case before it is annotated as
such. However, in choosing the best ORF out of all possible ORFs in
the contig, we take the 5'-most ORF that overlaps the 5'-most anchor
HSP of the *closest relative* (not majority). In this case our closest
relative could have disagreed with the majority.

Why shouldn't `findorf` just look at the closest relative? `findorf`
was written for use in predicting and annotating ORFs in wheat, where
some of the closest relatives (barley and brachypodium) are less
well-studied species than say, arabidopsis. For this reason, some of
the closest relative database information may be
incorrect. Incorporating all relative information into finding frame
and anchor HSPs will guard against this, but still leave us with some
incosistent cases that are best manually curated.