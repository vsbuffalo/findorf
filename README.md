# findorf: ORF prediction of de novo transcriptome assemblies

`findorf` is an ORF prediction and transcriptome contig annotation
tool designed to be non-model organism-friendly.

## How is `findorf` different?

There are many approaches to ORF annotation (see
[ORFPredictor](http://proteomics.ysu.edu/tools/docs/OrfPredictor_faq.html),
[Dragon TIS](http://cbrc.kaust.edu.sa/dts.), and
[MetWAMer](http://www.biomedcentral.com/1471-2105/9/381)). Below are
some key design differences of `findorf`:

 - Designed to work on transcriptome assemblies rather than on genome
   gene models.

 - Unlike most Translation Initiation Site (TIS) prediction software,
   `findorf` does not train on nucleotide signal matrices. This
   approach is common in computationally spliced gene models coming
   from genomic sequences after gene prediction.

 - `findorf` uses BLASTX hits against seperated relative databases and
   PFAM domain information to infer ORF start position.

 - `findorf` uses relative information from BLASTX to annotate
   frameshifts, premature stop codons, etc.

 - `findorf` is designed with non-model plant species in mind.

## Requirements

`findorf` requires:
 - [Python](http://python.org) (version >= 2.7)
 - [BioPython](http://biopython.org)
 - [BioRanges](https://github.com/vsbuffalo/BioRanges) - a small
   package for handling ranges one sequences

## Installation

Installation is easy: just clone or download this repository and enter
in the directory:

    python setup.py install

## `findorf`'s General Approach

### Using Seperate BLASTX Databases

First, `findorf` relies upon the idea of BLASTXing each contig against
a database of relative proteins seperately, rather than a combined
protein database. The motivation for this is that databases of
computationally predicted proteins (say from EST or cDNA sequences)
may include incorrectly predicted proteins. In non-model organisms,
closest relative may also be non-model and its protein data less
likely to be validated by multiple source, undergo multiple prediction
algorithms, or be externally validates with proteomic data. Thus, ORF
prediction that looks at just top hits is not robust against this.

`findorf` tries to prevent this by seperate BLASTXs against all
relatives. Key decisions such as annotating a case as having a
majority frameshift or internal stop codon are done by looking at what
many relatives' data say.

After ORF prediction, the output GTF file contains annotation on the
relative that determined the start codon position. This can be used as
independent validation that 5' sites are primarily being chosen by
closer relatives (since it's likely more distant relatives 5' regions
would have diverged more).

### Choosing a Start Site

After infering frame and noting any frameshifts or strand
inconsistencies, `findorf` then chooses its start site. There are two
methods `findorf` can use: **5'-most methionine**, and extension from
**5'-most HSP**. We recommend using the latter.

Both techniques start but enumerating every possible ORF including
overlapping cases. This is illustrated below:

              |--HSP1--|

       S---------6------E
            S----4------E           S----------5----------E
    s------------1------E  S-------------------2----------E   S-3-e
    |--M----M-----------*--M--------M---------------------*---M---|

    S: ORF start with start codon
    s: ORF start with missing start codon
    E: ORF end with stop codon
    e: ORF end with missing stop codon

Note that the contig above has six candidates, two with open ended
cases (indicated by lower case 's' and 'e'). Suppose HSP1 is our
5'-most HSP. Under the **5'-most HSP** rule, ORF 4 is chosen because
it overlaps the 5' most HSP (HSP1) with the least 5'
extension. Essentially, this method chooses the start site based on
whatever evidence we have.

If the **5'-most methionine** approach were used, we face an
ambiguity: should we chose the open ended case (1) or the closed case
(6)? If there is a non-open ended ORF (that is, there is a methionine
5' of the HSP), this is chosen. If there were two non-open ORF cases,
we would choose the one with the earliest methionine. But in this
approach, the handling of missing 5'-ends is less clear.

Between the two approaches there's a tradeoff between the cost of
possibly choosing an internal methionine and choosing the outermost
methionine and possibly predicting part of the UTR is coding
sequence. [Sparks and
Brendel](http://www.biomedcentral.com/1471-2105/9/381) (2008) show
that the 1st ATG approch works very well 94% specificity and
sensitivity (assuming complete 5' regions), but one may decide the
cost of predicting a more 5' start site is higher.

### Robust Against Mis-Assembly and Chimeric Contigs

A key feature of `findorf` is that it can also output contigs
sequences with the predicted ORF hardmasked (with "X"s). This allows
one to BLASTX these masked sequences and run `findorf` a second
iteration. If futher ORFs are found in these non-masked regions, it's
a candidate chimeric contig (as there shouldn't be homology with
protein sequences in the UTRs*). In these cases, we can use ORF
prediction in assembled transcriptome sequences to also judge the
incidence of misassembled contigs.

*: Note that there are cases when we this could occur, i.e. if a gene
 is interrupted with a nonsense mutation but still remains functional
 or an internal TIS start site was chosen accidentally.

### Including PFAM Domains

The ends of proteins can vary considerably; since `findorf` start site
choice is based on the 5'-most HSP with a sequence, there's also the
option to use PFAM domains to detect possible domains in novel protein
fusions. This is especially important in 5' regions were a lack of
N-terminus homology can confound ORF start prediction. 

HMMER's `hmmscan` can be used to annotate PFAM domains using HMM
approaches. We recommend running the command as such:

    hmmscan -E 0.001 --domE 1 --tblout <tblout> --domtblout <domtblout> -o <outfile> --noali <database> <infile>

The `<tblout>` file would then be passed to `findorf join` with
`--domain-hits`. **You also must tell `findorf` to use PFAM results in
the prediction process** with `-u` or `--use-pfam`. PFAM domains will
only affect ORF start site choice; they are not used for annotating
frameshifts or internal stop codons (due to the fact a PFAM domain 3'
of a stop codon could be due to a chimeric contig). Prediction cases
that are extended based on PFAM domains will have the key
`pfam_extended_5prime` set to `True` in the GTF file.

## Running `findorf`: Join

`findorf` first joins the contig sequence FASTA file with the results
of each separate BLASTX against relatives using the `join`
subcommand. This is to ensure that if prediction is run with different
parameters this step is not unnecessarily repeated.

    findorf join --ref contigs.fasta at:blast-a_thaliana_alt.xml bd:blast-b_distachyon_alt.xml \
      zm:blast-z_mays_alt.xml os:blast-o_sativa_alt.xml

Note that it is *highly* recommended organism abbreviation names are
provided (otherwise they'll be extracted from the basename). These are
then used throughout the predict stage as identifiers.

## Running `findorf`: Predict

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

## Prediction Methods: Pre-Prediction Attributes

`findorf` first gathers some information about a contig based on HSPs
from the relatives' seperate BLASTX results and PFAM domains. These
attributes are described below and are the requisite information for
ORF prediction.

### Relatives

If there are no BLASTX hits to a contig, ORF predict does not predict
an ORF. One could take these cases and use an *ab initio* procedure
based on coding potential, PFAM domains, etc.

### Inconsistent Strand

`findorf` also checks that HSPs are on the same strand (allowing for
different frames due to frameshift). This could be due do a local
translocation, mis-assembly (conjoined contigs), or overlap. In does
not predict in these degenerate cases and they are labelled in the GTF
with the `inconsistent_strand` key.

### Majority Frame

`findorf` determines the **majority frame** across all relatives. The
majority is calculated by the number of identities per each frame:
more identities in a certain frame, more evidence that this is the
correct frame. In practice, there is very little disagreementacross
relatives about the frame in our test data.

### Query and Subject Start Attributes

`findorf` annotates (in the GTF) the 5'-most HSP query start and query
end. In earlier versions, we experimented with using cutoffs of these
values to infer whether the 5'-end of a contig was missing. Currently
these threshold-based approaches are not being used, but one may wish
to use query and subject start positions as diagnostics. One approach
is to plot the query and subject start as points (one for every
contig) and color these by the type of ORF (missing 5', missing 3',
partial, full-length, etc).

             query start
              |   HSP
        |------------------------------------------| contig
              |||||||||||
      |.......|---------| subject
    subject
     start

## Prediction Methods: ORF Choice

Finally, with these necessary requisite attributes, `findorf` can
proceed and try to predict the ORF.

### 5'- and 3'-most Anchor HSPs

The first step is to gather the HSPs (from any relative) that are the
5'-most and 3'-most. These are subject to the e-value based filtering
(command line option `-evalue`, default 1e-5). The 5'-most HSP is what
determines ORF start position (more on this below).

### Majority Frameshift Prediction

If there's a majority frameshift, `findorf` does not use the majority
frame, but rather the frame of the 5'-most HSP. A contig with a
frameshift mutation could still be functional. In this case, the
prediction procedure operators as a ribosome would: it reads in the
frame of the 5'-most anchor HSP of the closest relative. These cases
are annotated in the GTF as: `majority_frameshift True`.

### PFAM 

If `--use-pfam` is set, `findorf` will then take all PFAM domains and
remove those with a frame that differs from the majority frame. **Note:
the specified e-value threshold does not apply to PFAM domain hits**
(since e-values are calculated differently). If you wish to filter
PFAM domains by e-value, you must do so via `awk` or the `hmmscan`
tool.

With these domain hits in the same frame as the majority frame,
`findorf` then checks if any is more 5' than the 5'-most anchor
HSP. If so, this is used as the 5'-most range during ORF prediction.

### ORF Enumeration and Overlap Finding

As illustrated in the "Choosing a Start Site" section above, `findorf`
enumerates all ORF possibilities, including the open-ended cases and
overlapping cases. This list of ORF candidates is then subset by those
overlapping the 5'-most HSP (or PFAM hit). Those that overlap are then
chosen according to the **5'-most HSP** or **5'-most methionine** rule
(discussed more above). Cases in which there is no start or stop codon
are annotated.

### Internal Stop Codon Annotation

To be robust against the possibility of chimeric contigs, an internal
(or premature) stop codon is annotated as the case in which there are
more relatives that have an HSP that:

1. Overlaps the ORF.  

2. Have an end position greater than the ORF stop codon + 60bp (known
as `buffer_bp` in the code):

                       ORF end
        ------------------|   buffer_bp
        ---------------------------|--------| HSP end

The overlap requirement protects against the case in which a chimeric
contig has an HSP (to the other chimeric mRNA) past the stop site.


