======
Katana
======

Command-line tool to soft-clip reads from amplicon-based sequence based on
specified primer locations.

.. image:: https://travis-ci.org/umich-brcf-bioinf/Katana.svg?branch=develop
    :target: https://travis-ci.org/umich-brcf-bioinf/Katana
    :alt: Build Status

.. image:: https://coveralls.io/repos/github/umich-brcf-bioinf/Katana/badge.svg?branch=develop
    :target: https://coveralls.io/github/umich-brcf-bioinf/Katana?branch=develop
    :alt: Coverage Status

The official repository is at:

https://github.com/umich-brcf-bioinf/Katana

--------
Overview
--------

In amplicon-based target panel sequencing, regions-of-interest are amplified by
specific pairs of primers; consequently the regions-of-interest typically
always start and end with these primer sequences, sequences which match the
reference sequence exactly and do not reflect the actual sample sequence. In
some panel designs, the amplicons may be tiled such that an amplicon of one
region of interest may overlap the primer region of different amplicon. In this
arrangement, the overlapping regions should enable detection of variants that
fall within that primer region. However, the presence of the primer sequences
will typically overwhelm the signature of true, low-frequency variants.


Katana matches each read to its corresponding primer pair based on start
position of the read. Katana then soft-clips the primer region from the edge of
the read sequence, rescuing the signal of true variants measured by overlapping
amplicons. The output is conceptually similar to hard-clipping the primers from
the original FASTQ reads based on sequence identity but with the advantage that
retaining the primers during alignment improves alignment quality.
::
                      [ primer REGION-OF-INTEREST primer ]
  input read sequence:  TGCATG AGTCTGATCTAGGTAGTT GACGTC
  output read sequence: tgcatg AGTCTGATCTAGGTAGTT gacgtc (lowercase = soft-clipped)


Tags are added to each output read to help explain how it was modified:
 - X0 : associated target_id
 - X1 : original cigar string
 - X2 : original reference start
 - X3 : original reference_end (not modified, but FYI)
 - X4 : why read would be excluded (appears only if --preserve_all_alignments)


Katana assumes that:
 - input bam is indexed
 - primers come in sense-antisense pairs
 - primer pairs are on the same chromosome
 - primer chromsomes match the bam regions
 - primer file is tab separated; the header line includes the following fields:
   - Customer TargetID
   - Chr
   - Sense Start
   - Antisense Start
   - Sense Sequence
   - Antisense Sequence
 - primer file sense and antisense start are specified in 1-based coordinates


-----------
Quick Start
-----------

1. **Install Katana (see INSTALL.rst):**
::
  $ pip install katana

2. **Get the examples directory:**
::
  $ git clone https://github.com/umich-brcf-bioinf/Katana

3. **Run Katana:**
::
  $ katana Katana/examples/primers.txt Katana/examples/chr10.pten.bam clipped.bam

This will read chr10.pten.bam and produce clipped.bam which contains reads
adjusted to soft-clip (exclude) their respective primer regions. Unmapped reads
or reads which do not match a known primer are excluded.


-----------
Katana help
-----------

::

  $ katana --help
  
  usage: katana primer_manifest input_bam output_bam
  
  Match each alignment in input BAM to primer, softclipping the primer region.
  
  positional arguments:
    primer_manifest       path to primer manifest (tab-separated text)
    input_bam             path to input BAM
    output_bam            path to output BAM
  
  
  optional arguments:
    -h, --help            show this help message and exit
    -V, --version         show program's version number and exit
    --preserve_all_alignments
                          Preserve all incoming alignments (even if they are 
                          unmapped, cannot be matched with primers, result in 
                          invalid CIGARs, etc.)

====

Email bfx-katana@umich.edu for support and questions.

UM BRCF Bioinformatics Core
