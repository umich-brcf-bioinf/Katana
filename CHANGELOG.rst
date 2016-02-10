Changelog
=========

0.1.1 (XX/XX/XXXX)
-----------------
- Fixed problems in BAM output:
   - Corrected next reference in paired reads
   - Excludes reads where CIGAR is entirely clipped
   - Unpairs reads which had no mate in input
- Added BAM tags to excluded reads (useful when --preserving_all_reads)
- Adjusted to improve performance (about 6x faster)
- Added support for pip install
- Added functional tests
- Added support for travis CI
- Added support for Python3
- Added support for pysam 0.8.3

0.1 (1/28/2016)
--------------
- Initial Release
