Changelog
=========

0.1.1 (XX/XX/XXXX)
-----------------
- Corrected problems in BAM output:
   - Corrected next reference in paired reads
   - Excludes reads where CIGAR is entirely clipped
   - Unpairs reads which had no mate in input
- Added support for pip installs
- Added BAM tags to excluded reads (useful when --preserving_all_reads)
- Adjusted to improve performance (about 6x faster)
- Added functional tests
- Added support for travis CI

0.1 (1/28/2016)
--------------
- Initial Release
