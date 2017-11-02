Changelog
=========

0.1.2 (11/2/2017)
-----------------
 - Adds/correctly updates MC tag
 - Fixed erroneous mate info when mate is filtered out
 
   - Correctly sets mate start pos to 0
   - Removes MC tag if present
 
 - Sanitizes BAM tag of primer names
 - Extended supported pysam versions to include 0.9-0.12 


0.1.1 (2/9/2016)
----------------
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
---------------
 - Initial Release
