#! /usr/bin/env python   
""" Basic parser for the Rhim Thunderbolts manifest file. """

from __future__ import print_function, absolute_import, division

import sys, os, re
import pybedtools
import pysam

from ampliconsoftclipper import __version__

class PrimerPairRecord(object):
    """Represents a primer pair derived from a MiSeq manifest. """

    def __init__(self, id, target_id, primer_set, chrom, original_start, original_end, converted_end, converted_start, genome_build, sense_start,
                 antisense_start, sense_sequence, antisense_sequence, sense_sequence_tailed_illumina, antisense_sequence_tailed_illumina):
        self.id = id
        self.target_id = target_id
        self.primer_set = primer_set
        self.chrom = 'chr'+chrom
        self.orginal_start = int(original_start)
        self.original_end = int(original_end)
        self.converted_start = int(converted_start)
        self.converted_end = int(converted_end)
        self.genome_build = genome_build
        self.sense_start = int(sense_start)
        self.antisense_start = int(antisense_start)
        self.sense_sequence = sense_sequence
        self.antisense_sequence = antisense_sequence
        self.sense_sequence_tailed_illumina = sense_sequence_tailed_illumina
        self.antisense_sequence_tailed_illumina = antisense_sequence_tailed_illumina
        
    
    @property
    def to_fasta(self):
        fa_header1 = ">{0}_{1}_{2}_F".format(self.id, self.target_id, self.sense_sequence)
        outline1 = "{0}\n{1}\n".format(fa_header1, self.sense_sequence.upper())
        fa_header2 = ">{0}_{1}_{2}_R".format(self.id, self.target_id, self.sense_sequence)
        outline1 = "{0}\n{1}\n".format(fa_header1, self.sense_sequence.upper())
        




def parse_thunderbolts_manifest(filename):
    """ Parses a Rhim manifest file to pull out primer data. Returns a dict. """
    primer_l = []
    datafile = file(filename)
    header = datafile.readline()
    for line in datafile.readlines():
        bits = line.strip().split("\t")
        pr = PrimerPairRecord(*bits[0:15])
        primer_l.append(pr)
    return primer_l


def get_region(primer_pair_record):
    return (primer_pair_record.chrom,
            primer_pair_record.sense_start,
            primer_pair_record.antisense_start)


def get_bedtools_region(primer_pair_record):
    regions = [(primer_pair_record.chrom,
                primer_pair_record.sense_start,
                primer_pair_record.antisense_start)]
    return pybedtools.BedTool(regions)


def get_reads_with_coords(bam_fn, chrom, start, end):
    bamfile = pysam.AlignmentFile(bam_fn, 'rb')
    read_l = []
    read_total = 0
    forw_total = 0
    rev_total = 0
    fmatch_total = 0
    rmatch_total = 0
    print("coords: {0}\t{1}\n".format(start, end))
    for read in bamfile.fetch(chrom, start, end):
        read_total += 1
        if (read.flag & 16 == 0):
            forw_total += 1
            if (read.reference_start == start):
                fmatch_total += 1
                read_l.append(read)
        else:
            rev_total += 1
            if (read.reference_end == end):
                rmatch_total += 1
                read_l.append(read)
    print("{0}\t{1}\t{2}\t{3}\t{4}".format(read_total, forw_total, rev_total, fmatch_total, rmatch_total))
    return read_l



def main():
    if (len(sys.argv) != 4):
        print("usage: {0} [thunderbolts manifest] [bam] [outfolder]".format(os.path.basename(sys.argv[0])))
        sys.exit()

    bam_fn = sys.argv[2]
#    bam_bedtool = pybedtools.BedTool(bam_fn)
    outfolder = sys.argv[3]

    thunderbolt_manifest = sys.argv[1]
    primers = parse_thunderbolts_manifest(thunderbolt_manifest)
    (chrom, start, end) = get_region(primers[185])

#    bam_bedtool.intersect(test_bed, f=0.90, r=True).saveas(outfolder+'/'+'testout.pten.bam')
    
    selected_reads = get_reads_with_coords(bam_fn, chrom, start-1, end)
    print(len(selected_reads))

    print("done.")

if __name__ == '__main__':
    main()

