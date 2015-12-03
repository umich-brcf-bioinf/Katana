#! /usr/bin/env python   
""" Basic parser for the Rhim Thunderbolts manifest file. """

from __future__ import print_function, absolute_import, division

import sys, os, re
import pybedtools
import pysam
from itertools import groupby

from ampliconsoftclipper import __version__


class Cigar(object):
    read_consuming_ops = ("M", "I", "S", "=", "X")
    ref_consuming_ops = ("M", "D", "N", "=", "X")

    def __init__(self, cigar_string):
        self.cigar = cigar_string

    def items(self):
        if self.cigar == "*":
            yield (0, None)
            raise StopIteration
        cig_iter = groupby(self.cigar, lambda c: c.isdigit())
        for g, n in cig_iter:
            yield int("".join(n)), "".join(next(cig_iter)[1])

    def __str__(self):
        return self.cigar

    def __repr__(self):
        return "Cigar('%s')" % self

    def __len__(self):
        """
        sum of MIS=X ops shall equal the sequence length.
        """
        return sum(l for l, op,in self.items() \
                               if op in Cigar.read_consuming_ops)

    def reference_length(self):
        return sum(l for l, op in self.items() \
                               if op in Cigar.read_consuming_ops)

    def mask_left(self, n_seq_bases, mask="S"):
        """
        Return a new cigar with cigar string where the first `n_seq_bases` are
        soft-masked unless they are already hard-masked.
        """
        cigs = list(self.items())
        new_cigs = []

        c, cum_len  = self.cigar, 0
        for i, (l, op) in enumerate(cigs):
            if op in Cigar.read_consuming_ops:
                cum_len += l
            if op == "H":
                cum_len += l
                new_cigs.append(cigs[i])
            elif cum_len < n_seq_bases:
                new_cigs.append(cigs[i])
            else:
                # the current cigar element is split by the masking.
                right_extra = cum_len - n_seq_bases
                new_cigs.append((l - right_extra, 'S'))
                if right_extra != 0:
                    new_cigs.append((right_extra, cigs[i][1]))
            if cum_len >= n_seq_bases: break
        else:
            pass

        new_cigs[:i] = [(l, op if op in "HS" else "S") for l, op in
                new_cigs[:i]]
        new_cigs.extend(cigs[i + 1:])
        return Cigar(Cigar.string_from_elements(new_cigs))


    @classmethod
    def string_from_elements(self, elements):
        return "".join("%i%s" % (l, op) for l, op in elements if l !=0)


    def mask_right(self, n_seq_bases, mask="S"):
        """
        Return a new cigar with cigar string where the last `n_seq_bases` are
        soft-masked unless they are already hard-masked.
        """
        return Cigar(Cigar(self._reverse_cigar()).mask_left(n_seq_bases, mask)._reverse_cigar())


    def _reverse_cigar(self):
        return Cigar.string_from_elements(list(self.items())[::-1])



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
    #primer_l = []
    primer_pos_d = {}
    datafile = file(filename)
    header = datafile.readline()
    for line in datafile.readlines():
        bits = line.strip().split("\t")
        ppr = PrimerPairRecord(*bits[0:15])
        primer_pos_d[(ppr.chrom, ppr.sense_start-1, '+')] = ppr
        primer_pos_d[(ppr.chrom, ppr.antisense_start, '-')] = ppr
        #primer_l.append(ppr)
    return primer_pos_d


def get_region(primer_pair_record):
    return (primer_pair_record.chrom,
            primer_pair_record.sense_start,
            primer_pair_record.antisense_start)


def get_bedtools_region(primer_pair_record):
    regions = [(primer_pair_record.chrom,
                primer_pair_record.sense_start,
                primer_pair_record.antisense_start)]
    return pybedtools.BedTool(regions)


def get_reads_with_coords(alignment_file, chrom, start, end):
    read_l = []
    read_total = 0
    forw_total = 0
    rev_total = 0
    fmatch_total = 0
    rmatch_total = 0
    print("coords: {0}\t{1}\n".format(start, end))
    for read in alignment_file.fetch(chrom, start, end):
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


def create_bam_from_reads_and_template(read_list, template_bam, new_fn):
    out_bam = pysam.AlignmentFile(new_fn, "wh", template=template_bam)
    for read in read_list:
        out_bam.write(read)
    out_bam.close()


def get_key_from_read(read):
    if (read.flag & 16 == 0):
        return (read.reference_name, read.reference_start, '+')
    else:
        return (read.reference_name, read.reference_end, '-')



def main():
    if (len(sys.argv) != 4):
        print("usage: {0} [thunderbolts manifest] [bam] [outfolder]".format(os.path.basename(sys.argv[0])))
        sys.exit()

    bam_fn = sys.argv[2]
    outfolder = sys.argv[3]

    thunderbolt_manifest = sys.argv[1]
    primer_pos_d = parse_thunderbolts_manifest(thunderbolt_manifest)
    #print(primer_pos_d)

    bamfile = pysam.AlignmentFile(bam_fn, 'rb')
    outfile = outfolder+"/"+"subset.sam"
    out_bam = pysam.AlignmentFile(outfile, "wh", template=bamfile)
    iter = bamfile.fetch()
    count = 0
    hit_count = 0
    sense_count = 0
    for r in iter:
        count += 1
        read_key = get_key_from_read(r)
        try:
            primer_pair_record = primer_pos_d[read_key]
            hit_count += 1
            if read_key[2] == '+':
                sense_count += 1
                original_cigar = r.cigarstring
                original_ref_start = int(r.reference_start)
                if original_cigar:
                    sense_clip_len = len(primer_pair_record.sense_sequence)
                    antisense_clip_len = len(primer_pair_record.antisense_sequence)
                    new_cigar = Cigar(original_cigar).mask_left(sense_clip_len).mask_right(antisense_clip_len).cigar
                    r.cigarstring = new_cigar
                    r.reference_start = original_ref_start + sense_clip_len
            elif read_key[2] == '-':
                original_cigar = r.cigarstring
                original_ref_start = int(r.reference_start)
                if original_cigar:
                    antisense_clip_len = len(primer_pair_record.antisense_sequence)
                    sense_clip_len = len(primer_pair_record.sense_sequence)
                    new_cigar = Cigar(original_cigar).mask_left(sense_clip_len).mask_right(antisense_clip_len).cigar
                    r.cigarstring = new_cigar
                    r.reference_start = original_ref_start + (sense_clip_len - 8)
            out_bam.write(r)
        except KeyError:
            #This read not an amplicon
            pass
    out_bam.close()
    print("total: {}\nhit: {}\nsense: {}".format(count, hit_count, sense_count))
    sys.exit()


    selected_reads = get_reads_with_coords(bamfile, chrom, start-1, end)
    print("Total reads: {0}".format(len(selected_reads)))
    
    create_bam_from_reads_and_template(selected_reads, bamfile, outfile)

    print("done.")

if __name__ == '__main__':
    main()

