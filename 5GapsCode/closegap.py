#!/usr/env/python

import re
import numpy
import math
import argparse


def main():
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group("required arguments")
    required.add_argument("-a", "--annotation", required=True, help="The annotation file to be used (expected .gtf format)")
    required.add_argument("-b", "--bed", required=True, help=".bed file containing the peaks from sequencing data")
    args = parser.parse_args()

    ann_file = open(args.annotation, 'r')
    seq_file = open(args.bed, 'r')
    results_file = open("merged.bed", 'w')

    ann_list, seq_list = parse(ann_file, seq_file)
    ann_chroms, seq_chroms = get_chroms(ann_list, seq_list)
    diffs = get_differences(ann_chroms, seq_chroms)

    # Determine threshold gap size
    numpyarr = numpy.array(diffs)
    numpyarr_sorted = numpy.sort(numpyarr)
    arrlen = numpy.size(numpyarr_sorted)
    threshold_index = int(math.floor(0.1*arrlen)) # 10% threshold
    threshold_gapsize = numpyarr_sorted[threshold_index]
    print "Threshold gap size is: {}".format(threshold_gapsize)

    # Merge the peaks according to the threshold gap size
    merged_list, processed_counter = join_gaps(seq_chroms, threshold_gapsize)
    print "Finished processing {} peaks!".format(processed_counter)

    # Write results to file
    for item in merged_list:
        start, end, chrom = item
        results_file.write("{}\t{}\t{}\n".format(chrom, start, end))

    ann_file.close()
    seq_file.close()
    results_file.close()


# Input: Annotation file (ann_file), Sequence/.bed file (seq_file)
# Output: Two lists containing tuples of (start, end, +/-), one for each file
def parse(ann_file, seq_file):
    ann_list = []
    seq_list = []
    ann_expr = re.compile("(chr[0-9XYM][0-9XYM]*)\t.*\t([0-9]*)\t([0-9]*)\t.*\t([+-])")
    seq_expr = re.compile("(chr[0-9X]*)\t([0-9]*)\t([0-9]*)")

    # Parse annotation file
    for line in ann_file:
        searchobj = re.search(ann_expr, line)
        chrom = searchobj.group(1)
        start = int(searchobj.group(2))
        end = int(searchobj.group(3))
        orientation = searchobj.group(4)

        ret_tup = (start, end, chrom, orientation)
        ann_list.append(ret_tup)

    for line in seq_file:
        searchobj = re.search(seq_expr, line)
        chrom = searchobj.group(1)
        start = int(searchobj.group(2))
        end = int(searchobj.group(3))
        
        ret_tup = (start, end, chrom)
        seq_list.append(ret_tup)
    
    return ann_list, seq_list


# Input: Parsed lists from annotation and sequence files
# Output: Dictionaries that map chromosome number to peak/annotation information
def get_chroms(ann_list, seq_list):
    # Separate into distinct dictionaries with chromosome keys and start/end
    # position tuples as values
    ann_chroms = {}
    seq_chroms = {}
    for item in ann_list:
        ann_s, ann_e, ann_c, _ = item
        if ann_chroms.get(ann_c) == None:
            ann_chroms[ann_c] = []
            ann_chroms[ann_c].append((ann_s, ann_e))
        else:
            ann_chroms[ann_c].append((ann_s, ann_e))
    for item in seq_list:
        seq_s, seq_e, seq_c = item
        if seq_chroms.get(seq_c) == None:
            seq_chroms[seq_c] = []
            seq_chroms[seq_c].append((seq_s, seq_e, False, seq_c))
        else:
            seq_chroms[seq_c].append((seq_s, seq_e, False, seq_c))

    # Order the distinct lists by start position
    for item in ann_chroms:
        ann_chroms[item] = sorted(ann_chroms[item], key=lambda x:x[0])
    for item in seq_chroms:
        seq_chroms[item] = sorted(seq_chroms[item], key=lambda x:x[0])

    # Note: The tuples in seq_chroms contain extra information needed for
    #       downstream processing: (start, end, within_gene_peak_bool, chr#)
    return ann_chroms, seq_chroms


# Input: Dictionaries containing chromosome positional information
# Output: Average and SD of gaps between peaks within annotated genes
def get_differences(ann_chroms, seq_chroms):
    # Get common chromosome list for easier comparison in the next step
    common_chr = []
    if len(seq_chroms) > len(ann_chroms):
        common_chr = [item for item in seq_chroms if ann_chroms.get(item) != None]
    else:
        common_chr = [item for item in ann_chroms if seq_chroms.get(item) != None]

    diffs = []

    # Do the gene peak finding per chromosome
    for chrom in common_chr:
        chr_ann_list = ann_chroms[chrom]
        chr_seq_list = seq_chroms[chrom]
        ann_index = 0
        seq_index = 0
        ann_size = len(chr_ann_list)
        seq_size = len(chr_seq_list)

        gene_peaks = []

        while True:
            if ann_index >= ann_size or seq_index >= seq_size:
                break

            ann_s, ann_e = chr_ann_list[ann_index]
            seq_s, seq_e, _, seq_c = chr_seq_list[seq_index]

            if seq_s < ann_s:
                seq_index = seq_index + 1
            elif seq_e <= ann_e:
                gene_peaks.append(chr_seq_list[seq_index])
                chr_seq_list[seq_index] = (seq_s, seq_e, True, seq_c)
                seq_index = seq_index + 1
            else:
                for diff in compute_differences(gene_peaks):
                    diffs.append(diff)
                del gene_peaks[:]
                ann_index = ann_index + 1

        if len(gene_peaks) > 0:
            for diff in compute_differences(gene_peaks):
                diffs.append(diff)
    return diffs


# Input: List of peaks
# Output: List of differences
# Note: Let the peaks be A-B-C-D. This function will compute distances AB, BC,
#       and CD. The 3 differences will be returned in a list.
def compute_differences(lst):
    size = len(lst)
    index = 0
    results = []
    while index < size-1:
        _, first_e, _, _ = lst[index]
        last_s, _, _, _ = lst[index+1]
        results.append(last_s - first_e)
        index = index + 1
    return results


# Input: Dictionary containing chromosome keys mapped to peak information
#        (data structure returned from get_chroms)
# Output: list of tuples containing merged peak information
def join_gaps(seq_chroms, threshold):
    total_merged_list = []
    processed_counter = 0

    # iterate through chromosomes
    keys = seq_chroms.keys()
    for key in keys:
        peak_list = seq_chroms[key]
        chr_merged_list = [] # stores tuples (start, end, chrom)
        index = 0
        num_peaks = len(peak_list)

        working_peak = peak_list[0]
        processed_counter = processed_counter + 1

        # process peaks in a linear fashion
        while index < num_peaks-1:
            p1_s, p1_e, within_gene1, _ = working_peak
            p2_s, p2_e, within_gene2, _ = peak_list[index+1]
            processed_counter = processed_counter + 1
            gap_d = p2_s - p1_e
            if within_gene1 == True and within_gene2 == True:
                working_peak = join(working_peak, peak_list[index+1])
            elif gap_d <= threshold:
                working_peak = join(working_peak, peak_list[index+1])
            else:
                wpeak_s, wpeak_e, _, chrom = working_peak
                chr_merged_list.append((wpeak_s, wpeak_e, chrom))
                working_peak = peak_list[index+1]

            index = index + 1

        wpeak_s, wpeak_e, _, chrom = working_peak
        chr_merged_list.append((wpeak_s, wpeak_e, chrom))

        # Add merged chromosome peaks to the total list
        for item in chr_merged_list:
            total_merged_list.append(item)

    return total_merged_list, processed_counter


# Input: 2 peak tuples (start, end, within_gene_peak_bool, chrom)
# Output: 1 working peak tuple (start, end, within_gene_peak_bool, chrom)
# Note: Assumes we are linking end of peak1 to start of peak2. Assumes peaks are
#       from the same chromosome
def join(peak1, peak2):
    p1_s, _, _, chrom = peak1
    _, p2_e, _, _ = peak2

    return (p1_s, p2_e, False, chrom)


if __name__ == "__main__":
    main()
