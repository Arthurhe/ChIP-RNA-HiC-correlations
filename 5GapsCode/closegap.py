#!/usr/env/python

import re
import numpy
import math
import argparse
import Queue


def main():
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group("required arguments")
    required.add_argument("-a", "--annotation", required=True, help="The annotation file to be used (expected .gtf format)")
    required.add_argument("-b", "--bed", required=True, help=".bed file containing the peaks from sequencing data")
    args = parser.parse_args()

    ann_file = open(args.annotation, 'r')
    seq_file = open(args.bed, 'r')
    #results_file = open("merged.bed", 'w')

    ann_list, seq_list, peak_counter = parse(ann_file, seq_file)
    #print "Total number of peaks in .bed file: {}".format(peak_counter)
    ann_chroms, seq_chroms = get_chroms(ann_list, seq_list)
    diffs = get_differences(ann_chroms, seq_chroms)

    # Determine threshold gap size
    numpyarr = numpy.array(diffs)
    numpyarr_sorted = numpy.sort(numpyarr)
    arrlen = numpy.size(numpyarr_sorted)
    threshold_index = int(math.floor(0.1 * arrlen))  # 10% threshold
    threshold_gapsize = numpyarr_sorted[threshold_index]
    #print "Threshold gap size is: {}".format(threshold_gapsize)

    # Merge the peaks according to the threshold gap size
    merged_list, processed_counter = join_gaps(seq_chroms, threshold_gapsize)
    #print "Finished processing {} peaks!".format(processed_counter)

    #if processed_counter == peak_counter:
    #    print "All peaks have been acounted for."
    #else:
    #    print "CAUTION: Some peaks were skipped!!"

    # Write results to file
    for item in merged_list:
        start, end, chrom = item
        print "{}\t{}\t{}".format(chrom, start, end)
        #results_file.write("{}\t{}\t{}\n".format(chrom, start, end))

    ann_file.close()
    seq_file.close()
    #results_file.close()


# Input: Annotation file (ann_file), Sequence/.bed file (seq_file)
# Output: Two lists containing tuples of (start, end, +/-), one for each file
def parse(ann_file, seq_file):
    peak_counter = 0
    ann_list = []
    seq_list = []
    ann_expr = re.compile("(chr[0-9XYM][0-9XYM]*)\t.*\t([0-9]*)\t([0-9]*)\t.*\t([+-])\t.\tgene_id \"(.*)\"; transcript_id")  # don't ask...
    seq_expr = re.compile("(chr[0-9X]*)\t([0-9]*)\t([0-9]*)")

    # Parse annotation file
    for line in ann_file:
        searchobj = re.search(ann_expr, line)
        chrom = searchobj.group(1)
        start = int(searchobj.group(2))
        end = int(searchobj.group(3))
        orientation = searchobj.group(4)
        gene_id = searchobj.group(5)

        ret_tup = (start, end, chrom, orientation, gene_id)
        ann_list.append(ret_tup)

    for line in seq_file:
        peak_counter = peak_counter+1
        searchobj = re.search(seq_expr, line)
        chrom = searchobj.group(1)
        start = int(searchobj.group(2))
        end = int(searchobj.group(3))

        ret_tup = (start, end, chrom)
        seq_list.append(ret_tup)

    return ann_list, seq_list, peak_counter


# Input: Parsed lists from annotation and sequence files
# Output: Dictionaries that map chromosome number to peak/annotation information
def get_chroms(ann_list, seq_list):
    # Separate into distinct dictionaries with chromosome keys and start/end
    # position tuples as values
    ann_chroms = {}
    seq_chroms = {}
    for item in ann_list:
        ann_s, ann_e, ann_c, _, ann_id = item
        if ann_chroms.get(ann_c) == None:
            ann_chroms[ann_c] = []
            ann_chroms[ann_c].append((ann_s, ann_e, ann_id))
        else:
            ann_chroms[ann_c].append((ann_s, ann_e, ann_id))
    for item in seq_list:
        seq_s, seq_e, seq_c = item
        if seq_chroms.get(seq_c) == None:
            seq_chroms[seq_c] = []
            seq_chroms[seq_c].append((seq_s, seq_e, [], seq_c))
        else:
            seq_chroms[seq_c].append((seq_s, seq_e, [], seq_c))

    # Order the distinct lists by start position
    for item in ann_chroms:
        ann_chroms[item] = sorted(ann_chroms[item], key=lambda x:x[0])
    for item in seq_chroms:
        seq_chroms[item] = sorted(seq_chroms[item], key=lambda x:x[0])

    # Note: The tuples in seq_chroms contain extra information needed for
    #       downstream processing: (start, end, [gene_id], chr#)
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

            ann_s, ann_e, ann_id = chr_ann_list[ann_index]
            seq_s, seq_e, _, seq_c = chr_seq_list[seq_index]

            if seq_s < ann_s:
                seq_index = seq_index + 1
            elif seq_e <= ann_e:
                gene_peaks.append(chr_seq_list[seq_index])
                id_list = []
                id_list.append(ann_id)

                # Check the later annotated genes for possible membership
                tmp_idx = ann_index + 1
                while tmp_idx < ann_size:
                    next_ann_s, next_ann_e, next_ann_id = chr_ann_list[tmp_idx]
                    if seq_s >= next_ann_s and seq_e <= next_ann_e:
                        id_list.append(next_ann_id)
                        tmp_idx = tmp_idx+1
                    else:
                        break

                chr_seq_list[seq_index] = (seq_s, seq_e, id_list, seq_c)

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
        #print key
        peak_list = seq_chroms[key]
        chr_merged_list = [] # stores tuples (start, end, chrom)
        index = 0
        num_peaks = len(peak_list)

        working_peak = peak_list[0]
        processed_counter = processed_counter + 1

        prev_stack = Queue.LifoQueue()
        tmp_stack = Queue.LifoQueue()

        # process peaks in a linear fashion
        while index < num_peaks-1:
            p1_s, p1_e, p1_gID, _ = working_peak
            p2_s, p2_e, p2_gID, _ = peak_list[index+1]
            processed_counter = processed_counter + 1
            gap_d = p2_s - p1_e
            numgeneID1 = len(p1_gID)
            numgeneID2 = len(p2_gID)

            # No ambiguous gene membership
            if numgeneID1 == 1 and numgeneID2 == 1:
                # same gene - merge no questions asked
                if p1_gID[0] == p2_gID[0]:
                    gID = p1_gID
                    working_peak = join(working_peak, peak_list[index+1], gID)
                else:
                    # different genes - merge depending on threshold
                    if gap_d <= threshold:
                        gID = p2_gID
                        working_peak = join(working_peak, peak_list[index+1], gID)
                    # dump working peak and reset
                    else:
                        wpeak_s, wpeak_e, _, chrom = working_peak
                        chr_merged_list.append((wpeak_s, wpeak_e, chrom))
                        working_peak = peak_list[index+1]
            # Ambiguous gene membership - do merging based on threshold
            else:
                # Case: unambiguous gene, ambiguous gene
                # put on prev_stack and reset working peak
                if numgeneID1 == 1 and numgeneID2 > 1:
                    prev_stack.put(working_peak)
                    working_peak = peak_list[index+1]
                # Case: ambiguous gene, unambiguous gene
                # continuously pop from prev_stack looking for the same gene
                # if same gene is found, merge everything and move on
                elif numgeneID1 > 1 and numgeneID2 == 1:
                    bool_merged = False
                    while not prev_stack.empty():
                        prev_peak = prev_stack.get()
                        _, _, prev_peak_IDS, _ = prev_peak
                        if prev_peak_IDS == p2_gID:
                            bool_merged = True
                            gID = p2_gID
                            working_peak = join(prev_peak, working_peak, gID)
                            working_peak = join(working_peak, peak_list[index+1], gID)
                        else:
                            tmp_stack.put(prev_peak)
                    assert prev_stack.empty() is True, "items still remain on prev_stack"
                # if same gene is not found, dump all peaks to file and check
                # if the immediatly previous peak is within the threshold
                    if bool_merged is False:
                        if tmp_stack.empty():
                            if gap_d <= threshold:
                                gID = p2_gID
                                working_peak = join(working_peak, peak_list[index+1], gID)
                            # dump working peak and reset
                            else:
                                wpeak_s, wpeak_e, _, chrom = working_peak
                                chr_merged_list.append((wpeak_s, wpeak_e, chrom))
                                working_peak = peak_list[index+1]
                        else:
                            # check if we can link by threshold
                            while tmp_stack.qsize() > 1:
                                s, e, _, chrom = tmp_stack.get()
                                chr_merged_list.append((s, e, chrom))
                            adjacent_peak = tmp_stack.get()
                            assert tmp_stack.empty() is True, "items still remain on tmp_stack"
                            if gap_d <= threshold:
                                gID = p2_gID
                                working_peak = join(working_peak, peak_list[index+1], gID)
                            # dump working peak and reset
                            else:
                                wpeak_s, wpeak_e, _, chrom = working_peak
                                chr_merged_list.append((wpeak_s, wpeak_e, chrom))
                                working_peak = peak_list[index+1]

                # Case: ambiguous gene, ambiguous gene
                # same ambiguous membership
                elif numgeneID1 > 1 and numgeneID2 > 1:
                    if p1_gID == p2_gID:
                        gID = p1_gID
                        working_peak = join(working_peak, peak_list[index+1], gID)
                    # differing ambiguous membership
                    else:
                        if gap_d <= threshold:
                            gID = p2_gID
                            working_peak = join(working_peak, peak_list[index+1], gID)
                        # dump working peak and reset
                        else:
                            wpeak_s, wpeak_e, _, chrom = working_peak
                            chr_merged_list.append((wpeak_s, wpeak_e, chrom))
                            working_peak = peak_list[index+1]

                # Case: anything paired with no gene
                # Handle all with threshold
                else:
                    if gap_d <= threshold:
                        gID = p2_gID
                        working_peak = join(working_peak, peak_list[index+1], gID)
                    # dump working peak and reset
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


# Input: 2 peak tuples (start, end, [geneIDs], chrom)
# Output: 1 working peak tuple (start, end, [geneIDs], chrom)
# Note: Assumes we are linking end of peak1 to start of peak2. Assumes peaks are
#       from the same chromosome
def join(peak1, peak2, gID):
    p1_s, _, _, chrom = peak1
    _, p2_e, _, _ = peak2

    return (p1_s, p2_e, gID, chrom)


if __name__ == "__main__":
    main()
