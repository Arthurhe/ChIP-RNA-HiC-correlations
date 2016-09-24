#!/usr/env/python

import re


def main():
    ann_file = open("./data/genes_cpy.txt", 'r')
    seq_file = open("./data/RNA/GM12_ENCFF439LIW_rRNA_minus_RNA_calledRegions_peaks.bed")

    ann_list, seq_list = parse(ann_file, seq_file)
    #compute_stats(ann_list, seq_list)

    for item in ann_list:
        print item

    for item in seq_list:
        print item

    ann_file.close()
    seq_file.close()


# Input: Annotation file (ann_file), Sequence/.bed file (seq_file)
# Output: Two lists containing tuples of (start, end, +/-), one for each file
def parse(ann_file, seq_file):
    ann_list = []
    seq_list = []
    ann_expr = re.compile("\t([0-9]*)\t([0-9]*)\t.*\t([+-])")
    seq_expr = re.compile("(chr[0-9X]*)\t([0-9]*)\t([0-9]*)")

    # Parse annotation file
    for line in ann_file:
        searchobj = re.search(ann_expr, line)
        start = int(searchobj.group(1))
        end = int(searchobj.group(2))
        orientation = searchobj.group(3)

        ret_tup = (start, end, orientation)
        ann_list.append(ret_tup)

    for line in seq_file:
        searchobj = re.search(seq_expr, line)
        chrom = searchobj.group(1)
        start = int(searchobj.group(2))
        end = int(searchobj.group(3))
        
        ret_tup = (start, end, chrom)
        seq_list.append(ret_tup)
    
    return ann_list, seq_list


if __name__ == "__main__":
    main()
