#!/usr/env/python

import re


def main():
    ann_file = open("/home/brian/genes_cpy.txt", 'r')
    seq_file = open("./data/GM12_ENCFF439LIW_rRNA_minus_RNA_calledRegions_peaks.bed")

    #ann_list, seq_list = parse(ann_file, seq_file)
    ann_list = [(2, 7), (11, 15), (33, 70)]
    seq_list = [(3, 4), (4, 5), (8, 9), (13, 14), (15, 17), (22, 27), (40, 50)]
    get_gene_peaks(ann_list, seq_list)
    #compute_stats(ann_list, seq_list)

    #for item in ann_list:
    #    print item

    #for item in seq_list:
    #    print item

    ann_file.close()
    seq_file.close()


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

def get_gene_peaks(ann_list, seq_list):
    ann_size = len(ann_list)
    seq_size = len(seq_list)

    ann_index = 0
    seq_index = 0

    while(True):
        if ann_index >= ann_size or seq_index >= seq_size:
            break

        # Unpack the lists
        seq_s, seq_e
        print "a: {}\ts: {}".format(ann_list[ann_index], seq_list[seq_index])

        ann_index = ann_index + 1
        seq_index = seq_index + 1
    pass


if __name__ == "__main__":
    main()
