#!/usr/env/python

import re


def main():
    #ann_file = open("/home/brian/genes_cpy.txt", 'r')
    #seq_file = open("./data/GM12_ENCFF439LIW_rRNA_minus_RNA_calledRegions_peaks.bed")

    #ann_list, seq_list = parse(ann_file, seq_file)
    ann_list = [(11, 15, "chr1"), (2, 7, "chr1"), (33, 70, "chr1"), (15, 500, "chr2")]
    seq_list = [(3, 4, "chr1"), (4, 5, "chr1"), (8, 9, "chr1"), (13, 14, "chr1"), (15, 17, "chr1"), (22, 27, "chr1"), (40, 50, "chr1"), (333, 667, "chr2"), (42, 56, "chr3")]
    ann_chroms, seq_chroms = get_chroms(ann_list, seq_list)
    compute_stats(ann_chroms, seq_chroms)
    #compute_stats(ann_list, seq_list)

    #for item in ann_list:
    #    print item

    #for item in seq_list:
    #    print item

    #ann_file.close()
    #seq_file.close()


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
        ann_s, ann_e, ann_c = item
        if ann_chroms.get(ann_c) == None:
            ann_chroms[ann_c] = []
            ann_chroms[ann_c].append((ann_s, ann_e))
        else:
            ann_chroms[ann_c].append((ann_s, ann_e))
    for item in seq_list:
        seq_s, seq_e, seq_c = item
        if seq_chroms.get(seq_c) == None:
            seq_chroms[seq_c] = []
            seq_chroms[seq_c].append((seq_s, seq_e))
        else:
            seq_chroms[seq_c].append((seq_s, seq_e))

    # Order the distinct lists by start position
    for item in ann_chroms:
        ann_chroms[item] = sorted(ann_chroms[item], key=lambda x:x[0])
    for item in seq_chroms:
        seq_chroms[item] = sorted(seq_chroms[item], key=lambda x:x[0])

    return ann_chroms, seq_chroms

# Input: Dictionaries containing chromosome positional information
# Output: Average and SD of gaps between peaks within annotated genes
def compute_stats(ann_chroms, seq_chroms):
    # Get common chromosome list for easier comparison in the next step
    common_chr = []
    if len(seq_chroms) > len(ann_chroms):
        common_chr = [item for item in seq_chroms if ann_chroms.get(item) != None]
    else:
        common_chr = [item for item in ann_chroms if seq_chroms.get(item) != None]

    # Do the gene peak finding per chromosome
    for chrom in common_chr:
        chr_ann_list = ann_chroms[chrom]
        chr_seq_list = seq_chroms[chrom]

        ann_index = 0
        seq_index = 0
        ann_size = len(chr_ann_list)
        seq_size = len(chr_seq_list)

        while True:
            if ann_index >= ann_size or seq_index >= seq_size:
                break

            ann_s, ann_e = chr_ann_list[ann_index]
            seq_s, seq_e = chr_seq_list[seq_index]

            if seq_s < ann_s:
                seq_index = seq_index + 1
            else:
                if seq_e <= ann_e:
                    print "+1"
                    seq_index = seq_index + 1
                    pass
                else:
                    ann_index = ann_index + 1

            print "asdf"


if __name__ == "__main__":
    main()
