#!/usr/bin/env python3

# Description: This script will take in two bed files (A, B) and output a
#              comparison matrix. Each interval in A will be normalized by
#              size and divided into n bins. Each bin will then contain a
#              coverage value describing the coverage of that region by
#              intervals in B.

import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("bedA")
    parser.add_argument("bedB")
    args = parser.parse_args()

    print(args.bedA)
    bedA = open(args.bedA, 'r')
    tmpFile = open("tmp.txt", 'w')

    peakId = 0
    for line in bedA:
        intervalSize = chunk(line, 14, tmpFile, peakId)
        print(intervalSize)
        peakId += 1

    tmpFile.close()

    # Call bedtools

    # TODO: delete outfile

    # TODO: divide by intervalSize to normalize by size

    bedA.close()
    pass


# Divide each given bed peak into n bins of even size, and write bins to file
def chunk(line, n, outFile, pId):
    line = line.rstrip()
    sline = line.split("\t")

    chrom = sline[0]
    start = int(sline[1])
    end = int(sline[2])

    d = (end - start + 1) // n
    i = 0

    # Handle the case where the bin size doesn't divide evenly by extending the
    # the last bin
    while i < n:
        if i == n-1:
            outFile.write("{}\t{}\t{}\t{}\n".format(chrom, start + i*d, end,
                                                    pId))
        else:
            outFile.write("{}\t{}\t{}\t{}\n".format(chrom, start + i*d,
                                                    start + (i+1)*d-1, pId))
        i += 1

    # TODO: not actually necessary
    return d


if __name__ == "__main__":
    main()
