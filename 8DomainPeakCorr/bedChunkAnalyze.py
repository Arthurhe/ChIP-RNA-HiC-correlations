#!/usr/bin/env python3

# Description: This script will take in two bed files (A, B) and output a
#              comparison matrix. Each interval in A will be normalized by
#              size and divided into n bins. Each bin will then contain a
#              coverage value describing the coverage of that region by
#              intervals in B.

import argparse
import subprocess
import sys
import numpy


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("bedA")
    parser.add_argument("bedPath")
    parser.add_argument("markerListFile")
    parser.add_argument("numBins")
    args = parser.parse_args()
    bins = int(args.numBins)

    bedA = open(args.bedA, 'r')
    tmpFile = open("splitBed.bed", 'w')

    # Divide every peak into chunks
    for line in bedA:
        chunk(line, bins, tmpFile)

    tmpFile.close()

    # I have the split bed file - compute the coverages of that file
    subprocess.call(["python3", "bedtoolsCoverageWrapper.py", "splitBed.bed",
                     args.bedPath, args.markerListFile, "-o", "splitCoverages.bed"])

    # Construct matrix from each additional column in the coverage file
    sCovFile = open("splitCoverages.bed", 'r')
    result = parse_coverage_file(sCovFile)
    sCovFile.close()

    # Create matrices from parsed results
    for i in range(len(result)):
        numSubBins = len(result[i][1])
        numBins = numSubBins // bins
        name = result[i][0]
        a = numpy.asarray(result[i][1])
        b = a.reshape((numBins, -1))
        print(name)
        print(b)

    bedA.close()
    return


# Divide each given bed peak into n bins of even size, and write bins to file
def chunk(line, n, outFile):
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
            outFile.write("{}\t{}\t{}\n".format(chrom, start + i*d, end))
        else:
            outFile.write("{}\t{}\t{}\n".format(chrom, start + i*d,
                                                start + (i+1)*d-1))
        i += 1

    return


def parse_coverage_file(covFile):
    sys.stderr.write("Parsing the coverage data...\n")
    header = covFile.readline()
    header = header.rstrip()
    numTab = header.count("\t")
    numMarkers = numTab - 2
    sheader = header.split("\t")
    coverages = [("", [])] * numMarkers

    for i in range(3, 3+numMarkers):
        coverages[i-3] = (sheader[i], [])

    for line in covFile:
        line = line.rstrip()
        sline = line.split("\t")
        for i in range(3, 3+numMarkers):
            coverages[i-3][1].append(sline[i])

    return coverages


if __name__ == "__main__":
    main()
