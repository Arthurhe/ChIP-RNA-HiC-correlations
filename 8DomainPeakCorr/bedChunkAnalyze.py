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
import os


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("bedA", help="bed to compare against (eg. domain file)")
    parser.add_argument("markerListFile", help="tab-separated-file describing the marker and a path name to that marker's bed file")
    parser.add_argument("numBins", help="Number of bins to divide each peak of bedA into")
    parser.add_argument("outputFormat", choices=['numpyArray', 'bed'], help="Format of the desired output. Specifying 'bed' will output the result of the bedtoolsCoverageWrapper. Specifying 'numpyArray' will reassemble the wrapper output and return a binary numpy array.")
    parser.add_argument("outfname", help="Name of the file output by bedtoolsCoverageWrapper")
    args = parser.parse_args()
    bins = int(args.numBins)

    curDir = os.path.dirname(os.path.abspath(__file__))
    curDir = curDir + "/"
    wrapperFunction = curDir + "bedtoolsCoverageWrapper.py"

    bedA = open(args.bedA, 'r')
    tmpFile = open("splitBed.bed", 'w')

    # Divide every peak into chunks
    for line in bedA:
        chunk(line, bins, tmpFile)

    tmpFile.close()

    # I have the split bed file - compute the coverages of that file
    subprocess.call(["python3", wrapperFunction, "splitBed.bed",
                     args.markerListFile, "-o", args.outfname])

    if args.outputFormat == "numpyArray":
        # Construct matrix from each additional column in the coverage file
        sCovFile = open(args.outfname, 'r')
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
            numpy.save(name, b)

    bedA.close()
    return


# Divide each given bed peak into n bins of even size, and write bins to file
def chunk(line, n, outFile):
    line = line.rstrip()
    sline = line.split("\t")

    chrom = sline[0]
    a = int(sline[1])
    b = int(sline[2])

    start = min(a, b)
    end = max(a, b)

    d = (end - start + 1) / n
    i = 0

    while i < n:
        outFile.write("{}\t{}\t{}\n".format(chrom, round(start + i*d),
                                            round(start + (i+1)*d-1)))
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
            coverages[i-3][1].append(float(sline[i]))

    return coverages


if __name__ == "__main__":
    main()
