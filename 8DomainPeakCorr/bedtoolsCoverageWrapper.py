#!/usr/bin/env python3

import subprocess
import argparse
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("bedA", help="bed to compare against (eg. domain file)")
    parser.add_argument("markerListFile", help="tab-separated-file describing the marker and a path name to that marker's bed file")
    parser.add_argument("-o", "--output", help="file to write output to. If not specified, print results to stdout")
    args = parser.parse_args()

    # Parse the markerListFile
    # The markerListFile is a 2-column, tab-separated file describing the tag
    # name and that tag's bed filepath. Eg.
    # CTCF  /some/path/to/ctcf.bed
    # SMC   /some/path/to/smc.bed
    # H3K4  /some/path/to/H3K4.bed

    f = open(args.markerListFile, 'r')

    # markers is a list of marker name/filepath tuples
    markers = []
    for line in f:
        line = line.rstrip()
        sline = line.split("\t")
        markers.append((sline[0], sline[1]))

    f.close()

    result = wrapper(args.bedA, args.bedPath, markers)

    # Decide how to return output
    if args.output is None:
        for thing in result:
            print(thing)
    else:
        outfile = open(args.output, 'w')
        for thing in result:
            outfile.write("{}\n".format(thing))
        outfile.close()


def wrapper(fa, bedPath, markers):
    # buff stores the coverage results
    buff = []
    for tupl in markers:
        markerName, path = tupl
        sys.stderr.write("Computing converage of: {}\n".format(path))
        buff = call_bedtools_filter(fa, path, buff, markerName)

    # Output formatting - Add the 'chr start end' columns to buff
    fileA = open(fa, 'r')
    buff[0] = "chrom\tstart\tend\t" + buff[0]
    i = 1
    for line in fileA:
        line = line.rstrip()
        sline = line.split("\t")
        buff[i] = sline[0] + "\t" + sline[1] + "\t" + sline[2] + "\t" + buff[i]
        i += 1

    return buff


# Bedtools works with filenames, not file descriptors
def call_bedtools_filter(fa, fb, buff, markerName):
    # Call bedtools coverage
    proc = subprocess.Popen(["bedtools", "coverage", "-a", fa, "-b", fb],
                            stdout=subprocess.PIPE)

    # Buff is essentially one column. This column contains as many rows as
    # there are peaks (rows) in bedA (+1). Each row in Buff will contain a tab
    # delimited string describing the coverage of that peak in bedA. The
    # first row is always a header describing the corresponding marker

    # Format the header
    if len(buff) == 0:
        buff.append(markerName)
    else:
        buff[0] = buff[0] + "\t" + markerName

    # Filter and format output
    i = 1
    while True:
        line = proc.stdout.readline()
        line = line.decode("utf-8")
        line = line.rstrip()
        sline = line.split("\t")

        if len(line) != 0:
            if len(buff) > i:
                buff[i] = buff[i] + "\t" + sline[-1]
            else:
                buff.append(sline[-1])
        else:
            break
        i += 1

    return buff


if __name__ == "__main__":
    main()
