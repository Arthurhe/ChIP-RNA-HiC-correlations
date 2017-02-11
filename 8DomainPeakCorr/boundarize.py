#!/usr/bin/env python3

# This script will take in a bed file (specifically, a bed file describing the
# called domains of a chromosome/genome) and shift the boundaries such that
# the new boundaries surround an original boundary edge.
# Example:
# original: start: 00   end: 10
#           start: 10   end: 20
# reworked: start: 00   end: 05
#           start: 05   end: 15
# The purpose is so that we can use bedtools to compute intersections between
# ChIP peaks and these new boundary edges to see which ChIP peaks are
# associated with a given domain boundary.

import argparse


def main():
    # Parse input
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group("required")
    required.add_argument("-i", "--infile", required=True)
    args = parser.parse_args()

    for asdf in boundarize(args.infile):
        print(asdf)


def boundarize(fname):
    infile = open(fname, 'r')
    prevLine = None

    for line in infile:
        line = line.rstrip()
        sline = line.split("\t")

        # we want integer division
        sline[2] = (int(sline[1]) + int(sline[2])) // 2
        sline[2] = str(sline[2])

        if prevLine is not None:
            sline[1] = prevLine[2]

        prevLine = sline

        newLine = '\t'.join(sline)
        yield newLine

    pass


if __name__ == "__main__":
    main()
