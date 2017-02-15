#!/usr/bin/env python3

import argparse


def main():
    # Handle input arguments
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group()
    required.add_argument("-o", "--original", required=True)
    required.add_argument("-s", "--shifted", required=True)
    required.add_argument("-c", "--correlated", required=True)
    args = parser.parse_args()

    shiftedBed = open(args.shifted, 'r')
    origBed = open(args.original, 'r')
    corrBed = open(args.correlated, 'rw')

    # Read in shifted boundaries, create dict of ((start, end), line #) tuples
    shiftedDict = read_shifted_bed(shiftedBed)

    # Read in original boundaries, create dict of (line #, (start, end)) tuples
    origDict = read_orig_bed(origBed)

    # Read intersected bed file and replace boundaries
    for line in corrBed:
        line = line.rstrip()
        sline = line.split("\t")
        start = sline[1]
        end = sline[2]
        coord = "{} {}".format(start, end)
        origLineNum = shiftedDict[coord]
        origPos = origDict[origLineNum]
        origPos = origPos.split()
        sline[1] = origPos[0]
        sline[2] = origPos[1]
        newStr = '\t'.join(sline)
        print(newStr)

    shiftedBed.close()
    origBed.close()
    corrBed.close()
    pass


def read_shifted_bed(bedFile):
    dictionary = {}
    counter = 1
    for line in bedFile:
        line = line.rstrip()
        sline = line.split("\t")
        start = sline[1]
        end = sline[2]
        coord = "{} {}".format(start, end)
        dictionary[coord] = counter
        counter += 1

    return dictionary


def read_orig_bed(bedFile):
    dictionary = {}
    counter = 1
    for line in bedFile:
        line = line.rstrip()
        sline = line.split("\t")
        start = sline[1]
        coord = "{} {}".format(start, start)
        dictionary[counter] = coord
        counter += 1

    return dictionary


if __name__ == "__main__":
    main()
