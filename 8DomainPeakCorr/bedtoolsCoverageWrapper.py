#!/usr/bin/env python3

import subprocess
import argparse
import os


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("bedA")
    parser.add_argument("bedPath")
    args = parser.parse_args()

    wrapper(args.bedA, args.bedPath)


def wrapper(fa, bedPath):
    # TODO: do filename filtering based on input list
    # FIXME: for now, just call coverage on all files

    # buff stores the coverage results
    buff = []
    for fb in os.listdir(bedPath):
        print("Computing converage of: {}".format(fb))
        buff = call_bedtools_filter(fa, "{}/{}".format(bedPath, fb), buff)

    # Add the 'chr start end' columns to buff
    fileA = open(fa, 'r')
    buff[0] = "chrom\tstart\tend\t" + buff[0]
    i = 1
    for line in fileA:
        line = line.rstrip()
        sline = line.split("\t")
        buff[i] = sline[0] + "\t" + sline[1] + "\t" + sline[2] + "\t" + buff[i]
        i += 1

    for thing in buff:
        print(thing)


# Bedtools works with filenames, not file descriptors
def call_bedtools_filter(fa, fb, buff):
    # Call bedtools coverage
    proc = subprocess.Popen(["bedtools", "coverage", "-a", fa, "-b", fb],
                            stdout=subprocess.PIPE)

    # Buff is essentially one column. This column contains as many rows as
    # there are peaks (rows) in bedA (+1). Each row in Buff will contain a tab
    # delimited string describing the coverage of that peak in bedA. The
    # first row is always a header describing the corresponding marker

    # Format the header
    if len(buff) == 0:
        buff.append(fb)
    else:
        buff[0] = buff[0] + "\t" + fb

    i = 1
    # Filter and format output
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
