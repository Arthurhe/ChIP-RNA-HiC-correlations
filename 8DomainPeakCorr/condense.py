#!/usr/bin/env python3

# Description: Take in a "multi-feature" matrix where each feature is M columns
#              wide and condense all features into a matrix of total width M.
#              This condensed matrix then provides an easy way to visualize
#              the relationships between each feature in the multi-feature
#              matrix.

import argparse
import numpy
import copy


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("matrFile", help="File containing the matrix")
    parser.add_argument("N", type=int, help="Number of features in the matrix")
    # parser.add_argument("colorFile", help="File containing feature color associations")
    args = parser.parse_args()

    n = args.N

    origMatr = numpy.loadtxt(args.matrFile)
    # origMatr = numpy.loadtxt("tmp.txt")
    a = condense(origMatr, n)

    print(a)
    return


def condense(origNumpyMatr, n):
    subarrays = numpy.hsplit(origNumpyMatr, n)
    subarraysCpy = copy.deepcopy(subarrays)

    tmpTotal = None
    for subarray in subarraysCpy:
        subarray[subarray > 0] = 1
        if tmpTotal is None:
            tmpTotal = subarray
        else:
            tmpTotal += subarray

    print(tmpTotal)
    a, b = numpy.where(tmpTotal > 1)

    # Do comparisons for overlapping positions
    changes = []
    for i in range(len(a)):
        j = 1
        maxval = 0
        maxj = 0
        for subarray in subarrays:
            val = subarray[a[i], b[i]]
            if val > maxval:
                maxj = j
                maxval = val
            j += 1

        # We found the max val at a given position among all subarrays
        changes.append((a[i], b[i], maxj))

    # Condense the subarrays into one final array
    counter = 1
    total = None
    for subarray in subarrays:
        subarray[subarray > 0] = counter
        counter += 1
        if total is None:
            total = subarray
        else:
            total = total + subarray

    # Apply changes to condensed array (maxval in overlapping positions)
    for change in changes:
        a, b, val = change
        total[a, b] = val

    return total


if __name__ == "__main__":
    main()
