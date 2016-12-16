#!/usr/bin/env python

import numpy
import argparse


def main():
    # parse input arguments
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group("required argument")
    required.add_argument("-i", "--infile", required=True, help="The HiC data file to transform")
    args = parser.parse_args()

    # read original matrix from file
    datamatr = numpy.loadtxt(args.infile)
    nrow, ncol = datamatr.shape

    # Transformed matrix zero-initialized
    transmat = numpy.zeros((nrow, ncol))
    tnrow, tncol = datamatr.shape

    # TODO: implement some threshold d beyond which we don't care about it

    # iterate through all elements in the matrix and transform
    it = numpy.nditer(datamatr, flags=['multi_index'])
    while not it.finished:
        currow, curcol = it.multi_index
        d = abs(currow - curcol)
        lindex = currow - d
        rindex = currow + d
        onedindex = (currow * ncol) + curcol

        # only compute if both left and right indices exist
        if lindex >= 0 and rindex < ncol:
            lval = datamatr.item(currow, lindex)
            rval = datamatr.item(currow, rindex)
            numer = lval - rval
            denom = lval + rval
            if denom == 0:
                denom = 1
            tval = numer / denom
            numpy.put(transmat, onedindex, tval)

        it.iternext()

    # print(transmat)
    # print(transmat.shape)
    numpy.savetxt("transformed.txt", transmat, fmt='%.13f', delimiter="\t")
    prepend("transformed.txt", str(transmat.shape[0]))

    squared = transmat * transmat
    # print(squared)
    # print(squared.shape)
    numpy.savetxt("squared.txt", squared, fmt='%.13f', delimiter="\t")
    prepend("squared.txt", str(squared.shape[0]))


def prepend(filename, line):
    '''
    Function to add the matrix dimensions to the beginning of the file, as
    expected by the scoring program
    '''
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)

if __name__ == "__main__":
    main()
