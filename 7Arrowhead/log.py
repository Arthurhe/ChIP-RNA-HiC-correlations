#!/usr/bin/env python

# Program to take every value in the given matrix and perform x = log_10(x)

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

    # do the log in a vectorized way (?)
    transmat = numpy.log(datamatr+1)

    print(datamatr)
    print(transmat)
    # Note, you may need to adjust the float precision of the output to save
    # save space.
    numpy.savetxt("log10.txt", transmat, fmt='%.13f', delimiter="\t")

if __name__ == "__main__":
    main()
