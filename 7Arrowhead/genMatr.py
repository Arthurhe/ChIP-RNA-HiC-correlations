#!/usr/bin/env python

import numpy
import sys


def main():
    matrsize = 500
    domstart = 200
    domend = 300

    for r in range(matrsize):
        for c in range(matrsize):
            if r >= domstart and r <= domend and c >= domstart and c <= domend:
                num = numpy.random.normal(10.0, 0.1)
                sys.stdout.write(str(num) + "\t")
            else:
                num = numpy.random.normal(1.0, 0.1)
                sys.stdout.write(str(num) + "\t")
        sys.stdout.write("\n")

if __name__ == "__main__":
    main()
