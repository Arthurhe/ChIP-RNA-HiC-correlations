#ifndef __LTRIANGLE_CPP__
#define __LTRIANGLE_CPP__

#include <cmath>
#include <cfloat>
#include <iostream>
#include <cassert>

#include "LTriangle.hpp"
#include "RCS.hpp"
#include "UTriangle.hpp"
#include "LVTriangle.hpp"

LTriangle::LTriangle(int size, RCS& rScores, RCS& cScores) {
    // sanity checks?

    // initialize pg
    pg = new TwoD_Array<float>(size, size);

    // Initialize everythhing to 0: fixes valgrind errors
    for(int i = 0; i < pg->getNumCols(); ++i) {
        for(int j = 0; j < pg->getNumRows(); ++j) {
            pg->at(j, i) = 0;
        }
    }

    // place holders for subtraction in DP recurrence
    int lastTop;
    int lastBot;
    
    // compute
    for(int r = 0; r < pg->getNumRows(); ++r) {
        for(int c = 0; c < pg->getNumCols(); ++c) {
            // compute dimensions of triangle
            int top = ceil(((float)(r + c))/2);
            int bottom = c;
            int right = 2*c - r;

            // recurrence
            if(r < c && right < pg->getNumCols()) {
                pg->at(r, c) = pg->at(r, c-1) +
                               rScores.query(bottom, bottom, right) -
                               cScores.query(c-1, lastTop, lastBot);
            }
            else if (r == c) {
                pg->at(r, c) = rScores.query(bottom, bottom, right);
            }
            else {
                pg->at(r,c) = nanf("");
            }

            lastTop = top;
            lastBot = bottom;
        }
    }
}

// Copy Constructor - make a new copy of the TwoD_array
LTriangle::LTriangle(const LTriangle& orig) {
    pg = new TwoD_Array<float>(orig.pg->getNumRows(), orig.pg->getNumCols());
    for(int r = 0; r < pg->getNumRows(); ++r) {
        for(int c = 0; c < pg->getNumCols(); ++c) {
            pg->at(r, c) = orig.pg->at(r, c);
        }
    }

    maxVal = orig.maxVal;
}

LTriangle::~LTriangle() {
    delete pg;
    pg = 0;
}

void LTriangle::subtract(UTriangle& upper) {
    // TODO: check that upper and lower are the same dimensions
    for(int r = 0; r < this->pg->getNumRows(); ++r) {
        for(int c = 0; c < this->pg->getNumCols(); ++c) {
            if(!std::isnan(pg->at(r, c))) {
                pg->at(r, c) = pg->at(r, c) - upper.getValue(r, c);
                if(r == 0 && c == 0) {
                    maxVal = pg->at(r, c);
                }
                else if (pg->at(r, c) > maxVal) {
                    maxVal = pg->at(r, c);
                }
            }
            else {
                pg->at(r, c) = nanf("");
            }
        }
    }
}

void LTriangle::normalize() {
    // Prevent divide by 0's only if 0 was the max value
    if(maxVal == 0) {
        maxVal = 1;
    }

    for(int r = 0; r < this->pg->getNumRows(); ++r) {
        for(int c = 0; c < this->pg->getNumCols(); ++c) {
            pg->at(r, c) = pg->at(r, c) / maxVal;
        }
    }
}

float LTriangle::getValue(int r, int c) {
    return pg->at(r, c);
}

void LTriangle::display() {
    (*pg).printOut();
}

void LTriangle::setCountPtr(LVTriangle& Lower) {
    lower = &Lower;
}

// Function to divide each number in the matrix by the number of elements in the
// triangle that it represents. This operation is performed in place.
void LTriangle::calcMean() {
    for(int r = 0; r < pg->getNumRows(); ++r) {
        for(int c = 0; c < pg->getNumCols(); ++c) {
            pg->at(r, c) = pg->at(r, c) / lower->getCount(r, c);
        }
    }
}

void LTriangle::write(std::ofstream& fname) {
    for (int r = 0; r < pg->getNumRows(); ++r) {
        for (int c = 0; c < pg->getNumCols(); ++c) {
            float val = pg->at(r, c);
            if(std::isnan(val)) {
                val = -1;
            }
            fname << val << "\t";
        }
        fname << std::endl;
    }
    fname << std::endl;
}

#endif
