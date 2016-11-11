#ifndef __LVTRIANGLE_CPP__
#define __LVTRIANGLE_CPP__

#include <cmath>
#include <cfloat>
#include <iostream>

#include "LVTriangle.hpp"
#include "RCS.hpp"

/**
 * Note: LVTriangle takes in 2 3D column sum matrices and outputs the mean of
 *       all values in the lower triangle.
 */

LVTriangle::LVTriangle(int size, RCS& rScores, RCS& cScores) {
    // sanity checks?

    // initialize pg
    pg = new TwoD_Array<float>(size, size);
    counts = new TwoD_Array<int>(size, size);

    // Initialize everythhing to 0: fixes valgrind errors
    for(int i = 0; i < pg->getNumCols(); ++i) {
        for(int j = 0; j < pg->getNumRows(); ++j) {
            pg->at(j, i) = 0;
            counts->at(j, i) = 0;
        }
    }

    compute(rScores, cScores);
    divide();
}

void LVTriangle::compute(RCS& rScores, RCS& cScores) {
    int lastTop;
    int lastBot;

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
                int removed = lastBot - lastTop + 1;
                int added = right - bottom + 1;
                counts->at(r, c) = counts->at(r, c-1) + added - removed;
            }
            else if (r == c) {
                pg->at(r, c) = rScores.query(bottom, bottom, right);
                counts->at(r, c) = 1;
            }
            else {
                pg->at(r,c) = nanf("");
            }

            lastTop = top;
            lastBot = bottom;
        }
    }
}

void LVTriangle::square() {
    for(int r = 0; r < this->pg->getNumRows(); ++ r) {
        for(int c = 0; c < this->pg->getNumCols(); ++ c) {
            this->pg->at(r, c) = this->pg->at(r, c) * this->pg->at(r, c);
        }
    }
}

void LVTriangle::subtract(LVTriangle& squared) {
    for(int r = 0; r < this->pg->getNumRows(); ++ r) {
        for(int c = 0; c < this->pg->getNumCols(); ++ c) {
            this->pg->at(r, c) = this->pg->at(r, c) - squared.pg->at(r, c);
        }
    }
}

void LVTriangle::divide() {
    for(int r = 0; r < pg->getNumRows(); ++r) {
        for(int c = 0; c < pg->getNumCols(); ++c) {
            if(r <= c) {
                pg->at(r, c) = pg->at(r, c) / counts->at(r, c);
            }
        }
    }
}

LVTriangle::~LVTriangle() {
    delete pg;
    delete counts;
    pg = 0;
    counts = 0;
}

void LVTriangle::display() {
    (*pg).printOut();
}

#endif
