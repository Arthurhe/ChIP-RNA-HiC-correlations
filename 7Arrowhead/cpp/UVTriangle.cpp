#ifndef __UVTRIANGLE_CPP__
#define __UVTRIANGLE_CPP__

#include <cmath>
#include <cfloat>
#include <iostream>

#include "UVTriangle.hpp"
#include "RCS.hpp"

/**
 * Note: UVTriangle takes in a 3D column sum matrix and outputs the mean of all
 *       values in the upper triangle.
 */

UVTriangle::UVTriangle(int size, RCS& sfxnmatr) {
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

    compute(sfxnmatr);
    divide();
}

void UVTriangle::compute(RCS& sfxnmatr) {
    // compute
    for(int r = 0; r < pg->getNumRows(); ++r) {
        for(int c = 0; c < pg->getNumCols(); ++c) {
            // compute dimensions of triangle
            // note: top < bottom
            int top = r;
            int bottom = floor(((float)(r + c))/2);
            int right = 2*c - r;

            // recurrence
            if(r < c && right < pg->getNumCols()) {
                pg->at(r, c) = pg->at(r, c-1) + sfxnmatr.query(c, top, bottom);
                counts->at(r, c) = counts->at(r, c-1) + bottom - top + 1;
            }
            else if (r == c) {
                pg->at(r, c) = sfxnmatr.query(c, top, bottom);
                counts->at(r, c) = 1;
            }
            else {
                pg->at(r,c) = nanl("");
            }
        }
    }
}

void UVTriangle::square() {
    for(int r = 0; r < this->pg->getNumRows(); ++ r) {
        for(int c = 0; c < this->pg->getNumCols(); ++ c) {
            this->pg->at(r, c) = this->pg->at(r, c) * this->pg->at(r, c);
        }
    }
}

void UVTriangle::subtract(UVTriangle& squared) {
    for(int r = 0; r < this->pg->getNumRows(); ++ r) {
        for(int c = 0; c < this->pg->getNumCols(); ++ c) {
            this->pg->at(r, c) = this->pg->at(r, c) - squared.pg->at(r, c);
        }
    }
}

void UVTriangle::divide() {
    for(int r = 0; r < pg->getNumRows(); ++r) {
        for(int c = 0; c < pg->getNumCols(); ++c) {
            if(r <= c) {
                pg->at(r, c) = pg->at(r, c) / counts->at(r, c);
            }
        }
    }
}

UVTriangle::~UVTriangle() {
    delete pg;
    delete counts;
    pg = 0;
    counts = 0;
}

void UVTriangle::display() {
    (*pg).printOut();
}

#endif
