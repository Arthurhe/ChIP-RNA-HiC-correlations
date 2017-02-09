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
    // TODO: find out a clever way to return both the sums and the variances
    // Right now, i just want cumulative variances, so don't divide
    // divide();
}

// Private function, only used by ctor
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

// public function, used to do (mean(X))^2
void UVTriangle::square() {
    for(int r = 0; r < this->pg->getNumRows(); ++ r) {
        for(int c = 0; c < this->pg->getNumCols(); ++ c) {
            this->pg->at(r, c) = this->pg->at(r, c) * this->pg->at(r, c);
        }
    }
}

// public function, used to do mean(X^2) - (mean(X))^2
void UVTriangle::subtract(UVTriangle& squared) {
    for(int r = 0; r < this->pg->getNumRows(); ++ r) {
        for(int c = 0; c < this->pg->getNumCols(); ++ c) {
            this->pg->at(r, c) = this->pg->at(r, c) - squared.pg->at(r, c);
        }
    }
}

// private function, used to compute the mean of the triangle, used by ctor
void UVTriangle::divide() {
    for(int r = 0; r < pg->getNumRows(); ++r) {
        for(int c = 0; c < pg->getNumCols(); ++c) {
            if(r <= c) {
                pg->at(r, c) = pg->at(r, c) / counts->at(r, c);
            }
        }
    }
}

// public function
float UVTriangle::getSum(int r, int c) {
    return pg->at(r, c);
}

// public function
int UVTriangle::getCount(int r, int c) {
    return counts->at(r, c);
}

// public function
int UVTriangle::getSize() {
    return pg->getNumRows();
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
