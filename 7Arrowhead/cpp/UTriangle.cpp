#ifndef __UTRIANGLE_CPP__
#define __UTRIANGLE_CPP__

#include <cmath>
#include <cfloat>

#include "UTriangle.hpp"
#include "RCS.hpp"
#include "UVTriangle.hpp"

UTriangle::UTriangle(int size, RCS& sfxnmatr) {
    // sanity checks?

    // initialize pg
    pg = new TwoD_Array<float>(size, size);

    // Initialize everythhing to 0: fixes valgrind errors
    for(int i = 0; i < pg->getNumCols(); ++i) {
        for(int j = 0; j < pg->getNumRows(); ++j) {
            pg->at(j, i) = 0;
        }
    }
    
    // compute
    for(int r = 0; r < pg->getNumRows(); ++r) {
        for(int c = 0; c < pg->getNumCols(); ++c) {
            // compute dimensions of triangle
            // note: top < bottom
            int top = r;
            int bottom = floor(((float)(r + c))/2);
            int right = 2*c - r;

            // recurrence
            // the right check is to not do calculations that cannot be done in
            // the lower triangle
            if(r < c && right < pg->getNumCols()) {
                pg->at(r, c) = pg->at(r, c-1) + sfxnmatr.query(c, top, bottom);
            }
            else if (r == c) {
                pg->at(r, c) = sfxnmatr.query(c, top, bottom);
            }
            else {
                pg->at(r,c) = nanl("");
            }
        }
    }
}

UTriangle::~UTriangle() {
    delete pg;
    pg = 0;
}

float UTriangle::getValue(int r, int c) {
    return this->pg->at(r, c);
}

void UTriangle::display() {
    (*pg).printOut();
}

void UTriangle::setCountPtr(UVTriangle& Upper) {
    upper = &Upper;
}

// Function to divide each number in the matrix by the number of elements in the
// triangle that it represents. This operation is performed in place.
void UTriangle::calcMean() {
    for(int r = 0; r < pg->getNumRows(); ++r) {
        for(int c = 0; c < pg->getNumCols(); ++c) {
            pg->at(r, c) = pg->at(r, c) / upper->getCount(r, c);
        }
    }
}

#endif
