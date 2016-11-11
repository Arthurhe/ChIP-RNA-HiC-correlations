#ifndef __UTRIANGLE_CPP__
#define __UTRIANGLE_CPP__

#include <cmath>
#include <cfloat>

#include "UTriangle.hpp"
#include "RCS.hpp"

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
            int bottom = floor(((float) r + (float) c)/2);

            // recurrence
            if(r < c) {
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

void UTriangle::display() {
    (*pg).printOut();
}

#endif
