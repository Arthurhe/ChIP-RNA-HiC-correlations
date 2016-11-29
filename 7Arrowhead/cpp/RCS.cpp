#ifndef __RCS_CPP__
#define __RCS_CPP__

#include <cmath>
#include <cstdlib>
#include <cfloat>
#include <cassert>

#include "TwoD_Array.hpp"
#include "ThreeD_Array.hpp"
#include "RCS.hpp"

RCS::RCS (TwoD_Array<float>& grid, char rc, char func) {
    // Applied function sanity check
    // s = sum, i = sign, v = variance
    assert(func == 's' || func == 'i');
    assert(rc == 'r' || rc == 'c');

    operation = func;
    rowcol = rc;
    pg = &grid;
}

RCS::~RCS() {
    pg = 0;
}


float RCS::sum(float a, float b) {
    return a + b;
}

float RCS::sign(float a, float b) {
    float bsign;

    if(b > 0) {
        bsign = 1;
    }
    else if (b < 0) {
        bsign = -1;
    }
    else {
        bsign = 0;
    }

    return a + bsign;
}

// Perform the query step here
// Note: Gets the sum [from, to] (inclusive)
float RCS::query (int rcNum, int from, int to) {
    float retval = 0;

    if(this->rowcol == 'c') {
        retval = pg->at(from, rcNum);
        for(int r = from+1; r <= to; ++r) {
            if(operation == 's') {
                retval = sum(retval, pg->at(r, rcNum));
            }
            else {
                retval = sign(retval, pg->at(r, rcNum));
            }
        }
    }
    else{
        retval = pg->at(rcNum, from);
        for(int c = from+1; c <= to; ++c) {
            if(operation == 's') {
                retval = sum(retval, pg->at(rcNum, c));
            }
            else {
                retval = sign(retval, pg->at(rcNum, c));
            }
        }
    }

    return retval;
}

// Function to print the computed matrix at a particular depth
void RCS::printOut(){
    pg->printOut();
}

#endif
