#ifndef __RCS_CPP__
#define __RCS_CPP__

#include <cstdlib>
#include <cfloat>
#include <cassert>

#include "TwoD_Array.hpp"
#include "ThreeD_Array.hpp"
#include "RCS.hpp"

RCS::RCS (TwoD_Array<float>& grid, char rc, char func) {
    // Applied function sanity check
    // s = sum, i = sign, v = variance
    assert(func == 's' || func == 'i' || func == 'v');

    // Create the pre-computed grid
    if(rc == 'r') {
        pg = new ThreeD_Array<float>(grid.getNumRows(), grid.getNumRows(), grid.getNumCols());
    }
    else if (rc == 'c') {
        pg = new ThreeD_Array<float>(grid.getNumCols(), grid.getNumRows(), grid.getNumCols());
    }
    else {
        std::cerr << "Must select either rows or columns" << std::endl;
        exit(1);
    }

    // Initialize all values to 0
    for(int d = 0; d < pg->getDepth(); ++d) {
        for(int i = 0; i < pg->getNumCols(); ++i) {
            for(int j = 0; j < pg->getNumRows(); ++j) {
                pg->at(d, j, i) = 0;
            }
        }
    }

    // Select function based on parameter
    float (*fxn)(float, float);
    if(func == 's') {
        fxn = &RCS::sum;
    }
    else if(func == 'i') {
    }
    else if(func == 'v') {
    }


    // Compute the sums using DP
    if(rc == 'r') {
        for(int d = 0; d < grid.getNumRows(); ++d) {
            RCS::calculate(d, 'r', grid, fxn);
        }
    }
    else {
        for(int d = 0; d < grid.getNumCols(); ++d) {
            RCS::calculate(d, 'c', grid, fxn);
        }
    }
}

RCS::~RCS() {
    delete pg;
    pg = 0;
}

// Function to compute the DP recurrence per row/column of the given grid
void RCS::calculate (int d, char rc, TwoD_Array<float>& grid, float(*f)(float, float)) {
    for(int r = 0; r < pg->getNumRows(); ++r) {
        for(int c = 0; c < pg->getNumCols(); ++c) {

            // Base Case
            if(r == 0 && rc == 'r') {
                pg->at(d, r, c) = grid.at(d, c);
            }
            else if (r == 0 && rc == 'c') {
                pg->at(d, r, c) = grid.at(c, d);
            }

            // Recurrence
            else if (r > 0 && c > 0 && c >= r && rc == 'r') {
                pg->at(d, r, c) = (*f)(pg->at(d, r-1, c-1), grid.at(d, c));
            }
            else if (r > 0 && c > 0 && c >= r && rc == 'c') {
                pg->at(d, r, c) = (*f)(pg->at(d, r-1, c-1), grid.at(c, d));
            }

            // else set to maxint?
            else {
                pg->at(d, r, c) = FLT_MAX;
            }
        }
    }
}

float RCS::sum(float a, float b) {
    return a + b;
}

float RCS::sign(float a, float b) {
    float asign;
    float bsign;
    (a >= 0) ? asign = 1 : asign = -1;
    (b >= 0) ? bsign = 1 : bsign = -1;

    return asign + bsign;
}

// Perform the query step here
// Note: Gets the sum [from, to] (inclusive)
float RCS::query (int rcNum, int from, int to) {
    // Do conversion calculation to query the pg
    assert(to >= from);
    int range = to - from;

    return pg->at(rcNum, range, to);
}

// Function to print the computed matrix at a particular depth
void RCS::printOut(int d){
    pg->printOut(d);
}

#endif
