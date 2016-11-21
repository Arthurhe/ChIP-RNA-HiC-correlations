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

    // Create the pre-computed grid
    if(rc == 'r' || rc == 'c') {
        pg = new TwoD_Array<float>(grid.getNumRows(), grid.getNumCols());
    }
    else {
        std::cerr << "Must select either rows or columns" << std::endl;
        exit(1);
    }

    // Initialize all values to 0
    for(int r = 0; r < pg->getNumCols(); ++r) {
        for(int c = 0; c < pg->getNumRows(); ++c) {
            pg->at(r, c) = 0;
        }
    }

    // Select function based on parameter
    float (*fxn)(float, float, float);
    if(func == 's') {
        fxn = &RCS::sum;
    }
    else if(func == 'i') {
        fxn = &RCS::sign;
    }

    // Compute the sums using DP
    if(rc == 'r') {
        RCS::calculate('r', func, grid, fxn);
    }
    else {
        RCS::calculate('c', func, grid, fxn);
    }
}

RCS::~RCS() {
    delete pg;
    pg = 0;
}

// Function to compute the DP recurrence
void RCS::calculate (char rc, char func, TwoD_Array<float>& grid, float(*f)(float, float, float)) {
    if(rc == 'r') {
        for(int r = 0; r < pg->getNumRows(); ++r) {
            for(int c = 0; c < pg->getNumCols(); ++c) {
                // Base Case
                if(r == c) {
                    if(func == 's') {
                        pg->at(r, c) = grid.at(r, c);
                    }
                    else if(func == 'i') {
                        pg->at(r, c) = (grid.at(r, c) >= 0) ? 1 : -1;
                    }
                }

                // Recurrence
                else if (r < c) {
                    pg->at(r, c) = (*f)(pg->at(r-1, c-1), grid.at(r, c));
                }

                // else set to maxint?
                else {
                    pg->at(r, c) = nanf("");
                }
            }
        }
    }
    else if(rc == 'c') {
        int lowerIndex;
        int lowerMod;
        float prevDPVal;
        float toSubtract;
        float toAdd;
        for(int r = 0; r < pg->getNumRows(); ++r) {
            for(int c = 0; c < pg->getNumCols(); ++c) {
                lowerIndex = floor(((float)(r+c))/2);
                lowerMod = (r+c) % 2;

                // base case
                if (r == c) {
                    if(func == 's') {
                        pg->at(r, c) = grid.at(r, c);
                    }
                    else if (func == 'i') {
                        pg->at(r, c) = (grid.at(r, c) >= 0) ? 1 : -1;
                    }
                }

                // recurrence
                else if (r < c && lowerMod != 0) {
                    prevDPVal = pg->at(r-1, c);
                    toSubtract = grid.at(r-1, c);
                    toAdd = 0;
                    pg->at(r, c) = (*f)(prevDPVal, toSubtract, toAdd);
                }
                else if (r < c && lowerMod == 0) {
                    prevDPVal = pg->at(r-1, c);
                    toSubtract = grid.at(r-1, c);
                    toAdd = grid.at(lowerIndex, c);
                    pg->at(r, c) = (*f)(prevDPVal, toSubtract, toAdd);
                }

                // else set to maxint?
                else {
                    pg->at(r, c) = nanf("");
                }
            }
        }
    }
}

float RCS::sum(float a, float b, float c) {
    return a - b + c;
}

float RCS::sign(float a, float b, float c) {
    float csign = (c >= 0) ? 1 : -1;

    return a - b + csign;
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
