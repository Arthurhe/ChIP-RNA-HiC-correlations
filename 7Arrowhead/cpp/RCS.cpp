#ifndef __DPRCS_CPP__
#define __DPRCS_CPP__

#include <cstdlib>
#include <climits>
#include <cassert>

#include "TwoD_Array.hpp"
#include "ThreeD_Array.hpp"
#include "RCS.hpp"

// Perform the precomputation step here
RCS::RCS (TwoD_Array<int>& grid, char rc) {
    // Create the pre-computed grid
    if(rc == 'r') {
        pg = new ThreeD_Array<int>(grid.getNumRows(), grid.getNumRows(), grid.getNumCols());
    }
    else if (rc == 'c') {
        pg = new ThreeD_Array<int>(grid.getNumCols(), grid.getNumRows(), grid.getNumCols());
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

    // Compute the sums using DP
    if(rc == 'r') {
        for(int d = 0; d < grid.getNumRows(); ++d) {
            RCS::calculate(d, 'r', grid);
        }
    }
    else {
        for(int d = 0; d < grid.getNumCols(); ++d) {
            RCS::calculate(d, 'c', grid);
        }
    }
}

// Function to compute the DP recurrence per row/column of the given grid
void RCS::calculate (int d, char rc, TwoD_Array<int>& grid) {
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
                pg->at(d, r, c) = pg->at(d, r-1, c-1) + grid.at(d, c);
            }
            else if (r > 0 && c > 0 && c >= r && rc == 'c') {
                pg->at(d, r, c) = pg->at(d, r-1, c-1) + grid.at(c, d);
            }

            // else set to maxint?
            else {
                pg->at(d, r, c) = INT_MAX;
            }
        }
    }
}

// Perform the query step here
// Note: Gets the sum [from, to] (inclusive)
int RCS::query (int rcNum, int from, int to) {
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
