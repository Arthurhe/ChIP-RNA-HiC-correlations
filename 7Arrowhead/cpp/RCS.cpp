// CSE 101 Winter 2016, PA 3
//
// Name: TODO put both partners' info if applicable
// PID: TODO
// Sources of Help: TODO
// Due: February 19th, 2016 at 11:59 PM

#ifndef __GRID_SUM_CPP__
#define __GRID_SUM_CPP__

#include "ThreeD_Array.hpp"
#include "RCS.hpp"

// Perform the precomputation step here
RCS::RCS (ThreeD_Array<int>& grid) {
    // TODO
    // Create the pre-computed grid
    // TODO: Fix the memory leak?
    pg = new ThreeD_Array<int>(1, grid.getNumRows(), grid.getNumCols());

    // TODO: remove when done
    // Initialize all values to 0 for debuging purposes
    for(int i = 0; i < pg->getNumCols(); ++i) {
        for(int j = 0; j < pg->getNumRows(); ++j) {
            pg->at(0, j, i) = 0;
        }
    }

    pg->at(0, 0, 0) = grid.at(0, 0, 0);

    // Compute the first row of elements in pre-computed grid pg O(N)
    for(int i = 1; i < grid.getNumCols(); ++i) {
        pg->at(0, 0, i) = pg->at(0, 0, i-1) + grid.at(0, 0, i);
    }

    // Compute the first column of elements in the pre-computed grid pg O(N)
    for(int i = 1; i < grid.getNumRows(); ++i) {
        pg->at(0, i, 0) = pg->at(0, i-1, 0) + grid.at(0, i, 0);
    }

    // Compute the rest of the grid O(N^2)
    for(int r = 1; r < pg->getNumRows(); ++r) {
        for(int c = 1; c < pg->getNumCols(); ++c) {
            pg->at(0, r, c) = pg->at(0, r-1, c) + pg->at(0, r, c-1) + grid.at(0, r, c) - pg->at(0, r-1, c-1);
        }
    }
}

// Perform the query step here
int RCS::query (int d, int x1, int y1, int x2, int y2) {
    // TODO

    // Case where one corner is (0, 0)
    if(x1 == 0 && y1 == 0) {
        return pg->at(0, x2, y2);
    }
    else if(x1 == 0) {
        return pg->at(0, x2, y2) - pg->at(0, x2, y1-1);
    }
    else if(y1 == 0) {
        return pg->at(0, x2, y2) - pg->at(0, x1-1, y2);
    }
    else if(x2 < pg->getNumRows() && y2 < pg->getNumCols()) {
        return pg->at(0, x2, y2) - pg->at(0, x1-1, y2) - pg->at(0, x2, y1-1) + pg->at(0, x1-1, y1-1);
    }
    // TODO: Check for out of bounds in the greater direction?
    return -1;
}

#endif
