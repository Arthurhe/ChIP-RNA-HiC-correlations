#ifndef __CUMVAR_CPP__
#define __CUMVAR_CPP__

#include <cmath>
#include <cfloat>
#include <iostream>

#include "CumVar.hpp"
#include "RCS.hpp"

/**
 * Note: CumVar takes in a 3D column sum matrix and outputs the mean of all
 *       values in the upper triangle.
 */

CumVar::CumVar(UVTriangle& upper, LVTriangle& lower) {
    // sanity checks?

    // initialize pg
    pg = new TwoD_Array<float>(upper.getSize(), upper.getSize());

    // Initialize everythhing to 0: fixes valgrind errors
    for(int i = 0; i < pg->getNumCols(); ++i) {
        for(int j = 0; j < pg->getNumRows(); ++j) {
            pg->at(j, i) = 0;
        }
    }
    
    for(int r = 0; r < pg->getNumRows(); ++r) {
        for(int c = 0; c < pg->getNumCols(); ++c) {
            float sum = upper.getSum(r, c) + lower.getSum(r, c);
            int count = upper.getCount(r, c) + lower.getCount(r, c);
            pg->at(r, c) = sum / count;
        }
    }
}

// public function, used to do (mean(X))^2
void CumVar::square() {
    for(int r = 0; r < this->pg->getNumRows(); ++ r) {
        for(int c = 0; c < this->pg->getNumCols(); ++ c) {
            this->pg->at(r, c) = this->pg->at(r, c) * this->pg->at(r, c);
        }
    }
}

// public function, used to do mean(X^2) - (mean(X))^2
void CumVar::subtract(CumVar& squared) {
    for(int r = 0; r < this->pg->getNumRows(); ++ r) {
        for(int c = 0; c < this->pg->getNumCols(); ++ c) {
            this->pg->at(r, c) = this->pg->at(r, c) - squared.pg->at(r, c);
        }
    }
}

// Scan through variances and get maximum value. Then divide every number by it
void CumVar::normalize() {
    // Find max value
    for(int r = 0; r < this->pg->getNumRows(); ++ r) {
        for(int c = 0; c < this->pg->getNumCols(); ++ c) {
            float val = this->pg->at(r, c);
            if(r == 0 && c == 0) {
                maxVal = val;
            }
            else if(val > maxVal) {
                maxVal = val;
            }
        }
    }

    // don't divide by zero
    if(maxVal == 0) {
        maxVal = 1;
    }

    // divide every value by maxVal
    for(int r = 0; r < this->pg->getNumRows(); ++ r) {
        for(int c = 0; c < this->pg->getNumCols(); ++ c) {
            this->pg->at(r, c) = this->pg->at(r, c) / maxVal;
        }
    }
}


CumVar::~CumVar() {
    delete pg;
    pg = 0;
}

void CumVar::display() {
    (*pg).printOut();
}

#endif
