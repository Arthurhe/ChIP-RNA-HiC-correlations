#ifndef __CUMVAR_CPP__
#define __CUMVAR_CPP__

#include <cmath>
#include <cfloat>
#include <iostream>
#include <cassert>

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
    // Find max and min value
    for(int r = 0; r < this->pg->getNumRows(); ++ r) {
        for(int c = 0; c < this->pg->getNumCols(); ++ c) {
            float val = this->pg->at(r, c);
            if(r == 0 && c == 0) {
                maxVal = val;
                minVal = val;
            }
            else if(val > maxVal) {
                maxVal = val;
            }
            else if(val < minVal) {
                minVal = val;
            }
        }
    }

    // don't divide by zero and don't subtract by a negative number
    if(maxVal == 0) {
        maxVal = 1;
    }
    assert(minVal >= 0);

    // divide every value by maxVal
    for(int r = 0; r < this->pg->getNumRows(); ++ r) {
        for(int c = 0; c < this->pg->getNumCols(); ++ c) {
            this->pg->at(r, c) = (this->pg->at(r, c) - minVal) / maxVal;
        }
    }
}

int CumVar::getSize() {
    return pg->getNumCols();
}

float CumVar::getValue(int r, int c) {
    return pg->at(r, c);
}

CumVar::~CumVar() {
    delete pg;
    pg = 0;
}

void CumVar::display() {
    (*pg).printOut();
}

void CumVar::write(std::ofstream& fname) {
    for (int r = 0; r < pg->getNumRows(); ++r) {
        for (int c = 0; c < pg->getNumCols(); ++c) {
            float val = pg->at(r, c);
            if(std::isnan(val)) {
                val = -1;
            }
            fname << val << "\t";
        }
        fname << std::endl;
    }
    fname << std::endl;
}

#endif
