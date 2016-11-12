#ifndef __CUMVAR_HPP__
#define __CUMVAR_HPP__

#include "TwoD_Array.hpp"
#include "UVTriangle.hpp"
#include "LVTriangle.hpp"

// Class to take in two triangle variance objects and return a cumulative mean
class CumVar{
public:
    // size - the triangle value matrix is NxN. size is the N
    // sfxn - scoring function matrix (from RCS.cpp)
    CumVar(UVTriangle& upper, LVTriangle& lower);
    ~CumVar();

    void display();

    void square();

    // Do the variance subtraction in place
    void subtract(CumVar& squared);

private:
    void divide();

    // precomputed grid
    TwoD_Array<float>* pg;
};

#endif
