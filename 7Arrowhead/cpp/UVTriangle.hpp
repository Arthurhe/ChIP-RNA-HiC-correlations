#ifndef __UVTRIANGLE_HPP__
#define __UVTRIANGLE_HPP__

#include "TwoD_Array.hpp"
#include "RCS.hpp"

class UVTriangle{
public:
    // size - the triangle value matrix is NxN. size is the N
    // sfxn - scoring function matrix (from RCS.cpp)
    UVTriangle(int size, RCS& sfxnmatr);
    ~UVTriangle();

    void display();

    void square();
    
    // Do the variance subtraction in place
    void subtract(UVTriangle& squared);

private:
    void compute(RCS& sfxnmatr);
    void divide();

    // precomputed grid
    TwoD_Array<float>* pg;
    TwoD_Array<int>* counts;
};

#endif
