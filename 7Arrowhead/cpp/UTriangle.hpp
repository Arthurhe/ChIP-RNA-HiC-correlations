#ifndef __UTRIANGLE_HPP__
#define __UTRIANGLE_HPP__

#include "TwoD_Array.hpp"
#include "RCS.hpp"

class UTriangle{
public:
    // size - the triangle value matrix is NxN. size is the N
    // sfxn - scoring function matrix (from RCS.cpp)
    UTriangle(int size, RCS& sfxnmatr);
    ~UTriangle();

    void display();

private:
    // precomputed grid
    TwoD_Array<float>* pg;
};

#endif
