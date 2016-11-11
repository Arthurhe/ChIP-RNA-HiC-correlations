#ifndef __LTRIANGLE_HPP__
#define __LTRIANGLE_HPP__

#include "TwoD_Array.hpp"
#include "RCS.hpp"

class LTriangle{
public:
    // size - the triangle value matrix is NxN. size is the N
    // sfxn - scoring function matrix (from RCS.cpp)
    LTriangle(int size, RCS& rScores, RCS& cScores);
    ~LTriangle();

    void display();

private:
    // precomputed grid
    TwoD_Array<float>* pg;
};

#endif
