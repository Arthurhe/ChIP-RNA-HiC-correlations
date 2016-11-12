#ifndef __LVTRIANGLE_HPP__
#define __LVTRIANGLE_HPP__

#include "TwoD_Array.hpp"
#include "RCS.hpp"

class LVTriangle{
public:
    // size - the triangle value matrix is NxN. size is the N
    // sfxn - scoring function matrix (from RCS.cpp)
    LVTriangle(int size, RCS& rScores, RCS& cScores);
    ~LVTriangle();

    void display();

    void square();

    // Getters for CumVar
    float getSum(int r, int c);
    int getCount(int r, int c);
    int getSize();
    
    // Do the variance subtraction in place
    void subtract(LVTriangle& squared);

private:
    void compute(RCS& rScores, RCS& cScores);
    void divide();

    // precomputed grid
    TwoD_Array<float>* pg;
    TwoD_Array<int>* counts;
};

#endif
