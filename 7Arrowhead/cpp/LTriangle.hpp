#ifndef __LTRIANGLE_HPP__
#define __LTRIANGLE_HPP__

#include <fstream>

#include "TwoD_Array.hpp"
#include "RCS.hpp"
#include "UTriangle.hpp"
#include "LVTriangle.hpp"

class LTriangle{
public:
    // size - the triangle value matrix is NxN. size is the N
    // sfxn - scoring function matrix (from RCS.cpp)
    LTriangle(int size, RCS& rScores, RCS& cScores);
    LTriangle(const LTriangle& orig);
    ~LTriangle();

    // Used to do L - U in place
    void subtract(UTriangle& upper);

    // Divide by largest value
    void normalize();

    void display();

    float getValue(int r, int c);
    void write(std::ofstream& fname);

    // Sets the pointers for lower so I can get the count
    void setCountPtr(LVTriangle& Lower);
    void calcMean();

private:
    // precomputed grid
    TwoD_Array<float>* pg;

    float maxVal;

    // Pointers to the UVTriangle and LVTriangle objects: They have matrices
    // containing the number of elements in each triangle. Used for computing
    // mean(sign(Upper/Lower)) in the filtering step.
    LVTriangle* lower;
};

#endif
