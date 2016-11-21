#ifndef __RCS_HPP__
#define __RCS_HPP__

#include "ThreeD_Array.hpp"
#include "TwoD_Array.hpp"

class RCS {
public:
    // The constructor does the pre-computation
    RCS(TwoD_Array<float>& grid, char rc, char func);
    ~RCS();

    float query(int rcNum, int from, int to);

    void printOut(int d);

private:
    // The DP algorithm implementation
    void calculate(char rc, char func, TwoD_Array<float>& grid, float (*f)(float, float, float));

    // Scoring functions
    static float sum(float a, float b, float c);
    static float sign(float a, float b, float c);

    // pg is the precomputed grid
    TwoD_Array<float> * pg;
};

#endif
