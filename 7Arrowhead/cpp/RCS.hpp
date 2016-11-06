#ifndef __DPRCS_HPP__
#define __DPRCS_HPP__

#include "ThreeD_Array.hpp"
#include "TwoD_Array.hpp"

class RCS {
public:
    // The constructor does the pre-computation
    RCS(TwoD_Array<int>& grid, char rc);

    int query(int rcNum, int from, int to);

    void printOut(int d);

private:
    void calculate(int d, char rc, TwoD_Array<int>& grid);

    // pg is the precomputed grid
    ThreeD_Array<int> * pg;
};

#endif
