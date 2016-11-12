#ifndef __AGGREGATE_HPP__
#define __AGGREGATE_HPP__

#include <fstream>

#include "TwoD_Array.hpp"
#include "CumVar.hpp"
#include "LTriangle.hpp"

// Class to take in two triangle variance objects and return a cumulative mean
class Aggregate{
public:
    Aggregate(LTriangle& Ssum, LTriangle& Ssign, CumVar& Svar);
    ~Aggregate();

    void display();
    void write(std::ofstream& fname);

private:
    // precomputed grid
    TwoD_Array<float>* pg;
};

#endif
