#ifndef __AGGREGATE_HPP__
#define __AGGREGATE_HPP__

#include <fstream>

#include "TwoD_Array.hpp"
#include "CumVar.hpp"
#include "LTriangle.hpp"
#include "UTriangle.hpp"

// Class to take in two triangle variance objects and return a cumulative mean
class Aggregate{
public:
    Aggregate(LTriangle& Ssum, LTriangle& Ssign, CumVar& Svar);
    Aggregate(const Aggregate& original);
    ~Aggregate();

    void filter1(UTriangle& meanUSign, LTriangle& meanLSign, CumVar& Svar, float var_t, float sign_t);
    // void filter2(UTriangle& meanUSign, LTriangle& meanLSign, CumVar& Svar);
    void toBinary();

    void display();
    void write(std::ofstream& fname);

private:
    // precomputed grid
    TwoD_Array<float>* pg;
    
};

#endif
