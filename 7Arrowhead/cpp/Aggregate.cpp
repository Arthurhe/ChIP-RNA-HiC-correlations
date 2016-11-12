#ifndef __AGGREGATE_CPP__
#define __AGGREGATE_CPP__

#include <cmath>
#include <cfloat>
#include <iostream>
#include <fstream>

#include "Aggregate.hpp"
#include "CumVar.hpp"
#include "LTriangle.hpp"

Aggregate::Aggregate(LTriangle& Ssum, LTriangle& Ssign, CumVar& Svar) {
    // TODO: Check that all the input matrices are the same size? might not have to

    // initialize pg
    pg = new TwoD_Array<float>(Svar.getSize(), Svar.getSize());

    // Initialize everythhing to 0: fixes valgrind errors
    for(int i = 0; i < pg->getNumCols(); ++i) {
        for(int j = 0; j < pg->getNumRows(); ++j) {
            pg->at(j, i) = 0;
        }
    }

    // Sum the values from all matrices
    for(int r = 0; r < pg->getNumRows(); ++r) {
        for(int c = 0; c < pg->getNumCols(); ++c) {
            pg->at(r, c) = Ssum.getValue(r, c) +
                           Ssign.getValue(r, c) +
                           Svar.getValue(r, c);
        }
    }
}

Aggregate::~Aggregate() {
    delete pg;
    pg = 0;
}

void Aggregate::display() {
    (*pg).printOut();
}

void Aggregate::write(std::ofstream& fname) {
    for (int r = 0; r < pg->getNumRows(); ++r) {
        for (int c = 0; c < pg->getNumCols(); ++c) {
            fname << pg->at(r,c) << "\t";
        }
        fname << std::endl;
    }
    fname << std::endl;
}

#endif
