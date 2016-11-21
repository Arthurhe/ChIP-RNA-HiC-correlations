#ifndef __AGGREGATE_CPP__
#define __AGGREGATE_CPP__

#include <cmath>
#include <iostream>
#include <fstream>

#include "Aggregate.hpp"
#include "CumVar.hpp"
#include "UTriangle.hpp"
#include "LTriangle.hpp"

static const float var_t = 0.2;
static const float usign_t1 = -0.5;
static const float lsign_t1 = 0.5;

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

// Copy constructor, probably not necessary
Aggregate::Aggregate(const Aggregate& orig) {
    pg = new TwoD_Array<float>(orig.pg->getNumRows(), orig.pg->getNumCols());
    for(int r = 0; r < orig.pg->getNumRows(); ++r) {
        for(int c = 0; c < orig.pg->getNumCols(); ++c) {
            pg->at(r, c) = orig.pg->at(r, c);
        }
    }
}

Aggregate::~Aggregate() {
    delete pg;
    pg = 0;
}

// Filter by Svar and Mean(Sign(Upper/Lower))
void Aggregate::filter1(UTriangle& meanUSign, LTriangle& meanLSign, CumVar& Svar) {
    for(int r = 0; r < pg->getNumRows(); ++r) {
        for(int c = 0; c < pg->getNumCols(); ++c) {
            if(Svar.getValue(r, c) < var_t) {
                pg->at(r, c) = 0;
            } else if(meanUSign.getValue(r, c) < usign_t1) {
                pg->at(r, c) = 0;
            } else if(meanLSign.getValue(r, c) > lsign_t1) {
                pg->at(r, c) = 0;
            }
            else if(pg->at(r, c) == -1) {
                pg->at(r, c) = 0;
            }
        }
    }
}

// TODO: not used yet. just try to get one level threshold working first
void Aggregate::filter2(UTriangle& meanUSign, LTriangle& meanLSign, CumVar& Svar) {
    for(int r = 0; r < pg->getNumRows(); ++r) {
        for(int c = 0; c < pg->getNumCols(); ++c) {
            if(Svar.getValue(r, c) < var_t) {
                pg->at(r, c) = 0;
            } else if(meanUSign.getValue(r, c) < usign_t1) {
                pg->at(r, c) = 0;
            } else if(meanLSign.getValue(r, c) > lsign_t1) {
                pg->at(r, c) = 0;
            }
        }
    }
}

// Filter by Svar and Mean(Sign(Upper/Lower))
void Aggregate::toBinary() {
    for(int r = 0; r < pg->getNumRows(); ++r) {
        for(int c = 0; c < pg->getNumCols(); ++c) {
            if(pg->at(r, c) > 0) {
                pg->at(r,c ) = 1;
            }
        }
    }
}

void Aggregate::display() {
    (*pg).printOut();
}

void Aggregate::write(std::ofstream& fname) {
    for (int r = 0; r < pg->getNumRows(); ++r) {
        for (int c = 0; c < pg->getNumCols(); ++c) {
            float val = pg->at(r, c);
            if(std::isnan(val)) {
                val = -1;
            }
            fname << val << "\t";
        }
        fname << std::endl;
    }
    fname << std::endl;
}


#endif
