#ifndef __TESTUVTRIANGLE_CPP__
#define __TESTUVTRIANGLE_CPP__

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <fstream>
#include <sstream>

#include "TwoD_Array.hpp"
#include "RCS.hpp"
#include "UVTriangle.hpp"

// first input is the normal matrix, second input is the squared matrix

int main(int argc, char * argv[]) {
    if (argc != 3) {
        std::cerr << "Invalid number of arguments; expecting 2 for file name" << std::endl;
        exit(1);
    }

    std::ifstream input1 (argv[1], std::ios::in);

    // grab the grid size
    int size1;
    std::string strSize1;
    if (getline(input1, strSize1)) {
        std::stringstream stream(strSize1);
        stream >> size1; 
    }
    else {
        std::cerr << "Unable to open file '" << argv[1] << "'" << std::endl;
        exit(1);
    }

    TwoD_Array<float> grid1(size1, size1);

    for (int i = 0; i < size1; i++) {
        std::string row;
        if (getline(input1, row)) {
            std::stringstream stream(row);
            for (int j = 0; j < size1; j++) {
                stream >> grid1.at(i, j);
            }
        }
        else {
            std::cerr << "input1 file was incorrectly formatted per size" << std::endl;
            exit(1);
        }
    }

    std::ifstream input2 (argv[2], std::ios::in);

    // grab the grid size
    int size2;
    std::string strSize2;
    if (getline(input2, strSize2)) {
        std::stringstream stream(strSize2);
        stream >> size2; 
    }
    else {
        std::cerr << "Unable to open file '" << argv[2] << "'" << std::endl;
        exit(1);
    }

    TwoD_Array<float> grid2(size2, size2);

    for (int i = 0; i < size2; i++) {
        std::string row;
        if (getline(input2, row)) {
            std::stringstream stream(row);
            for (int j = 0; j < size2; j++) {
                stream >> grid2.at(i, j);
            }
        }
        else {
            std::cerr << "input2 file was incorrectly formatted per size" << std::endl;
            exit(1);
        }
    }

    grid1.printOut();
    grid2.printOut();

    RCS cScores(grid1, 'c', 's');
    RCS c2Scores(grid2, 'c', 's');
    UVTriangle sum = UVTriangle(size1, cScores);
    UVTriangle sum2 = UVTriangle(size2, c2Scores);
    sum.square();
    sum2.subtract(sum);

    sum2.display();
    
    return 0;
}

#endif
