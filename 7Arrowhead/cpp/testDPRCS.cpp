#ifndef __TESTDPRCS_CPP__
#define __TESTDPRCS_CPP__

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <fstream>
#include <sstream>

#include "TwoD_Array.hpp"
#include "ThreeD_Array.hpp"
#include "RCS.hpp"

int main(int argc, char * argv[]) {
    if (argc != 2) {
        std::cerr << "Invalid number of arguments; expecting 1 for file name" << std::endl;
        exit(1);
    }

    std::ifstream input (argv[1], std::ios::in);

    // grab the grid size
    int size;
    std::string strSize;
    if (getline(input, strSize)) {
        std::stringstream stream(strSize);
        stream >> size; 
    }
    else {
        std::cerr << "Unable to open file '" << argv[1] << "'" << std::endl;
        exit(1);
    }

    TwoD_Array<int> grid(size, size);

    for (int i = 0; i < size; i++) {
        std::string row;
        if (getline(input, row)) {
            std::stringstream stream(row);
            for (int j = 0; j < size; j++) {
                stream >> grid.at(i, j);
            }
        }
        else {
            std::cerr << "Input file was incorrectly formatted per size" << std::endl;
            exit(1);
        }
    }

    grid.printOut();

    RCS rcs1(grid, 'r');
    rcs1.printOut(0);
    rcs1.printOut(1);

    std::cout << rcs1.query(0, 1, 2) << std::endl;;
    std::cout << rcs1.query(0, 1, 3) << std::endl;;
    std::cout << rcs1.query(0, 2, 3) << std::endl;;
    std::cout << rcs1.query(0, 0, 3) << std::endl;;
    std::cout << rcs1.query(0, 0, 2) << std::endl;;

    // RCS rcs2(grid, 'c');
    // rcs2.printOut(0);
    // rcs2.printOut(1);
    
    return 0;
}

#endif
