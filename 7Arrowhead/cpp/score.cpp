#ifndef __SCORE_CPP__
#define __SCORE_CPP__

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <fstream>
#include <sstream>

#include "TwoD_Array.hpp"
#include "RCS.hpp"
#include "UVTriangle.hpp"
#include "LVTriangle.hpp"
#include "UTriangle.hpp"
#include "LTriangle.hpp"

// first input is the normal matrix, second input is the squared matrix

int main(int argc, char * argv[]) {
    if (argc != 3) {
        std::cerr << "Invalid number of arguments; expecting 2 for file name" << std::endl;
        exit(1);
    }

    /**
     * ====================================================================
     * Read in the normal matrix file
     * ====================================================================
     */
    std::cout << "Reading in the normal matrix file" << std::endl;
    std::ifstream norm_in (argv[1], std::ios::in);

    // grab the grid size
    int size1;
    std::string strSize1;
    if (getline(norm_in, strSize1)) {
        std::stringstream stream(strSize1);
        stream >> size1; 
    }
    else {
        std::cerr << "Unable to open file '" << argv[1] << "'" << std::endl;
        exit(1);
    }

    TwoD_Array<float> norm_matr(size1, size1);

    for (int i = 0; i < size1; i++) {
        std::string row;
        if (getline(norm_in, row)) {
            std::stringstream stream(row);
            for (int j = 0; j < size1; j++) {
                stream >> norm_matr.at(i, j);
            }
        }
        else {
            std::cerr << "norm_in file was incorrectly formatted per size" << std::endl;
            exit(1);
        }
    }

    /**
     * ====================================================================
     * Read in the squared matrix file
     * ====================================================================
     */
    std::cout << "Reading in the squared matrix file" << std::endl;
    std::ifstream squared_in (argv[2], std::ios::in);

    // grab the grid size
    int size2;
    std::string strSize2;
    if (getline(squared_in, strSize2)) {
        std::stringstream stream(strSize2);
        stream >> size2; 
    }
    else {
        std::cerr << "Unable to open file '" << argv[2] << "'" << std::endl;
        exit(1);
    }

    TwoD_Array<float> squared_matr(size2, size2);

    for (int i = 0; i < size2; i++) {
        std::string row;
        if (getline(squared_in, row)) {
            std::stringstream stream(row);
            for (int j = 0; j < size2; j++) {
                stream >> squared_matr.at(i, j);
            }
        }
        else {
            std::cerr << "squared_in file was incorrectly formatted per size" << std::endl;
            exit(1);
        }
    }

    /**
     * ====================================================================
     * Compute U_sum
     * ====================================================================
     */
    std::cout << "Computing Usums" << std::endl;
    RCS normColSums(norm_matr, 'c', 's');
    UTriangle Usums(size1, normColSums);
    Usums.display();

    /**
     * ====================================================================
     * Compute L_sum
     * ====================================================================
     */
    std::cout << "Computing Lsums" << std::endl;
    RCS normRowSums(norm_matr, 'r', 's');
    LTriangle Lsums(size1, normRowSums, normColSums);
    Lsums.display();

    /**
     * ====================================================================
     * Compute U_sign
     * ====================================================================
     */
    std::cout << "Computing Usigns" << std::endl;
    RCS normColSigns(norm_matr, 'c', 'i');
    UTriangle Usigns(size1, normColSigns);
    Usigns.display();

    /**
     * ====================================================================
     * Compute L_sign
     * ====================================================================
     */
    std::cout << "Computing Lsigns" << std::endl;
    RCS normRowSigns(norm_matr, 'r', 'i');
    LTriangle Lsigns(size1, normRowSigns, normColSigns);
    Lsigns.display();
    
    /**
     * ====================================================================
     * Compute U_var
     * ====================================================================
     */
    std::cout << "Computing Uvars" << std::endl;
    RCS squaredColSums(squared_matr, 'c', 's');
    UVTriangle utmp = UVTriangle(size1, normColSums);
    UVTriangle Uvar = UVTriangle(size2, squaredColSums);
    utmp.square();
    Uvar.subtract(utmp);

    Uvar.display();

    /**
     * ====================================================================
     * Compute L_var
     * ====================================================================
     */
    std::cout << "Computing Lvars" << std::endl;
    RCS squaredRowSums(squared_matr, 'c', 's');
    LVTriangle ltmp = LVTriangle(size1, normRowSums, normColSums);
    LVTriangle Lvar = LVTriangle(size2, squaredRowSums, squaredColSums);
    ltmp.square();
    Lvar.subtract(ltmp);

    Lvar.display();

    std::cout << "Done!" << std::endl;
    
    return 0;
}

#endif
