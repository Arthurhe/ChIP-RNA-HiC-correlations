#ifndef __SCORE_CPP__
#define __SCORE_CPP__

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <fstream>
#include <sstream>
#include <ctime>

#include "TwoD_Array.hpp"
#include "RCS.hpp"
#include "UVTriangle.hpp"
#include "LVTriangle.hpp"
#include "UTriangle.hpp"
#include "LTriangle.hpp"
#include "CumVar.hpp"
#include "Aggregate.hpp"

// first input is the normal matrix, second input is the squared matrix

int main(int argc, char * argv[]) {
    // Get start time
    std::time_t start = std::time(0);

    if (argc != 3) {
        std::cerr << "Invalid number of arguments; expecting 2 for file name" << std::endl;
        exit(1);
    }

    /**
     * ====================================================================
     * Read in the normal matrix file
     * ====================================================================
     */
    std::time_t now = std::time(0);
    std::cout << "[Elapsed: " << now - start << "s]\tReading in the normal matrix file" << std::endl;
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
    now = std::time(0);
    std::cout << "[Elapsed: " << now - start << "s]\tReading in the squared matrix file" << std::endl;
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
    now = std::time(0);
    std::cout << "[Elapsed: " << now - start << "s]\tComputing Usums" << std::endl;
    RCS normColSums(norm_matr, 'c', 's');
    UTriangle Usums(size1, normColSums);
    // Usums.display();

    /**
     * ====================================================================
     * Compute L_sum
     * ====================================================================
     */
    now = std::time(0);
    std::cout << "[Elapsed: " << now - start << "s]\tComputing Lsums" << std::endl;
    RCS normRowSums(norm_matr, 'r', 's');
    LTriangle Lsums(size1, normRowSums, normColSums);
    // Lsums.display();

    /**
     * ====================================================================
     * Compute U_sign
     * ====================================================================
     */
    now = std::time(0);
    std::cout << "[Elapsed: " << now - start << "s]\tComputing Usigns" << std::endl;
    RCS normColSigns(norm_matr, 'c', 'i');
    UTriangle Usigns(size1, normColSigns);
    // Usigns.display();

    /**
     * ====================================================================
     * Compute L_sign
     * ====================================================================
     */
    now = std::time(0);
    std::cout << "[Elapsed: " << now - start << "s]\tComputing Lsigns" << std::endl;
    RCS normRowSigns(norm_matr, 'r', 'i');
    LTriangle Lsigns(size1, normRowSigns, normColSigns);
    // Lsigns.display();

    /**
     * ====================================================================
     * Compute S_sum = L_sum - U_sum
     * ====================================================================
     */
    now = std::time(0);
    std::cout << "[Elapsed: " << now - start << "s]\tComputing S_sum = L_sum - U_sum" << std::endl;
    LTriangle Ssum = LTriangle(Lsums);
    Ssum.subtract(Usums);
    Ssum.normalize();

    // TODO: replace with a copy constructor thing
    // Lsums.subtract(Usums);
    // LTriangle* Ssum = &Lsums;
    // Ssum->normalize();

    std::ofstream sumfile;
    sumfile.open("SsumMatr.txt");
    Ssum.write(sumfile);
    sumfile.close();

    //Ssum->display();

    /**
     * ====================================================================
     * Compute S_sign
     * ====================================================================
     */
    now = std::time(0);
    std::cout << "[Elapsed: " << now - start << "s]\tComputing S_sign = L_sign - U_sign" << std::endl;
    LTriangle Ssign = LTriangle(Lsigns);
    Ssign.subtract(Usigns);
    Ssign.normalize();

    // TODO: replace with a copy constructor thing
    // Lsigns.subtract(Usigns);
    // LTriangle* Ssign = &Lsigns;
    // Ssign->normalize();

    std::ofstream signfile;
    signfile.open("SsignMatr.txt");
    Ssign.write(signfile);
    signfile.close();


    //Ssign->display();

    /**
     * ====================================================================
     * Compute S_var
     * ====================================================================
     */
    now = std::time(0);
    std::cout << "[Elapsed: " << now - start << "s]\tComputing S_var" << std::endl;
    UVTriangle untmp = UVTriangle(size1, normColSums);
    LVTriangle lntmp = LVTriangle(size1, normRowSums, normColSums);
    CumVar cvtmp = CumVar(untmp, lntmp);

    RCS squaredColSums(squared_matr, 'c', 's');
    RCS squaredRowSums(squared_matr, 'r', 's');
    UVTriangle ustmp = UVTriangle(size2, squaredColSums);
    LVTriangle lstmp = LVTriangle(size2, squaredRowSums, squaredColSums);
    CumVar Svar = CumVar(ustmp, lstmp);

    cvtmp.square();
    Svar.subtract(cvtmp);
    Svar.normalize();

    std::ofstream varfile;
    varfile.open("SvarMatr.txt");
    Svar.write(varfile);
    varfile.close();

    //Svar.display();

    /**
     * ====================================================================
     * Aggregate Scores to make S'
     * ====================================================================
     */
    now = std::time(0);
    std::cout << "[Elapsed: " << now - start << "s]\tAggregating scores to make S'" << std::endl;
    Aggregate Stot = Aggregate(Ssum, Ssign, Svar);
    //Stot.display();
    std::ofstream outfile;
    outfile.open("CumScoreMatr.txt");
    Stot.write(outfile);
    outfile.close();

    /**
     * ====================================================================
     * Filter Scores to make S
     * ====================================================================
     */
    now = std::time(0);
    std::cout << "[Elapsed: " << now - start << "s]\tFiltering scores to make S" << std::endl;
    Usigns.setCountPtr(untmp);
    Lsigns.setCountPtr(lntmp);
    Usigns.calcMean();
    Lsigns.calcMean();
    Stot.filter1(Usigns, Lsigns, Svar, 0.2, 0.5);
    std::ofstream filtered;
    filtered.open("FilteredScoreMatr.txt");
    Stot.write(filtered);
    filtered.close();

    /**
     * ====================================================================
     * Create binary matrix for connected component analysis
     * ====================================================================
     */
    now = std::time(0);
    std::cout << "[Elapsed: " << now - start << "s]\tCreating binary matrix" << std::endl;
    Aggregate binaryMatr = Aggregate(Stot);
    binaryMatr.toBinary();
    std::ofstream binary;
    binary.open("BinaryMatrix.txt");
    binaryMatr.write(binary);
    binary.close();

    now = std::time(0);
    std::cout << "[Elapsed: " << now - start << "s]\tDone!" << std::endl;
    
    return 0;
}

#endif
