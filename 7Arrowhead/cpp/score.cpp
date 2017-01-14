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

    Aggregate sign01 = Aggregate(Stot);
    Aggregate sign02 = Aggregate(Stot);
    Aggregate sign03 = Aggregate(Stot);
    Aggregate sign04 = Aggregate(Stot);
    Aggregate sign05 = Aggregate(Stot);
    Aggregate sign06 = Aggregate(Stot);
    Aggregate sign07 = Aggregate(Stot);
    Aggregate sign08 = Aggregate(Stot);
    Aggregate sign09 = Aggregate(Stot);

    sign01.filter1(Usigns, Lsigns, Svar, 0.2, 0.1);
    sign02.filter1(Usigns, Lsigns, Svar, 0.2, 0.2);
    sign03.filter1(Usigns, Lsigns, Svar, 0.2, 0.3);
    sign04.filter1(Usigns, Lsigns, Svar, 0.2, 0.4);
    sign05.filter1(Usigns, Lsigns, Svar, 0.2, 0.5);
    sign06.filter1(Usigns, Lsigns, Svar, 0.2, 0.6);
    /*sign07.filter1(Usigns, Lsigns, Svar, 0.2, 0.7);
    sign08.filter1(Usigns, Lsigns, Svar, 0.2, 0.8);
    sign09.filter1(Usigns, Lsigns, Svar, 0.2, 0.9);*/

    std::ofstream filtered01;
    filtered01.open("filtered01.txt");
    sign01.write(filtered01);
    filtered01.close();

    std::ofstream filtered02;
    filtered02.open("filtered02.txt");
    sign02.write(filtered02);
    filtered02.close();

    std::ofstream filtered03;
    filtered03.open("filtered03.txt");
    sign03.write(filtered03);
    filtered03.close();

    std::ofstream filtered04;
    filtered04.open("filtered04.txt");
    sign04.write(filtered04);
    filtered04.close();

    std::ofstream filtered05;
    filtered05.open("filtered05.txt");
    sign05.write(filtered05);
    filtered05.close();

    std::ofstream filtered06;
    filtered06.open("filtered06.txt");
    sign06.write(filtered06);
    filtered06.close();

    /*std::ofstream filtered07;
    filtered07.open("filtered07.txt");
    sign07.write(filtered07);
    filtered07.close();

    std::ofstream filtered08;
    filtered08.open("filtered08.txt");
    sign08.write(filtered08);
    filtered08.close();

    std::ofstream filtered09;
    filtered09.open("filtered09.txt");
    sign09.write(filtered09);
    filtered09.close();*/

    /**
     * ====================================================================
     * Create binary matrix for connected component analysis
     * ====================================================================
     */
    now = std::time(0);
    std::cout << "[Elapsed: " << now - start << "s]\tCreating binary matrix" << std::endl;
    Aggregate binaryMatr1 = Aggregate(sign01);
    Aggregate binaryMatr2 = Aggregate(sign02);
    Aggregate binaryMatr3 = Aggregate(sign03);
    Aggregate binaryMatr4 = Aggregate(sign04);
    Aggregate binaryMatr5 = Aggregate(sign05);
    Aggregate binaryMatr6 = Aggregate(sign06);

    binaryMatr1.toBinary();
    std::ofstream binary1;
    binary1.open("BinaryMatrix01.txt");
    binaryMatr1.write(binary1);
    binary1.close();

    binaryMatr2.toBinary();
    std::ofstream binary2;
    binary2.open("BinaryMatrix02.txt");
    binaryMatr2.write(binary2);
    binary2.close();

    binaryMatr3.toBinary();
    std::ofstream binary3;
    binary3.open("BinaryMatrix03.txt");
    binaryMatr3.write(binary3);
    binary3.close();

    binaryMatr4.toBinary();
    std::ofstream binary4;
    binary4.open("BinaryMatrix04.txt");
    binaryMatr4.write(binary4);
    binary4.close();

    binaryMatr5.toBinary();
    std::ofstream binary5;
    binary5.open("BinaryMatrix05.txt");
    binaryMatr5.write(binary5);
    binary5.close();

    binaryMatr6.toBinary();
    std::ofstream binary6;
    binary6.open("BinaryMatrix06.txt");
    binaryMatr6.write(binary6);
    binary6.close();

    now = std::time(0);
    std::cout << "[Elapsed: " << now - start << "s]\tDone!" << std::endl;
    
    return 0;
}

#endif
