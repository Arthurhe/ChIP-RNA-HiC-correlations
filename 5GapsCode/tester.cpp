#include <fstream>
#include "closegap.hpp"

int main(int argc, char* argv[]) {
    //closegap(5);

    std::ifstream infile0;
    std::ifstream infile1;
    infile0.open(argv[1]);
    infile1.open(argv[2]);

    get_stats(infile0, infile1);

    infile0.close();
    infile1.close();
}
