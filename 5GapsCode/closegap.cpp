#ifndef __CLOSEGAP_CPP__
#define __CLOSEGAP_CPP__

#include <iostream>
#include <regex>

void closegap(int x) {
    std::cout << "Hello world!" << std::endl;
    std::cout << x << std::endl;
}

/**Input: annotation filestream, sequence filestream
 * Output: stats for start and end times
 * Method: Take all start positions and put them in a priority queue sorted by
 *         position. Pop off the queue into an intermediate variable and compare
 *         intermediate variable with number in bed file to determine nearest
 *         start. Collect differences between .bed start and annotation start
 *         in an accumulator and return statistics (object?).
 * Notes: "Start" and "End" in the annotation file are orientation (+/-)
 *        dependent and must be accounted for
 */
void get_stats(std::istream& ann_fstream, std::istream& seq_fstream) {
    std::regex expr("[0-9]*\t[0-9]*\t.*\t[+-]");
    std::cmatch m;

    char line[450];

    // Parse annotation file for start/end positions
    while(true) {
        ann_fstream.getline(line, 450);
        if(ann_fstream.eof()) {
            break;
        }

        // Search for the regex expr in the line
        std::regex_search(line, m, expr);

        char* substr;
        for(auto x:m) {
            
            std::cout << x;
        }
        std::cout << std::endl;
        
        
    }
}

#endif
