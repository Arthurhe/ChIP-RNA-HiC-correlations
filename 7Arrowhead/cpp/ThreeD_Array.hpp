#ifndef __THREED_ARRAY_HPP__
#define __THREED_ARRAY_HPP__

#include <cassert>

// This class will create an N * M * O array of type T for you since 
// multi dimensional arrays are slow in C++
template <typename T>
class ThreeD_Array {
public:
    // o represents depth
    // n represents number of rows,
    // m represents number of columns
    ThreeD_Array (int o, int n, int m): o(o), n(n), m(m) {
        linear_array = new T*[o];
        for(int i = 0; i < o; ++i) {
            linear_array[i] = new T[n * m];
        }
    }

    ~ThreeD_Array() {
        delete [] linear_array;
        linear_array = 0;
    }

    // you can use 'at' to both set and retrieve values in the 3D_Array
    T& at (int d, int r, int c) {
        // Boundary check
        assert(d < o && d >= 0);
        assert(r < n && r >= 0);
        assert(c < m && c >= 0);
        return linear_array[d][(r*m) + c];
    }

    // returns number of rows
    int getNumRows() {
        return n;
    }

    // returns number of columns
    int getNumCols() {
        return m;
    }

    // returns the depth
    int getDepth() {
        return o;
    }

    void printOut(int depth) {
        assert(depth < o);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                std::cout << at(depth, i,j) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

private:
    T** linear_array;
    int n;
    int m;
    int o;
};

#endif
