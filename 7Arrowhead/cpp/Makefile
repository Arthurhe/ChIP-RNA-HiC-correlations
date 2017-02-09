CC=g++
FLAGS=-std=c++0x -O3 -ftree-vectorize -ftree-vectorizer-verbose=7 -I./
# -fopt-info-vec-missed

HEADERS=$(wildcard *.hpp)
SOURCES=$(wildcard *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)

score.o: score.cpp $(HEADERS)
	$(CC) $(FLAGS) -g -c score.cpp

score: score.o RCS.o UTriangle.o LTriangle.o UVTriangle.o LVTriangle.o CumVar.o Aggregate.o
	$(CC) $(FLAGS) -g -o score.out score.o RCS.o UTriangle.o LTriangle.o UVTriangle.o LVTriangle.o CumVar.o Aggregate.o

testDPRCS.o: testDPRCS.cpp $(HEADERS)
	$(CC) $(FLAGS) -g -c testDPRCS.cpp

testDPRCS: testDPRCS.o RCS.o
	$(CC) $(FLAGS) -g -o testDPRCS.out testDPRCS.o RCS.o

testUTriangle.o: testUTriangle.cpp $(HEADERS)
	$(CC) $(FLAGS) -g -c testUTriangle.cpp

testUTriangle: RCS.o UTriangle.o testUTriangle.o
	$(CC) $(FLAGS) -g -o testUTriangle.out testUTriangle.o RCS.o UTriangle.o

testLTriangle.o: testLTriangle.cpp $(HEADERS)
	$(CC) $(FLAGS) -g -c testLTriangle.cpp

testLTriangle: RCS.o LTriangle.o testLTriangle.o
	$(CC) $(FLAGS) -g -o testLTriangle.out testLTriangle.o RCS.o LTriangle.o

testUVTriangle.o: testUVTriangle.cpp $(HEADERS)
	$(CC) $(FLAGS) -g -c testUVTriangle.cpp

testUVTriangle: RCS.o UVTriangle.o testUVTriangle.o
	$(CC) $(FLAGS) -g -o testUVTriangle.out testUVTriangle.o RCS.o UVTriangle.o

testLVTriangle.o: testLVTriangle.cpp $(HEADERS)
	$(CC) $(FLAGS) -g -c testLVTriangle.cpp

testLVTriangle: RCS.o LVTriangle.o testLVTriangle.o
	$(CC) $(FLAGS) -g -o testLVTriangle.out testLVTriangle.o RCS.o LVTriangle.o

%.o: %.cpp
	$(CC) $(FLAGS) -g -c -o $@ $<

clean:
	rm -f *.o
	rm -f *.out
