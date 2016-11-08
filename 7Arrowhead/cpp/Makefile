CC=g++
FLAGS=-std=c++0x -I./

HEADERS=$(wildcard *.hpp)
SOURCES=$(wildcard *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)

testDPRCS.o: testDPRCS.cpp $(HEADERS)
	$(CC) $(FLAGS) -g -c testDPRCS.cpp

testDPRCS: testDPRCS.o RCS.o
	$(CC) $(FLAGS) -g -o testDPRCS.out testDPRCS.o RCS.o

testUTriangle.o: testUTriangle.cpp $(HEADERS)
	$(CC) $(FLAGS) -g -c testUTriangle.cpp

testUTriangle: RCS.o UTriangle.o testUTriangle.o
	$(CC) $(FLAGS) -g -o testUTriangle.out testUTriangle.o RCS.o UTriangle.o

%.o: %.cpp
	$(CC) $(FLAGS) -g -c -o $@ $<

clean:
	rm -f *.o
	rm -f *.out
