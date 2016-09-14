# Makefile for the RNA gap project

CC=g++
FLAGS=-std=c++0x -I./

HEADERS=$(wildcard *.hpp)
SOURCES=$(wildcard *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)

tester: tester.o closegap.o
	$(CC) $(FLAGS) -g -o tester tester.o closegap.o

tester.o: tester.cpp $(HEADERS)
	$(CC) $(FLAGS) -g -c tester.cpp

closegap.o: closegap.cpp $(HEADERS)
	$(CC) $(FLAGS) -g -c closegap.cpp

clean:
	rm -f *.o
	rm -f tester
