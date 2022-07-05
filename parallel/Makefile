CCOMPILE=g++
CPPFLAGS= -O2

all: run

run: run.cpp
	$(CCOMPILE) -g -std=c++11 -fopenmp -lpthread run.cpp $(CPPFLAGS)  -o run

clean:
	-rm run
