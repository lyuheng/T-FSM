CCOMPILE=g++
CPPFLAGS=-O3 -march=native

all: run

run: run.cpp
	$(CCOMPILE) -std=c++14 run.cpp intersection/computesetintersection.cpp $(CPPFLAGS)  -o run

clean:
	-rm run
