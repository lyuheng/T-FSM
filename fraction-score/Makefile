CCOMPILE=g++
CPPFLAGS=-O3 -march=native
SRCS = run.cpp intersection/computesetintersection.cpp

all: run

run: run.cpp 
	$(CCOMPILE) -std=c++14 -fopenmp -lpthread  $(CPPFLAGS) $(SRCS) -o run

clean:
	-rm run
