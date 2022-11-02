CC=g++
FLAGS=-O3 -march=native -fopenmp -lpthread
DEPS=grami.h global.h gmatch.h

all: run

%.o: %.cpp $(DEPS)
	$(CC) -o $@ $< $(CFLAGS)

run: run.o computesetintersection.o	
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	-rm run







CC=gcc
CFLAGS=-I.
DEPS = hellomake.h
OBJ = hellomake.o hellofunc.o 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

hellomake: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)