CC=g++
CFLAGS=-std=c++14 -O3 -march=native
OBJDIR=objs/
OBJS=$(OBJDIR)intersection.o 
INC=canonical.h decompose.h gmatch.h grami.h graph.h leapfrogjoin.h pretty_print.h setting.h types.h

all: run

run: $(OBJDIR)run.o $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@ 

$(OBJDIR)intersection.o: intersection/computesetintersection.cpp intersection/computesetintersection.h 
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)run.o: run.cpp $(INC)
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf run $(OBJDIR)


MKDIR_P = mkdir -p
