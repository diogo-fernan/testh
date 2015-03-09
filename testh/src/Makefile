# http://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/

# sudo bash
# export ARCHFLAGS='-arch i386 -arch x86_64'
# make

EXEC = ../bin/test

BDIR = ../bin
IDIR = ../include
SDIR = ../src
ODIR = obj

CC 		= gcc
CFLAGS 	= -I$(IDIR)
LIBS 	= -lm

_DEPS = func.h reg.h dist.h est.h gen.h stat.h rng.h proc.h util.h io.h batt.h sim.h
DEPS  = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = func.o reg.o dist.o est.o gen.o stat.o util.o io.o rng.o proc.o batt.o sim.o test.o 
OBJ  = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -g -c -Wall -o $@ $< $(CFLAGS)

$(EXEC): $(OBJ)
	$(CC) -g -Wall -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

run:
	./$(EXEC) $(h)

clean:
	rm -f  *~ $(ODIR)/*.o $(IDIR)/*~ 
	rm -f $(EXEC)

# check discussions below for input variables to make
# http://stackoverflow.com/questions/2826029/passing-additional-variables-from-command-line-to-make