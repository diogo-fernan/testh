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

_DEPS = regression.h distributions.h estimators.h generators.h statistics.h rng.h process.h utilities.h io.h batteries.h
DEPS  = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = regression.o distributions.o estimators.o generators.o statistics.o utilities.o io.o rng.o process.o batteries.o test.o 
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