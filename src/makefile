# You may need to alter the CPP variable to point to your C++ compiler and 
# the link command to include your thread and other required libraries.

CPP = g++
CC = gcc

CFLAGS = -O3 -Wall -DNDEBUG

SOURCES = me-PCR.cpp stsmatch.cpp fasta-io.cpp util.cpp stsmatch.h fasta-io.h util.h makefile

all : me-PCR

clean: 
	- rm *.o
	- rm me-PCR

# For AIX, you must use -Xlinker -bmaxdata:0x80000000 

me-PCR : me-PCR.cpp stsmatch.o fasta-io.o util.o
	$(CPP) $(CFLAGS) -lm -lpthread -o me-PCR me-PCR.cpp stsmatch.o fasta-io.o util.o


stsmatch.o : stsmatch.cpp stsmatch.h util.h
	$(CPP) $(CFLAGS) -c stsmatch.cpp

fasta-io.o : fasta-io.cpp fasta-io.h util.h
	$(CPP) $(CFLAGS) -c fasta-io.cpp

util.o : util.cpp util.h
	$(CPP) $(CFLAGS) -c util.cpp

distrib:
	- rm -rf distrib
	mkdir distrib
	mkdir distrib/doc
	mkdir distrib/bin
	mkdir distrib/src
	mkdir distrib/test
	cp $(SOURCES) distrib/src/.
