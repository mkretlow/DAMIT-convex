MAKE=make
CC=gcc
LD=ld
AR=ar

OPTFLAGS=-O3
CFLAGS=$(OPTFLAGS) -Wall
LDFLAGS=-lm
MY_LIBS=blmatrix_direct.o phasec_direct.o matrix_direct.o memory.o bright_direct.o cross_product.o\
        dot_product.o

all: lcgenerator

libs: $(MY_LIBS)

lcgenerator: lcgenerator.c $(MY_LIBS)
	$(CC) $(CFLAGS) -o $@ $< $(MY_LIBS) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm *.o
	
