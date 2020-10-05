MAKE=make
CC=gcc
LD=ld
AR=ar

OPTFLAGS=-O3
CFLAGS=$(OPTFLAGS) -Wall
LDFLAGS=-lm
MY_LIBS=trifac.o blmatrix_cg.o conv_cg.o phasec_cg.o areanorm_cg.o\
	memory.o dot_product.o bright_all.o bright_cg.o\
	matrix_cg.o dbright_all.o f1dim.o df1dim.o bright_cg_deriv.o\
	frprmn.o dlinmin.o mnbrak.o dbrent.o conv_cg_deriv.o

all: conjgradinv

libs: $(MY_LIBS)

conjgradinv: conjgradinv.c $(MY_LIBS)
	$(CC) $(CFLAGS) -o $@ $< $(MY_LIBS) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm *.o
	
