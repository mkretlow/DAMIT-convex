MAKE=make
CC=gcc
LD=ld
AR=ar

OPTFLAGS=-O3
CFLAGS=$(OPTFLAGS) -Wall
LDFLAGS=-lm
MY_LIBS=trifac.o areanorm.o sphfunc.o ellfit.o covsrt.o ludcmp.o lubksb.o mrqmin.o mrqcof.o\
        curv.o blmatrix.o conv.o gauss.o phasec.o matrix.o bright.o memory.o\
	dot_product.o

all: convexinv period_scan

libs: $(MY_LIBS)

convexinv: convexinv.c $(MY_LIBS)
	$(CC) $(CFLAGS) -o $@ $< $(MY_LIBS) $(LDFLAGS)

period_scan: period_scan.c $(MY_LIBS)
	$(CC) $(CFLAGS) -o $@ $< $(MY_LIBS) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm *.o
	
