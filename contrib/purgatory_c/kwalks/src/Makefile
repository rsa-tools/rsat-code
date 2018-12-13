CC=gcc                          # used C-compiler
CFLAGS=-O3 -Wall		# release C-Compiler flags
LD=gcc                          # used linker
LFLAGS=-O3                	# linker flags
LIBS=-L. -lm                    # used libraries

all: lkwalk

clean:	
	rm -f *.o
	rm -f lkwalk

lkwalk: lkwalk.o tools.o
	$(LD) $(LDFLAGS) lkwalk.o tools.o -o ../bin/lkwalk $(LIBS)	

lkwalk.o: lkwalk.c tools.h
	$(CC) -c $(CFLAGS) lkwalk.c -o lkwalk.o

tools.o: tools.c tools.h
	$(CC) -c $(CFLAGS) tools.c -o tools.o

