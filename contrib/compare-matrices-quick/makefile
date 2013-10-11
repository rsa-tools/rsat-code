################################################################
## Compile compare-matrices-quick
##
## Author: Sebastien JAEGER <jaeger@ciml.univ-mrs.fr>

#CC      = g++ 
CC      = gcc
CCFLAGS = -O3 -Wall
SRC     = $(wildcard *.cpp) 
CSRC    = $(wildcard *.c)
OBJS    = $(SRC:.cpp=.o) $(CSRC:.c=.o)
APP     = compare-matrices-quick

compile: $(OBJS)
	$(CC) $(OBJS) -o $(APP) -lm

%.o: %.cpp
	$(CC) -c $(CCFLAGS) $<

%.o: %.c
	$(CC) -c $(CCFLAGS) $<

clean:
	rm -f *.o $(APP)

all: clean compile

################################################################
## Test compare-matrices-quick (this C version)

test: test_quick test_slow

LTH_W=3
LTH_NCOR2=0.7
MATRIX_DB=${RSAT}/data/motif_databases/JASPAR/jaspar_core_vertebrates_2009_10.tf
QUERY_MATRICES=DemoCompMat.txt
RESULT=result
test_quick:
	time compare-matrices-quick -o ${RESULT}_quick.tab -file1 ${QUERY_MATRICES} -file2 ${MATRIX_DB} -lth_ncor2 ${LTH_NCOR2} -lth_w ${LTH_W}
	@echo "	-o ${RESULT}_quick.tab "
	text-to-html -i ${RESULT}_quick.tab -o ${RESULT}_quick.html
	@echo "	-o ${RESULT}_quick.html "

## Test compare-matrices (the original RSAT Perl script)
V=1
test_slow:
	time compare-matrices -v ${V} -file1 ${QUERY_MATRICES} -file2 ${MATRIX_DB} -format transfac \
		-lth Ncor2 ${LTH_NCOR2} -lth w ${LTH_W} \
		-mode profiles \
		-return matrix_id,cor,Ncor,Ncor1,Ncor2,w,offset,strand \
		-o ${RESULT}_slow.tab 
	@echo "	-o ${RESULT}_slow.tab "
