CC      = g++
CCFLAGS = -O3 -Wall
SRC		= $(wildcard *.cpp)
OBJS    = $(SRC:.cpp=.o)
	
compile: info-gibbs

info-gibbs:	$(OBJS)
	$(CC) $(OBJS) -o info-gibbs

%.o: %.cpp
	$(CC) -c $(CCFLAGS) $<
	
clean:
	rm -f *.o info-gibbs *.ic

all: clean compile
	
zip:
	@cd ../; zip -q -r info-gibbs.zip info-gibbs/README info-gibbs/Makefile \
	info-gibbs/*.cpp info-gibbs/*.h

test:
	time ./info-gibbs -i data/LexA.fa -w 8 --iter=100 --collect -v 1
	time ./info-gibbs -i data/LexA.fa --iter=100 -w 8 --bgfile=data/MM2.bg	
	time ./info-gibbs -i data/LexA.fa --iter=100 --bgfile=data/MM2.bg --seedmatrix=data/random.mat	

