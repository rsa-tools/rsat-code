CC      = gcc
CCFLAGS = -Wall -O3
INC     = -I../lib
OBJS    = main.o count.o ../lib/fasta.o ../lib/utils.o ../lib/markov.o ../lib/binomial.o 
APP     = word-analysis
LIBS    = -lm

$(APP):	$(OBJS)
	$(CC) $(LIBS) $(OBJS) -o $(APP)

%.o: %.c
	$(CC) -c $(CCFLAGS) $(INC) $< -o $@

clean:
	rm -f *.o $(APP)
	rm -f ../lib/*.o $(APP)

all: 
	$(MAKE) clean
	$(MAKE) $(APP)
