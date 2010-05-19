CC      = gcc
CCFLAGS = -Wall -O3
INC     = -I../lib
OBJS    = main.o count.o ../lib/fasta.o ../lib/utils.o
APP     = word-analysis

$(APP):	$(OBJS)
	$(CC) $(LIBS) $(OBJS) -o $(APP)

%.o: %.c
	$(CC) -c $(CCFLAGS) $(INC) $< -o $@

clean:
	rm -f *.o $(APP)

all: 
	clean $(APP)
