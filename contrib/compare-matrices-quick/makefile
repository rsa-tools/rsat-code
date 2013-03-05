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
	$(CC) $(OBJS) -o $(APP)

%.o: %.cpp
	$(CC) -c $(CCFLAGS) $<

%.o: %.c
	$(CC) -c $(CCFLAGS) $<

clean:
	rm -f *.o $(APP)

all: clean compile
