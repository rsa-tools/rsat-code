CC      = g++ 
CCFLAGS = -O3 -Wall -stdlib=libc++
SRC     = $(wildcard *.cpp) 
OBJS    = $(SRC:.cpp=.o)
APP     = matrix-scan-quick

compile: $(OBJS)
	$(CC) $(OBJS) -lc++ -o $(APP)

%.o: %.cpp
	$(CC) -c $(CCFLAGS) $<

%.o: %.c
	$(CC) -c $(CCFLAGS) $<

clean:
	rm -f *.o $(APP)

all: clean compile

