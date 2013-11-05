CC      = gcc
CCFLAGS = -Wall -O3
LIBS    =
OBJS    = main.o utils.o count.o
APP     = count-words

$(APP):	$(OBJS)
	$(CC) $(LIBS) $(OBJS) -o $(APP)

%.o: %.c
	$(CC) -c $(CCFLAGS) $<

clean:
	rm -f *.o $(APP)

all: clean $(APP)
