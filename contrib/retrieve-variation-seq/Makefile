CC = gcc
CFLAGS = -g -Wall -std=gnu11

retrieve-variation-seq:
	@echo ""
	@echo "Compiling retrieve-variation-seq"
	$(CC) $(CFLAGS) -o retrieve-variation-seq main.c
	@echo "	retrieve-variation-seq"

install:
	@echo ""
	@echo "Installing retrieve-variation-seq"
	rsync -ruptl retrieve-variation-seq ../../bin/retrieve-variation-seq

clean :
	rm -f retrieve-variation-seq

all:  clean retrieve-variation-seq install


