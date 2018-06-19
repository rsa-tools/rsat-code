CC = gcc
CFLAGS = -g -std=gnu11

variation-scan:
	@echo ""
	@echo "Compiling retrieve-variation-seq"
	$(CC) $(CFLAGS) -o variation-scan main.c

install:
	@echo ""
	@echo "Installing retrieve-variation-seq"
	rsync -ruptl variation-scan ../../bin/variation-scan

clean :
	rm -f variation-scan


all:  clean variation-scan install
