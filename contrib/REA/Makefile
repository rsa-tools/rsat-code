#========================================================================
#  Makefile for REA.c
# 
#  Use: make REA
#
#  Authors: Victor Jimenez and Andres Marzal
#           Dept. de Informatica, Universitat Jaume I, Castellon, Spain
#           {vjimenez,amarzal}@inf.uji.es
#
#  Program REA - Version 1.1 - July 1999
#========================================================================


# Flags for debug:
CFLAGS = -g -I. -DDEBUG -DDEBUG_LOAD_GRAPH
CFLAGS = -g -I. -DDEBUG

# Flags for experimentation (maximum optimization level):
CFLAGS = -O6 -Wall -I. 

OBJS = loadgraph.o dijkstra.o chronometer.o 

#========================================================================
# Main program
#========================================================================

REA: REA.c REA.h $(OBJS)
	gcc REA.c $(CFLAGS) $(OBJS) -o $@

#========================================================================
# Auxiliary modules
#========================================================================

loadgraph.o: loadgraph.c loadgraph.h REA.h
	gcc loadgraph.c $(CFLAGS) -c -o $@

dijkstra.o: dijkstra.c dijkstra.h REA.h
	gcc dijkstra.c $(CFLAGS) -c -o $@

chronometer.o: chronometer.c chronometer.h
	gcc chronometer.c $(CFLAGS) -c -o $@


