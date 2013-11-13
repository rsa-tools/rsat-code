/*========================================================================

  File: loadgraph.c

  Program REA - Version 1.1 - July 1999

  ========================================================================

  This module contains the function for reading a graph from a file
  into the data structures of the REA.

  ========================================================================

    Copyright (C) 1999 Victor Jimenez and Andres Marzal

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program (file COPYING); if not, write to the Free
    Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

    You can contact the authors at:

           Victor Jimenez and Andres Marzal
           Dept. de Informatica, Universitat Jaume I, Castellon, Spain
           {vjimenez,amarzal}@inf.uji.es

    Updated information about this software will be available at the
    following www address:

           http://terra.act.uji.es/REA
 
  ========================================================================= */

#include <stdio.h>
#include <stdlib.h>

#include <REA.h>

#define MAX_LINE_LENGTH 100


/*============================================================================*/
void LoadGraph(Graph *graph, char *fileName) {
  /* Reads a (multi)graph from file in time O(number of arcs) */
  FILE *file;
  char line[MAX_LINE_LENGTH];
  int initialNodeName, finalNodeName, nodeName, sourceName, destName,
      numberArcs;
  COST_TYPE cost;
  Node *node, *source, *dest;
  Arc  *arc;
#ifdef DEBUG_LOAD_GRAPH
  PtrArc *ptrArc;
#endif

  if ((file = fopen (fileName, "r")) == NULL) {
    perror ("Graph file does not exist");
    exit (1);
  }
  
  /* Scans the file header looking for the number of nodes,  */
  /* the number of arcs, the initial node and the final node */

  graph->numberNodes = graph->numberArcs = initialNodeName = finalNodeName = -1;

  while (fgets (line, MAX_LINE_LENGTH , file) != NULL) {
    if (line[0] == 'n')
      sscanf(line, "n%d\n", &graph->numberNodes);
    else if (line[0] == 'm')
      sscanf(line, "m%d\n", &graph->numberArcs);
    else if (line[0] == 's')
      sscanf(line, "s%d\n", &initialNodeName);
    else if (line[0] == 't')
      sscanf(line, "t%d\n", &finalNodeName);
    /* The declarations for n, m, s and t are supossed to be before any arc
       declaration. This restriction can be removed by removing the following
       else-if, but then this loop will consume more time */
    else if (line[0] == 'a') break;
  }

  if (graph->numberNodes <= 0) {
    fprintf (stderr, "\nIncorrect number of nodes in file %s\n", fileName);
    exit (1);
  }
  if (graph->numberArcs <= 0) {
    fprintf (stderr, "\nIncorrect number of arcs in file %s\n", fileName);
    exit (1);
  }
  if (initialNodeName <= 0 || initialNodeName > graph->numberNodes) {
    fprintf (stderr, "\nIncorrect initial node in file %s\n", fileName);
    exit (1);
  }
  if (finalNodeName <= 0 || finalNodeName > graph->numberNodes) {
    fprintf (stderr, "\nIncorrect final node in file %s\n", fileName);
    exit (1);
  }

#ifdef DEBUG_LOAD_GRAPH
  printf ("NumberNodes=%d\nNumberArcs=%d\nInitialNode=%d\nFinalNode=%d\n",
          graph->numberNodes, graph->numberArcs,
          initialNodeName, finalNodeName);
#endif

  /* Allocates memory for nodes and arcs */
  
  graph->node = malloc (sizeof (graph->node[0])*graph->numberNodes);
  if (graph->node == NULL) {
    perror("Not enough memory to load graph\n");
    exit(1);
  }
  graph->arc = malloc (sizeof (graph->arc[0])*graph->numberArcs);
  if (graph->arc == NULL) {
    perror("Not enough memory to load graph\n");
    exit(1);
  }
  graph->arcIn = malloc(sizeof ((graph->arcIn)[0])*graph->numberArcs);
  if (graph->arcIn == NULL) {
    perror("Not enough memory to load graph\n");
    exit(1);
  }
  
  /* Sets the pointers to the initial and final nodes */
  graph->finalNode   = graph->node + (finalNodeName - 1);
  graph->initialNode = graph->node + (initialNodeName - 1);
  
  
  /* Initializes for each node the counters of input and output arcs, */
  /* the heap of candidate paths, the best path and the node name     */
  for (node = graph->node, nodeName = 1; nodeName <= graph->numberNodes;
       node++, nodeName++) {
    node->numberArcsIn       = 0;
    node->numberArcsOut      = 0;
    node->name               = nodeName;
    node->bestPath           = CreatePath (node, NULL, INFINITY_COST);
    node->bestPath->nextPath = NULL;
    node->heap.elt           = NULL;
  }

  /* Scans the file looking for arcs and counts incoming and outgoing arcs */
  /* for each node                                                         */
  rewind (file);
  numberArcs = 0;
  while (fgets (line, MAX_LINE_LENGTH , file) != NULL) {
    if (line[0] == 'a') {
      if (sscanf(line, "a%d %d %f\n", &sourceName, &destName, &cost) != 3) {
	fprintf (stderr, "\nError in file %s in line:\n%s\n", fileName, line);
	exit (1);
      }
      if (sourceName <= 0 || sourceName > graph->numberNodes) {
	fprintf (stderr, "\nIncorrect node in file %s in line:\n%s\n",
                 fileName, line);
	exit (1);
      }
      if (destName <= 0 || destName > graph->numberNodes) {
	fprintf (stderr, "\nIncorrect node in file %s in line:\n%s\n",
                 fileName, line);
	exit (1);
      }
      graph->node[sourceName-1].numberArcsOut++;
      graph->node[destName-1].numberArcsIn++;
      numberArcs++;
    }
  }
  if (numberArcs != graph->numberArcs) {
    fprintf (stderr, "\nThe number of arcs found in %s differs from m\n",
             fileName);
    exit (1);
  }
    
  /* Initializes the pointer to the first output arc for every node */
  numberArcs = 0;
  for (node = graph->node, nodeName = 1; nodeName <= graph->numberNodes;
       node++, nodeName++) {
    node->firstArcOut = graph->arc + numberArcs;
    numberArcs += node->numberArcsOut;
  }

  /* Initializes the pointer to the first input arc for every node */
  numberArcs = 0;
  for (node = graph->node, nodeName = 1; nodeName <= graph->numberNodes;
       node++, nodeName++) {
    node->firstArcIn = graph->arcIn + numberArcs;
    numberArcs += node->numberArcsIn;
  }

  /* Scans again the file stores the arcs in memory */
  rewind (file);
  for (node = graph->node, nodeName = 1; nodeName <= graph->numberNodes;
       node++, nodeName++) {
    node->numberArcsIn  = 0;
    node->numberArcsOut = 0;
  }
  while (fgets (line, MAX_LINE_LENGTH , file) != NULL) {
    if (line[0] == 'a') {
      sscanf(line, "a%d %d %f\n", &sourceName, &destName, &cost);
      source = graph->node + (sourceName - 1);
      dest   = graph->node + (destName - 1);
      arc    = source->firstArcOut + source->numberArcsOut;
      arc->source = source;
      arc->dest   = dest;
      arc->cost   = cost;
      dest->firstArcIn[dest->numberArcsIn] = arc;
      source->numberArcsOut++;
      dest->numberArcsIn++;
    }
  }
  
  fclose (file);

#ifdef DEBUG_LOAD_GRAPH
  for (node = graph->node, nodeName = 1; nodeName <= graph->numberNodes;
       node++, nodeName++) {
    printf ("Node=%d:\tNumberOutputArcs=%d:\n",
            node->name, node->numberArcsOut);
    for (arc = node->firstArcOut, numberArcs = 1;
         numberArcs <= node->numberArcsOut; arc++, numberArcs++)
      printf ("\t\t(%i, %i, %f)\n", arc->source->name,
              arc->dest->name, arc->cost);
    printf ("\tNumberInputArcs=%d:\n", node->numberArcsIn);
    for (ptrArc = node->firstArcIn, numberArcs = 1;
         numberArcs <= node->numberArcsIn; ptrArc++, numberArcs++)
      printf ("\t\t(%i, %i, %f)\n", (*ptrArc)->source->name,
	      (*ptrArc)->dest->name, (*ptrArc)->cost);
  }
#endif

}

