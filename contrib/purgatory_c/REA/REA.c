/*========================================================================

  File: REA.c

  Program REA - Version 1.1 - July 1999

  ========================================================================

  This module contains an implementation of the Recursive Enumeration
  Algorithm (REA) that enumerates (by increasing weight) the N shortest
  paths in weighted graphs. The algorithm is described in: 

    "Computing the K Shortest Paths: a New Algorithm and an
    Experimental Comparison" by Victor Jimenez and Andres Marzal,
    3rd Workshop on Algorithm Engineering, London, July 1999.
    To be published by Springer-Verlag in the LNCS series.

  The sets of candidate paths are implemented using binary heaps.

  ========================================================================

    Copyright (C) 1999 Victor Jimenez and Andres Marzal

    MODIFIED BY Jean-Noël Monette in October 2005 (jmonette@info.ucl.ac.be)
    
    MODIFIED BY Jean-Noël Monette in September 2006 (jmonette@info.ucl.ac.be)

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

#define VERSION "1.1"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>  
#include <unistd.h>  

#ifdef DEBUG
#include <assert.h>
#endif

#include <REA.h>
#include <loadgraph.h>
#include <dijkstra.h>
#include <chronometer.h>


/*============================================================================
  GLOBAL VARIABLES
  ============================================================================*/

/* A buffer for heap elements is dinamically allocated by the main function  */
/* with size the number of arcs in the graph (which is the worst case bound) */
/* New heaps are initialized using this buffer by the function CreateHeap    */
Path **heapsBuffer;
Path **firstFreeHeap;

/* For experimental purposes we avoid to use the malloc system for new paths */
/* New paths are taken from this buffer by the function CreatePath           */
#define MAX_PATHS 10000000
Path pathsBuffer[MAX_PATHS];
int  numberUsedPaths = 0;
Path *firstFreePath = pathsBuffer;


/*============================================================================*/
__inline__ Path *CreatePath (Node *node, Path *backPath, COST_TYPE cost) {

  Path   *newPath;
  
#ifdef DEBUG
  assert (node != NULL);
#endif
  
  if (++numberUsedPaths > MAX_PATHS) {
    fprintf (stderr, "Redefine MAX_PATHS");
    exit(1);
  }
  newPath = firstFreePath++;
  
  newPath->backPath = backPath;
  newPath->cost     = cost;
  newPath->lastNode = node;
  
  return newPath;
  
}


/*============================================================================*/
void PrintPath (Path *path) {
  /* Prints the path in reverse order (with the final node first). */
  Path *backPath;
  
#ifdef DEBUG
  assert (path != NULL); 
#endif
 if (path->cost < INFINITY_COST) {
  for (backPath = path; backPath != NULL; backPath = backPath->backPath)
    printf ("%i-", backPath->lastNode->name);
  printf (" \t(Cost: %f)", path->cost);
 }
  
}



#define COST(_i_) (H->elt[_i_]->cost)

/*============================================================================*/
__inline__ void CreateHeap (Heap *H, int dimension) {
  /* Allocates memory for the heap elements from the heaps buffer 
     and initializes the number of element inside the heap to 0 */
#ifdef DEBUG
  assert (dimension >= 0);
  H->dimension = dimension;
#endif
  H->size = 0;
  H->elt = firstFreeHeap;
  firstFreeHeap += dimension;
}


/*============================================================================*/
__inline__ void PreInsertInHeap(Heap *H, Path *elt) {
  /* Inserts a new element in the heap in time O(1)
     without preserving the heap property */
#ifdef DEBUG
  assert (H->size < H->dimension);
#endif
  
  H->elt[H->size++] = elt;
}


/*============================================================================*/
__inline__ void BuildHeap(Heap *H) {
  /* Obtains the heap property in time linear with the heap size */
  
  register int i;
  for (i=H->size/2; i>=0; i--) {
    register Path * item = H->elt[i];
    register int j = 2*i+1;
    while (j<H->size) {
      if (j<H->size-1 && COST(j+1) < COST(j)) j++;
      if (item->cost <= COST(j)) break;
      H->elt[(j-1)/2] = H->elt[j];
      j = 2*j+1;
    }
    H->elt[(j-1)/2] = item;
  }
}


/*============================================================================*/
__inline__ Path * DeleteBestInHeap (Heap *H) {
  /* Returns the element with minimum cost in the heap and deletes
     it from the heap, restoring the heap property in worst case time
     logarithmic with the heap size */
  Path *best;
  
  if (H->size == 0) return (NULL);
  best = H->elt[0];
  H->elt[0] = H->elt[--H->size];
  {
    register Path * item = H->elt[0];
    register int j = 1;
    while (j<H->size) {
      if (j<H->size-1 && COST(j+1) < COST(j)) j++;
      if (item->cost <= COST(j)) break;
      H->elt[(j-1)/2] = H->elt[j];
      j = 2*j+1;
    }
    H->elt[(j-1)/2] = item;
  }
  return (best);
}

/*============================================================================*/
__inline__ Path * BestInHeap (Heap *H) {
     /* Returns the element with minimum cost in the heap without deleting it */

  if (H->size == 0) return (NULL);
  return H->elt[0];
}

/*============================================================================*/
__inline__ void UpdateBestInHeap (Heap *H, Path *elt) {
     /* Deletes the element with minimum cost in the heap and inserts a new one
	preserving the heap property, in worst case time logarithmic with the
	heap size */

  register Path * item = elt;
  register int j = 1;
  while (j<H->size) {
    if (j<H->size-1 && COST(j+1) < COST(j)) j++;
    if (item->cost <= COST(j)) break;
    H->elt[(j-1)/2] = H->elt[j];
    j = 2*j+1;
  }
  H->elt[(j-1)/2] = item;
}


/*============================================================================*/
__inline__ void InitializeCandidateSet (Node *node) {
  /* The set of candidates is  initialized with the best path from
     each  predecessor node, except  the  one from which  the best
     path at the current node  comes. If the current node is the 
     initial node in  the graph, all its predecessor nodes (if some
     exists) provide a candidate */

  PtrArc *ptrArc;
  int    numberArcs;
  Path   *path,
         *newCand;
  
  CreateHeap (&node->heap, node->numberArcsIn);
  
  for (ptrArc = node->firstArcIn, numberArcs = node->numberArcsIn;
       numberArcs != 0; ptrArc++, numberArcs--)
    if ( (path = (*ptrArc)->source->bestPath) != node->bestPath->backPath) {
      /* It is important to compare pointers and not costs or node names,
         because  several candidates could came from the same node (if
         the graph has parallel arcs) and/or with the same cost. */
      newCand = CreatePath (node, path, path->cost + (*ptrArc)->cost);
      PreInsertInHeap (&node->heap, newCand);
    }
  BuildHeap(&node->heap);
}

   
/*============================================================================*/
Path* NextPath (Path *path) {
  /* Central routine of the Recursive Enumeration Algorithm: computes the
     next path from the initial node to the same node in which the argument
     path ends, assuming that it has not been computed before.
  */
     
  Node      *node         = NULL;
  Path      *backPath     = NULL,
            *nextBackPath = NULL,
            *bestCand     = NULL,
            *newCand      = NULL;
  COST_TYPE arcCost;
  
#ifdef DEBUG
  assert (path != NULL);
  assert (path->nextPath == NULL);
  assert (path->lastNode != NULL);
#endif
  
  node = path->lastNode;
  
#ifdef DEBUG
  printf ("\nComputing next path at node %i\n", node->name);
#endif
  
  if (node->heap.elt == NULL) 
    /* This is done here instead of in Dijkstra function, so that it is */
    /* only done for nodes in which alternative paths are required.     */
    InitializeCandidateSet (node);
  
  if ((backPath = path->backPath) != NULL) {
    nextBackPath = backPath->nextPath;
    if (nextBackPath == NULL)
      nextBackPath = NextPath (backPath);
    if (nextBackPath != NULL) {
      arcCost = path->cost - backPath->cost;
      newCand = CreatePath (node, nextBackPath, nextBackPath->cost + arcCost);
    }
  }
  
  if (newCand == NULL)
    bestCand = DeleteBestInHeap (&node->heap);
  else {
    bestCand = BestInHeap (&node->heap);
    if (bestCand != NULL && bestCand->cost < newCand->cost)
      UpdateBestInHeap (&node->heap, newCand);
    else
      bestCand = newCand;
  }
  
  if (bestCand == NULL)
    return NULL;
  
  /* Adds the best candidate to the list of paths ending in the same node */
  bestCand->nextPath = NULL;
  path->nextPath     = bestCand;
  
#ifdef DEBUG
  printf ("\nNext path at node %i:\t", node->name);
  PrintPath (bestCand);
#endif
  
  return bestCand;
  
}

/*============================================================================*/
void Copyright () {
  printf ("\n======================================================================\n");
  printf ("REA version %s, Copyright (C) 1999 Victor Jimenez and Andres Marzal\n",
	  VERSION);
  printf ("REA comes with ABSOLUTELY NO WARRANTY.\n");
  printf ("This is free software, and you are welcome to redistribute it under\n");
  printf ("certain conditions; see the README and COPYING files for more details.\n");
}

/*============================================================================*/
void Help (char *program) {
  printf ("======================================================================\n");
  printf ("\nUse: %s GRAPH_FILE NUMBER_OF_PATHS [-paths] [-tdijkstra]\n", program);
  printf ("     Optional arguments:\n");
  printf ("       -paths Print the sequence of nodes for each path\n");
  printf ("       -tdijkstra Cumulate also the time of Dijkstra's algorithm\n");
  printf ("See the README file for more details.\n");
  printf ("======================================================================\n");
  exit (1);
}

/*============================================================================*/
int main (int argc, char **argv) {
 
  Graph  graph;
  Path   *path;
  int    i, numberPaths = 1;
  time_t date;
  struct tm *localDate;
  char   hostName[200] = "";
  int    showPaths = 0;
  int    measureDijkstra = 0;
  float  *cumulatedSeconds;
#ifdef DEBUG
  Node   *node;
  int    numberNodes;
#endif
  
  /****************** Prints the copyright notice  **************************/
  Copyright ();

  /************** Reads the command line parameters *************************/
  if (argc < 3) Help (argv[0]);
  numberPaths = atoi (argv[2]);
  for (i=3; i < argc; i++) {
    if (strcmp(argv[i], "-paths") == 0) showPaths = 1;
    else if (strcmp(argv[i], "-tdijkstra") == 0)  measureDijkstra = 1;
    else Help(argv[0]);
  }

  /****************** Prints experimental trace information *****************/
  gethostname (hostName, 200);
  date = time (NULL);
  localDate = localtime (&date);
  printf ("======================================================================\n");
  printf ("CommandLine: ");
  for (i = 0; i < argc; i++)
    printf (" %s", argv[i]);
  printf ("\nHostname: %s", hostName);
  printf ("\nDate: %s", asctime(localDate));
  printf ("======================================================================\n");

  
  /******************* Reads the graph from file ****************************/
  
  LoadGraph(&graph, argv[1]);
  
  /********** Allocates memory for heaps ************************************/
  /* Memory is allocated for the worst case: one element per arc */
  heapsBuffer = (Path **) malloc(graph.numberArcs*sizeof(*heapsBuffer));
  if (heapsBuffer == NULL) {
    perror("Not enough memory for heaps.\n");
    exit(1);
  }
  firstFreeHeap = heapsBuffer;
  
  /********** Allocates memory for time counters *****************************/
  cumulatedSeconds = malloc (sizeof(cumulatedSeconds[0])*numberPaths);
  if (cumulatedSeconds == NULL) {
    perror("Not enough memory for time counters.\n");
    exit(1);
  }
  
  /****************** Computes the shortest path tree ***********************/
  if (measureDijkstra == 1)
    ClockReset ();

  Dijkstra (&graph);

  if (measureDijkstra == 1)
    cumulatedSeconds[0] =  ClockTotal();
  else {
    cumulatedSeconds[0] =  0;
    ClockReset ();
  }

#ifdef DEBUG
  for (node = graph.node, numberNodes = graph.numberNodes; numberNodes != 0; 
       node++, numberNodes--) { 
    printf ("\nBest Path for FinalNode=%i: ", node->name);
    PrintPath (node->bestPath);
  }
#endif
  
  /******************** Computes the K shortest paths ***********************/
  i = 2;
  path = graph.finalNode->bestPath;
  while (i <= numberPaths && path != NULL) {
    printf("\nN=%i:\t", i-1);
    if (showPaths == 1) PrintPath (path);
    printf (" \t(CumulatedSeconds: %.2f)",0.5);
    path = NextPath (path);
    cumulatedSeconds[i-1] = ClockTotal();
    i++;
  }
  if (path != NULL){
    printf("\nN=%i:\t", i-1);
    if (showPaths == 1) PrintPath (path);
    printf (" \t(CumulatedSeconds: %.2f)",0.5);
  }
  
  /************ Prints the computed paths and time counters ******************/
  /*i = 1;
  path = graph.finalNode->bestPath;
  while (i <= numberPaths && path != NULL && path->cost < INFINITY_COST) {
    printf ("\nN=%i:\t", i);
    if (showPaths == 1) PrintPath (path);
    printf (" \t(CumulatedSeconds: %.2f)", (float) cumulatedSeconds[i-1]);
    path = path->nextPath;
    i++;
    }*/

  printf ("\nTotalExecutionTime: %.2f\n", (float) cumulatedSeconds[i-2]);
  return (0);
}
