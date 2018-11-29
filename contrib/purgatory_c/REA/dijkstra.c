/*========================================================================

  File: dijkstra.c

  Program REA - Version 1.1 - July 1999

  ========================================================================

  This module contains an implementation of Dijkstra's algorithm for
  computing the shortest path from the initial node to every other node
  in a graph.

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

#include <stdlib.h>
#include <stdio.h>

#ifdef DEBUG
#include <assert.h>
#endif

#include <REA.h>



/* The following data structures and fuctions implement a binary
   heap including the possibility of improving the priority of an
   element that is in the heap. The are used by Dijkstra function. 
*/ 
typedef struct {
  COST_TYPE cost;      /* Priority                               */
  int       name;      /* Identifier in the range 0..dimension-1 */
  int       position;  /* Position in the array of elements      */
} DijkstraHeapElement;

typedef struct  {
#ifdef DEBUG
  int                 dimension; /* Allocated space            */
#endif
  int                 size;      /* Current number of elements */
  DijkstraHeapElement *info;     /* Information fields, indexed by name */
  DijkstraHeapElement **elt;     /* Pointers to the information fields  */
} DijkstraHeap;


#define NOT_IN_HEAP -1

#define PRIORITY(_i_) (heap->elt[_i_]->cost)


/*============================================================================*/
__inline__ void CreateDijkstraHeap (DijkstraHeap * heap, int dimension) {

  int i;
#ifdef DEBUG   
  heap->dimension = dimension;
#endif
  
  heap->size = 0;
  /* Space is allocated for dimension+1 elements so that they are inserted
     starting at position 1 (there is not special reason for it) */
  heap->elt = malloc ((dimension+1)*sizeof(heap->elt[0]));
  if (heap->elt == NULL) {
    perror("Not enough memory for Dijkstra heap\n");
    exit(1);
  }
  heap->info = malloc(dimension*sizeof(heap->info[0]));
  if (heap->info == NULL) {
    perror("Not enough memory for Dijkstra heap\n");
    exit(1);
  }
  for (i = 0; i < dimension; i++) 
    heap->info[i].position = NOT_IN_HEAP;
}


/*============================================================================*/
__inline__ void FreeDijkstraHeap (DijkstraHeap * heap) {
  
  free(heap->info);
  free(heap->elt);
}


/*============================================================================*/
__inline__ void InsertInDijkstraHeap (DijkstraHeap * heap, COST_TYPE cost, int name) {
  
  int position, parent;
  DijkstraHeapElement *elt;

#ifdef DEBUG
  assert (heap->size < heap->dimension);
#endif

  heap->size++;
  position = heap->size;
  parent   = position>>1; /* position/2 */
  while (parent >= 1 && PRIORITY(parent) > cost) {
    heap->elt[position]           = heap->elt[parent];
    heap->elt[position]->position = position;
    position = parent;
    parent  = position>>1; /* position/2 */
  }
  elt = heap->info + name;
  heap->elt[position] = elt;
  elt->position       = position;
  elt->cost           = cost;
  elt->name           = name;

}


/*============================================================================*/
__inline__ COST_TYPE DeleteBestInDijkstraHeap (DijkstraHeap * heap, int * name) {
  
  COST_TYPE bestCost;
  DijkstraHeapElement *item;
  int position, parent;
  
#ifdef DEBUG  
  assert (heap->size > 0);
#endif
  
  *name = heap->elt[1]->name;
  heap->elt[1]->position = NOT_IN_HEAP;
  bestCost = PRIORITY(1);
  
  item = heap->elt[heap->size];
  heap->size--;
  
  position = 2;
  while (position<=heap->size) {
    if (position<=heap->size-1 && PRIORITY(position+1) < PRIORITY(position)) position++;
    if (item->cost <= PRIORITY(position)) break;
    parent = position>>1;   /* position/2 */
    heap->elt[parent]           = heap->elt[position];
    heap->elt[parent]->position = parent;
    position = position<<1; /* 2*position */
  }
  parent = position>>1;  /* position/2 */
  heap->elt[parent]           = item;
  heap->elt[parent]->position = parent;
  
  return (bestCost);
}


/*============================================================================*/
__inline__ int BelongsToDijkstraHeap(DijkstraHeap *heap, int name) {
  
  return (heap->info[name].position != NOT_IN_HEAP);
}


/*============================================================================*/
__inline__ void DecreaseCostInDijkstraHeap(DijkstraHeap *heap, COST_TYPE cost, int name) {
  
  int parent, position;

  position = heap->info[name].position;

#ifdef DEBUG 
  assert (position != NOT_IN_HEAP);
  assert (cost < PRIORITY(position));
#endif

  parent   = position>>1;  /* position/2 */
  while (parent >= 1 && PRIORITY(parent) > cost) {
    heap->elt[position]           = heap->elt[parent];
    heap->elt[position]->position = position;
    position = parent;
    parent = position>>1;  /* position/2 */
  }
  heap->elt[position]           = heap->info + name;
  heap->elt[position]->position = position;
  heap->elt[position]->cost = cost;

}

/*============================================================================*/
void Dijkstra (Graph *graph) {
  /* Computes the tree with the best path from the initial node to
     every node in the graph. Works for graphs with positive weigths.
     node->bestPath was initialized by LoadGraph for every node.
  */
  
  COST_TYPE    bestCost;   
  int          indexNode, numberArcs;
  Node         *node;
  DijkstraHeap heap;
  Arc          *arc;
  
  CreateDijkstraHeap (&heap, graph->numberNodes);
  
#ifdef DEBUG
  assert (graph->initialNode != NULL);
#endif
  
  graph->initialNode->bestPath->cost = 0;
  indexNode = graph->initialNode->name - 1;
  InsertInDijkstraHeap (&heap, 0, indexNode);
  
  while ( heap.size > 0 ) {
    bestCost = DeleteBestInDijkstraHeap (&heap, &indexNode);
    node = graph->node + indexNode;
    for (arc = node->firstArcOut, numberArcs = node->numberArcsOut;
         numberArcs != 0; arc++, numberArcs--) {
      Node *dest = arc->dest;
      if (bestCost + arc->cost < dest->bestPath->cost) {
	indexNode = dest->name - 1;
	if (BelongsToDijkstraHeap (&heap, indexNode)) 
	  DecreaseCostInDijkstraHeap (&heap, bestCost + arc->cost, indexNode); 
	else 
	  InsertInDijkstraHeap (&heap, bestCost + arc->cost, indexNode);
	dest->bestPath->cost     = bestCost + arc->cost;
	dest->bestPath->backPath = node->bestPath;
      }
    }
  }
  
  FreeDijkstraHeap (&heap);   

}
