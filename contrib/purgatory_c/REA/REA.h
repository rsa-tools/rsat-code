/*========================================================================

  File: REA.h

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


#ifndef _REA_H_INCLUDED

#define COST_TYPE float     /* You can change this type but you should verify
                             its use in printf's, scanf's and castings */

#define INFINITY_COST 100000000 /* Instead of INT_MAX, to avoid overflow  */
                                /* when adding some cost to INFINITY_COST */


/* The following data structures represent a graph in memory providing
   access to incoming and outgoing arcs for each node, and also allow
   to represent multiple shortest paths from the initial node to each node. 
*/

typedef struct Path {
  COST_TYPE   cost;      /* Path cost                                         */
  struct Node *lastNode; /* Node in which this path ends                      */
  struct Path *backPath; /* Prefix path, ending in a predecessor node         */
  struct Path *nextPath; /* Next path in the list of computed paths from the  */
                         /* initial node to the node in which this paths ends */
} Path;   

typedef struct Arc {
  COST_TYPE   cost;     /* Arc cost         */
  struct Node *source;  /* Source node      */
  struct Node *dest;    /* Destination node */
} Arc;  

typedef struct Arc *PtrArc;

typedef struct Heap {
#ifdef DEBUG
  int dimension; /* Allocated space                                         */
#endif
  int size;      /* Current number of elements                              */
  Path **elt;    /* Elements are pointers to paths, allocated by CreateHeap */
} Heap;

typedef struct Node {
  int    name;          /* An integer in the range 1..graph.numberNodes */
  int    numberArcsOut; /* Number of arcs that leave the node           */
  Arc    *firstArcOut;  /* Pointer to an element in graph.arc           */
  int    numberArcsIn;  /* Number of arcs that reach the node           */
  PtrArc *firstArcIn;   /* Pointer to an element in graph.arcIn         */
  Path   *bestPath;     /* First path in the list of computed paths     */
                        /* from the initial node to this node           */
  Heap   heap;          /* Set of candidate paths (binary heap)         */
} Node;

typedef struct Graph {
  int    numberNodes;  /* Number of nodes in the graph                        */
  int    numberArcs;   /* Number of arcs in the graph                         */
  Node   *initialNode; /* The node from which all paths depart                */
  Node   *finalNode;   /* The node to which the K shortest paths are required */
  Node   *node;        /* Nodes, with node named i in position i-1            */
  Arc    *arc;         /* Arcs sorted by source node                          */
  PtrArc *arcIn;       /* Pointers to arcs, sorted by destination node        */
}  Graph;

extern Path *CreatePath (Node *node, Path *backPath, COST_TYPE cost);

#define _REA_H_INCLUDED
#endif
