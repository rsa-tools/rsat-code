#include <stdio.h>
#include <stdlib.h>

int n; /* The n number of nodes */
int directed; /*Specifies whether the network is directed or not. In case it is not directed the shortest path from node A to node B is the same as the shortest path from B to A*/
int size = 15000;
float dist[15000][15000]; /* dist[i][j] is the length of the edge between i and j if
			it exists, or 0 if it does not */
int rev[15000][15000];

void printDist() {
	int i, j;
	for (i = 0; i < n; ++i)
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j)
			printf("%.2f\t", dist[i][j]);
		printf("\n");
	}
	printf("\n");
}
void printRev() {
	int i, j;
	printf("    ");
	for (i = 0; i < n; ++i)
		printf("%4c", 'A' + i);
	printf("\n");
	for (i = 0; i < n; ++i) {
		printf("%4c", 'A' + i);
		for (j = 0; j < n; ++j)
		
			printf("%4d", rev[i][j]);
		printf("\n");
	}
	printf("\n");
}

 /*ReconstructPath ( from vertex i to j )        # finds the path from i to j
    If (shortest path from i to j represents the edge between them) Then
        Return a path composed of i and j
    Else
        Let k be the intermediate vertex of the path from i to j
        Return the path: [ ReconstructPath(i,k) ] concatenated with [ ReconstructPath(k,j) without first vertex k ]
   End If*/
   
   
void reconstructPath(int a, int b) {
    if (rev[a][b] > 0) {
      reconstructPath(a, rev[a][b]);
      printf("%d,", rev[a][b]);
      reconstructPath(rev[a][b], b);
    }
}

/*
	floyd_warshall()

	after calling this function dist[i][j] will the the minimum distance
	between i and j if it exists (i.e. if there's a path between i and j)
	or 0, otherwise
*/
void floyd_warshall() {
	int i, j, k;
	
	//setup variable
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            rev[i][j]=-1;
        }//End for
    }//End for
 
	
	
	
	for (k = 0; k < n; ++k) {
// 			printDist();

		for (i = 0; i < n; ++i) {
			for (j = 0; j < n; ++j) {
				/* If i and j are different nodes and if 
					the paths between i and k and between
					k and j exist, do */
				if ((dist[i][k] * dist[k][j] != 0) && (i != j)) {
					/* See if you can't get a shorter path
						between i and j by interspacing
						k somewhere along the current
						path */
					if ((dist[i][k] + dist[k][j] < dist[i][j]) || (dist[i][j] == 0))  {
						dist[i][j] = dist[i][k] + dist[k][j];
						rev[i][j] = k;
					}
				}
			}
		}
	}


	for (i = 0; i < n; ++i) {
	  j = 0;
	  if (directed == 0) {
	    j = i+1;
	  }
	  for (; j < n; ++j) {
	    if (j != i) {
	      printf("%d %d ", i,j);
	      if (dist[i][j] > 0) {
	        printf ("%d,",i);
	        reconstructPath(i,j);
	        printf ("%d",j);
	        printf (" %.2f", dist[i][j]);
	      } else {
	        printf ("No route");
	      }
	        printf ("\n");
	    }
	  }
	}
	
}

int main(int argc, char *argv[]) {
        if( argc < 2 ) {
          printf("One argument expected (input file).\n");
          printf("The input file must be as follows\n");
          printf("1st line : number of nodes in the network\n");
          printf("2d  line : 1 if the graph is directed, 0 if it is not directed\n");
          printf("3d->n+2 lines : The adjacency matrix with weights. 0 if no edge between two nodes\n");
          exit(0);
        }
	FILE *fin = fopen(argv[1], "r");
	fscanf(fin, "%d", &n);
	fscanf(fin, "%d", &directed);
	int i, j;
	for (i = 0; i < n; ++i)
		for (j = 0; j < n; ++j)
			fscanf(fin, "%f", &dist[i][j]);
	fclose(fin);

	floyd_warshall();

	return 0;
}

