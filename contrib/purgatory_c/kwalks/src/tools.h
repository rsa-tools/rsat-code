/*************************************************************************/
/* Authors     : Jerome Callut (jerome.callut@uclouvain.be)              */
/*               Pierre Dupont (pierre.dupont@uclouvain.be)              */
/* Institution : Departement d'Ingenierie Informatique                   */
/*               Universite catholique de Louvain                        */
/* Ref.        : P. Dupont, J. Callut, G. Dooms, J.-N. Monette, and      */
/*               Y. Deville. Relevant subgraph extraction from random    */
/*               walks in a graph. Technical Report RR 2006-07, INGI,    */
/*               UCL, 2006.                                              */
/* Version     : 1.04                                                    */
/* Modified    : 09/10/2010                                              */
/* License     : GPL 3.0                                                 */ 
/*               GNU GENERAL PUBLIC LICENSE                              */
/*               Version 3, 29 June 2007                                 */
/*               For details, see the file gpl-3.0.txt in this bundle    */
/*************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <search.h>
#include <time.h>
#include <getopt.h>

/*CONSTANTS FOR TOOLS.C*/
#define BUF_LINE_LEN   150000    /* Size of the buffer for file reading */
#define FNAME_LEN      1000      /* Maximum size of file names*/ 
#define PROB_EPS       1e-4      /* Precision to check if probas sum up to 1*/         
#define VERSION        1.04      /* Current software version*/

/*ARITHMETIC MACROS*/
#define log2(x) (log10((double)(x))/log10(2.0))
#define min(x,y) x < y ? x : y
#define ABS(x) (((x) < 0) ? -(x) : (x))


/*****************************DATA STRUCTURES*****************************/

typedef struct SparseDim{
    int    nbElems;
    int    absorbing;
    double starting;
    double sum;
    struct SparseElem *first;
    struct SparseElem *last;
    struct Elem *firstAbs;
    struct Elem *lastAbs;
}SparseDim;

typedef struct SparseElem{
    int    l;
    int    c;
    double val;
    double ept;
    double cur;
    double diff;    
    struct SparseElem *nextL;
    struct SparseElem *nextC;
}SparseElem;

typedef struct Elem{
    struct SparseElem *selem;
    struct Elem *next;
}Elem;

typedef struct DegStat{
    int    minDeg;
    int    maxDeg;
    double meanDeg;
}DegStat;

/****************MEMORY(DE)ALLOCATIONS AND INITIALIZATION*****************/

/* Allocate a vector of doubles */
double *vec_alloc(int n);

/* Allocate a vector of integers */
int *int_vec_alloc(int n);

/* Allocate a matrix of doubles */
double **mat_alloc(int n, int k);

/* Allocate a matrix of integers */
int **int_mat_alloc(int n, int k);

/* Allocate a matrix of char */
char **char_mat_alloc(int n,int k);

/* Allocate a tridimensional array of doubles */
double ***mat3_alloc(int n, int k, int v);

/* Allocate a tridimensional array of integers */
int ***int_mat3_alloc(int n, int k, int v);

/* Deallocate a matrix of doubles*/
void free_mat(double **matrix, int dim1);

/* Deallocate a matrix of integers */
void free_int_mat(int **matrix, int dim1);

/* Deallocate a tridimensional array of doubles */
void free_mat3(double ***matrix, int dim1, int dim2);

/* Deallocate a tridimensional array of integers*/
void free_int_mat3(int ***matrix, int dim1, int dim2);

/* Initialize a vector of doubles with zeros */
void init_vec(double *vector, int dim1);

/* Initialize a vector of integers with zeros */
void init_int_vec(int *vector, int dim);

/* Initialize a matrix of doubles with zeros */
void init_mat(double **matrix, int dim1, int dim2);

/* Initialize a matrix of integers with zeros */
void init_int_mat(int **matrix, int dim1, int dim2);

/* Initialize a tridimensional matrix of integers with zeros */
void init_int_mat3(int ***matrix, int dim1, int dim2,int dim3);

/*Initialize an array of SparseDim*/
void init_SparseDims(SparseDim* sdim,int dim1);

/****************************I/O OPERATIONS*******************************/

/*Tell if a file contains sparse data*/
int isSparseFile(char* fname);

/* Read I in input file */
int readVec(char* fname,double* I,int m);

/* Read P in input file */
int readMat(char* fname,double** P,int m);

/* Read matrix of double in input file */
int readMat_sparse(char* fname,SparseDim* rows,SparseDim* cols,int m,DegStat* degStat);

/* Read the number of nodes from the adjacency matrix file */
int readNbNodes(char* fname);

/* Write a vector of double in a file*/
void writeVec(double* I,int m,char* fname);

/* Write a vector of double in a sparse file*/
void writeVec_sparse(double* I,int m,char* fname);

/* Write a matrix of double in a file*/
void writeMat(double** P,int m,char* fname);

/*Write the MEPT in a sparse file*/
void writeMEPT_sparse(char* fname,SparseDim* rows,int m);

/*Write the absolute MEPT difference in a sparse file*/
void writeDiff_sparse(char* fname,SparseDim* rows,int m);

/*Export a matrix of MEPN to the graphviz dot format*/
void writeDot(char* fname,SparseDim* rows,int m,int undirected);

/****************************Miscellaneous*******************************/

/* Print command line help for the limited kwalks*/
void printHelpKwalk();

/* Display a vector of integers */
void show_int_vec(int* vec,int n);

/* Display a matrix of integers */
void show_int_mat(int** mat,int l,int c);

/* Display a vector of doubles */
void show_vec(double* vec,int n);

/* Display a matrix of doubles */
void show_int_mat(int** mat,int l,int c);

/* Display a matrix of doubles*/
void show_mat(double** mat,int l,int c);

/* Display an array of SparseDims */
void show_SparseDims(SparseDim* sdim,int m,char mode);

/*Tokenize a list of filenames separated by commas*/
int tokenizeFnames(char* list,char** fnames);

/*Add an extension to a file name*/
void addExtension(char* fname,char* extension,char* target);

/*Tokenize a list of integers separated by :*/
int tokenizeInt(char* ilist,int* olist);

/*Tokenize a list of doubles separated by :*/
int tokenizeDouble(char* ilist,double* olist);

/*Return the number of groups and the max number of nodes per group*/
void getGroupsInfo(char* list,int *nbGroups, int *maxNPG);

/*Tokenize groups of relevant nodes*/
int tokenizeKGroups(char* list,int** kgroups);

/*Tokenize initial probabilities of groups*/
int tokenizeKGProbas(char* list,double** kgprobas);

/* Check whether the maximal size of the reading buffer is reached */
void check_fgets(char *s, int limit, int iline, char *filename);

/****************************Matrix operations****************************/

/* Set a vector with element i equal 1 and others equal 0 */
void setE_i(double* vec,int m,int i);

/* Compute the sum of a vector of doubles*/
double vecSum(double* vec,int m);

/* Compute the sum of a vector of doubles*/
int int_vecSum(int* vec,int m);

/*Divide all elements of a vector of doubles by a number a*/
void vecNormalize(double* vec,int m,double a);

/*Divide all elements of a matrix of doubles by a number a*/
void matNormalize(double** P,int m,double a);

/*Divide all elements of a sparse matrix by a number a*/
void matNormalize_sparse(SparseDim* rows,int m,double a);

/*Compute the position of the maximum element in a vector of doubles*/
int max_vec_pos(double* vec,int m);

/*Enforce the stochasticity of a matrix of positive double */
void stochMat(double** P,int m);

/*Enforce the stochasticity of a sparse matrix of positive double*/
void stochMat_sparse(SparseDim* rows,int m);

/*Enforce the stochasticity of a vector of positive double */
double stochVector(double* vec,int m);

/* Build a forward lattice up to wlen*/
void buildLatticeForward(double** lat,int wlen,int* kgroup,double* kgproba,double** P,int m);

/* Build a forward lattice up to wlen*/
void buildLatticeForward_sparse(double** lat,int wlen,int* kgroup,double* kgproba,SparseDim* cols,int m,int verbose);
    
/* Build a backward lattice up to wlen*/
void buildLatticeBackward(double** lat,int wlen,double** P,int m);

/* Build a backward lattice up to wlen*/
void buildLatticeBackward_sparse(double** lat,int wlen,SparseDim* rows,int m,int verbose);

/*Compute the cumulated absorption probability*/
double cumulatedAbsorptionProba(double** latB,int wlen,int* kgroup,double* kgproba);

/*Set all the groups being absorbing except the starting group*/
void setAbsStates(double** P,double** PStart,int m,int** kgroups,int nbGroups,int start);

/*Set all the groups being absorbing except the starting group*/
void setAbsStates_sparse(SparseDim* rows,SparseDim* cols,int m,int** kgroups,double** kgprobas,int nbGroups,int start);

/*Compute the absolute EPT difference : |E_ij - E_ji|*/
void diffEPT(double** E,int m);

/*Compute the absolute EPT difference : |E_ij - E_ji|*/
void diffEPT_sparse(SparseDim* rows,SparseDim* cols,int m);

/*Test if a matrix is symmetric*/
int isSymmetric(double** P,int m);

/*Test if a sparse matrix is symmetric*/
int isSymmetric_sparse(SparseDim* rows,SparseDim* cols,int m);

/*Get the number of non-nul elements in a matrix*/
int matrixNNZ(double** P,int m);

/*Get the number of non-nul elements in a matrix*/
void matrixNNZ_sparse(SparseDim* rows,int m,int* nnz,int* nnr);

/*Copy edge passage times to edge value for inflation*/
void copyForInflation(SparseDim* rows,int undirected,int m);

/*Copy a vetcor into another (must be allocated)*/
void copy_mat(double** src,double** dest,int l,int c);

/*************************************************************************/
