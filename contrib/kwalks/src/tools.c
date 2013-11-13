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
/* Modified    : 09/10/2010 by P.Dupont                                  */
/*               - Fix overflow issue with fgets                         */
/* License     : GPL 3.0                                                 */ 
/*               GNU GENERAL PUBLIC LICENSE                              */
/*               Version 3, 29 June 2007                                 */
/*               For details, see the file gpl-3.0.txt in this bundle    */
/*************************************************************************/

#include "tools.h"

/****************MEMORY(DE)ALLOCATIONS AND INITIALIZATION*****************/

/* Allocate a vector of doubles */
double *vec_alloc(int n)
{
    double *v;
    v=(double *) calloc(n,sizeof(double));
    if(v==NULL) {
        fprintf(stderr,"could not allocate memory");
        exit(EXIT_FAILURE);
    }
    return v;
}

/* Allocate a vector of integers */
int *int_vec_alloc(int n)
{
    int *v;
    v=(int *) calloc(n,sizeof(int));
    if(v==NULL) {
        fprintf(stderr,"could not allocate memory");
        exit(EXIT_FAILURE);
    }
    return v;
}

/* Allocate a matrix of doubles */
double **mat_alloc(int n, int k)
{
    int i;
    double **mat;
    mat=(double **) calloc(n,sizeof(double *));
    if(mat == NULL) {
        fprintf(stderr,"could not allocate memory");
        exit(EXIT_FAILURE);
    }
    for(i=0; i<n; i++) {
        mat[i]=(double *) calloc(k,sizeof(double));
        if(mat[i] == NULL) {
            fprintf(stderr,"could not allocate memory");
            exit(EXIT_FAILURE);
        }
    }
    return mat;
}

/* Allocate a matrix of integers */
int **int_mat_alloc(int n, int k)
{
    int i, **mat;
    mat=(int **) calloc(n,sizeof(int *));
    if(mat == NULL) {
        fprintf(stderr,"could not allocate memory");
        exit(EXIT_FAILURE);
    }
    for(i=0; i<n; i++) {
        mat[i]=(int *) calloc(k,sizeof(int));
        if(mat[i] == NULL) {
            fprintf(stderr,"could not allocate memory");
            exit(EXIT_FAILURE);
        }
    }
    return mat;
}

/* Allocate a matrix of char */
char **char_mat_alloc(int n,int k){
    char **mat;
    mat=(char **) calloc(n,sizeof(char *));
    if(mat == NULL) {
        fprintf(stderr,"could not allocate memory");
        exit(EXIT_FAILURE);
    }
    return mat;    
}    

/* Allocate a tridimensional array of doubles */
double ***mat3_alloc(int n, int k, int v)
{
  int i,j;
  double ***mat;
  mat=(double ***) calloc(n,sizeof(double **));
  if(mat == NULL) {
    fprintf(stderr,"could not allocate memory");
    exit(EXIT_FAILURE);
  }
  for(i=0; i<n; i++) {
    mat[i]=(double **) calloc(k,sizeof(double *));
    if(mat[i] == NULL) {
      fprintf(stderr,"could not allocate memory");
      exit(EXIT_FAILURE);
    }
    for(j=0; j<k; j++) {
      mat[i][j]=(double *) calloc(v,sizeof(double));
      if(mat[i][j] == NULL) {
	fprintf(stderr,"could not allocate memory");
	exit(EXIT_FAILURE);
      }
    }
  }
  return mat;
}

/* Allocate a tridimensional array of integers */
int ***int_mat3_alloc(int n, int k, int v)
{
    int i,j;
    int ***mat;
    mat=(int ***) calloc(n,sizeof(int **));
    if(mat == NULL) {
        fprintf(stderr,"could not allocate memory");
        exit(EXIT_FAILURE);
    }
    for(i=0; i<n; i++) {
        mat[i]=(int **) calloc(k,sizeof(int *));
        if(mat[i] == NULL) {
            fprintf(stderr,"could not allocate memory");
            exit(EXIT_FAILURE);
        }
        for(j=0; j<k; j++) {
            mat[i][j]=(int *) calloc(v,sizeof(int));
            if(mat[i][j] == NULL) {
                fprintf(stderr,"could not allocate memory");
                exit(EXIT_FAILURE);
            }
        }
    }
    return mat;
}


/* Deallocate a matrix of doubles*/
void free_mat(double **matrix, int dim1)
{
    int i;
    
    for (i=0; i<dim1; i++)
        free(matrix[i]);
    free(matrix);
}

/* Deallocate a matrix of integers */
void free_int_mat(int **matrix, int dim1)
{
    int i;
    
    for (i=dim1-1; i>=0; i--)
        free(matrix[i]);
    free(matrix);
}

/* Deallocate a tridimensional array of doubles */
void free_mat3(double ***matrix, int dim1, int dim2)
{
  int i;

  for (i=0; i<dim1; i++) 
    free_mat(matrix[i],dim2);
  free(matrix);
}

/* Deallocate a tridimensional array of integers */
void free_int_mat3(int ***matrix, int dim1, int dim2)
{
  int i;

  for (i=0; i<dim1; i++) 
    free_int_mat(matrix[i],dim2);
  free(matrix);
}

/* Initialize a vector of doubles with zeros */
void init_vec(double *vector, int dim1)
{
    int i;
    for(i=0; i < dim1; i++)
        vector[i] = 0.0;
}

/* Initialize a vector of integers with zeros */
void init_int_vec(int *vector, int dim)
{
    int i;
    
    for (i=0; i<dim; i++)
        vector[i] = 0;
}

/* Initialize a matrix of doubles with zeros */
void init_mat(double **matrix, int dim1, int dim2)
{
    int i, j;
    
    for (i=0; i<dim1; i++)
        for(j=0; j<dim2; j++)
            matrix[i][j] = 0.0;
}

/* Initialize a matrix of integers with zeros */
void init_int_mat(int **matrix, int dim1, int dim2)
{
    int i, j;
    
    for (i=0; i<dim1; i++)
        for(j=0; j<dim2; j++)
            matrix[i][j] = 0;
}


/* Initialize a matrix of integers with zeros */
void init_int_mat3(int ***matrix, int dim1, int dim2,int dim3)
{
    int i,j,k;
    
    for (i=0; i<dim1; i++)
        for(j=0; j<dim2; j++)
            for(k=0; k<dim3; k++)
                matrix[i][j][k] = 0;
}

/*Initialize an array of SparseDim*/
void init_SparseDims(SparseDim* sdim,int dim1){
    int i;
    
    for(i=0;i<dim1;i++){
        sdim[i].nbElems   = 0;
        sdim[i].starting  = 0;
        sdim[i].absorbing = 0;
        sdim[i].sum       = 0;
        sdim[i].first     = NULL;
        sdim[i].last      = NULL;
        sdim[i].firstAbs  = NULL;
        sdim[i].lastAbs   = NULL;        
    }    
}    

/****************************I/O OPERATIONS*******************************/

/*Tell if a file contains sparse data*/
int isSparseFile(char* fname){
    FILE* inFile = fopen(fname,"r");
    char  line[BUF_LINE_LEN];
    int   limit=BUF_LINE_LEN-1;
    char* sym;
    
    if(!inFile){
        printf("File %s does not exist\n",fname);
        exit(EXIT_FAILURE);
    }    
    
    /*Read the first line*/
    fgets(line,BUF_LINE_LEN,inFile);
    check_fgets(line,limit,1,fname);
    sym = (char*)strtok(line," \n");
    sym = (char*)strtok(NULL, " \n");
    
    fclose(inFile);
    
    /*Several elements on the first line --> dense format*/
    return (sym) ? 0 : 1;
}  

/* Read a vector of double in input file */
int readVec(char* fname,double* I,int m){
    
    FILE* inFile = fopen(fname,"r");
    int   i,nnz=0;
    float val;
    
    if(!inFile){
        printf("File %s does not exist\n",fname);
        exit(EXIT_FAILURE);
    }       
    
    for(i=0;i<m;i++){
        fscanf(inFile,"%e",&val);
        I[i] = (double)val;
        if (val != 0) nnz++;
    }
    fclose(inFile);
    return nnz;
}

/* Read matrix of double in input file */
int readMat(char* fname,double** P,int m){
    
    FILE* inFile = fopen(fname,"r");
    int   i,j,nnz=0;
    float val;
    
    if(!inFile){
        printf("File %s does not exist\n",fname);
        exit(EXIT_FAILURE);
    }    
    
    for(i=0;i<m;i++){
        for(j=0;j<m;j++){
            fscanf(inFile,"%e",&val);
            P[i][j] = (double)val;
            if (val != 0) nnz++;
        }
    }
    
    fclose(inFile);
    return nnz;
}

/* Read matrix of double in input file */
int readMat_sparse(char* fname,SparseDim* rows,SparseDim* cols,int m,DegStat* degStat){
    FILE*       inFile = fopen(fname,"r");
    char        line[BUF_LINE_LEN];
    int         limit=BUF_LINE_LEN-1;
    char*       sym;
    int         from,to;
    int         nnz=0;
    int         sumDeg=0;
    float       val;
    char        *ctxL,*ctxE;
    int         iline=0;
        
    if(!inFile){
        printf("File %s does not exist\n",fname);
        exit(EXIT_FAILURE);
    }    
    
    /*Initialize the degree statistics*/
    degStat->minDeg  = m;
    degStat->maxDeg  = 0;
    degStat->meanDeg = 0.0;
    
    /**SPARSE FILE**/
    if(isSparseFile(fname)){
        /*Skip the first line*/
        fgets(line,BUF_LINE_LEN,inFile);
	check_fgets(line,limit,iline++,fname);		/* Check buffer overflow */
        while(fgets(line,BUF_LINE_LEN,inFile)){
	  check_fgets(line,limit,iline++,fname);	/* Check buffer overflow */
	  /*Fetch the "from" node number*/
	  if ((sym=(char*)strtok_r(line," \n",&ctxL))){
	    from  = atoi(sym)-1;
	    /*Scan all the "to" nodes*/
	    while((sym=(char*)strtok_r(NULL, " \n",&ctxL))){
	      to  = atoi((char*)strtok_r(sym,":",&ctxE))-1;
	      val = atof((char*)strtok_r(NULL,":",&ctxE));
	      /*Add the element in the row "from"*/
	      if (!rows[from].last){
		rows[from].first = malloc(sizeof(SparseElem));
		rows[from].last  = rows[from].first;
	      }
	      else{
		rows[from].last->nextC  = malloc(sizeof(SparseElem));
		rows[from].last         = rows[from].last->nextC;
	      }
	      /*Put the value in the element*/
	      rows[from].last->l     = from;
	      rows[from].last->c     = to;
	      rows[from].last->val   = val;
	      rows[from].last->ept   = 0;
	      rows[from].last->cur   = 0;
	      rows[from].last->diff  = 0;
	      rows[from].last->nextL = NULL;
	      rows[from].last->nextC = NULL;
	      rows[from].nbElems    += 1;
	      rows[from].sum        += val;
                    
	      /*Add the element in the column "to"*/
	      if(!cols[to].last){
		cols[to].first  = rows[from].last;
		cols[to].last   = rows[from].last;
	      }
	      else{
		cols[to].last->nextL = rows[from].last;
		cols[to].last = rows[from].last;
	      }
	      cols[to].nbElems += 1;
	      cols[to].sum += val;
	      /*Update the number of non-nul elements*/
	      nnz += 1;                    
	    }/*End while(to node)*/
	    if(rows[from].nbElems > degStat->maxDeg)
	      degStat->maxDeg = rows[from].nbElems;
	    if(rows[from].nbElems < degStat->minDeg)
	      degStat->minDeg = rows[from].nbElems;                
	    sumDeg += rows[from].nbElems;                
	  }/*End if(line contains tokens)*/    
        }/*End iterate on the lines*/    
       degStat->meanDeg = (float)sumDeg/(float)m;
    }/*End if(sparse)*/
    /**DENSE FILE**/
    else{
        for(from=0;from<m;from++){
            for(to=0;to<m;to++){
                fscanf(inFile,"%e",&val);    
                if(val != 0){
                    /*Add the element in the row "from"*/
                    if (!rows[from].last){
                        rows[from].first = malloc(sizeof(SparseElem));
                        rows[from].last  = rows[from].first;
                    }
                    else{
                        rows[from].last->nextC  = malloc(sizeof(SparseElem));
                        rows[from].last         = rows[from].last->nextC;
                    }
                    /*Put the value in the element*/
                    rows[from].last->l     = from;
                    rows[from].last->c     = to;
                    rows[from].last->val   = val;
                    rows[from].last->ept   = 0;
                    rows[from].last->cur   = 0;
                    rows[from].last->diff  = 0;                    
                    rows[from].last->nextL = NULL;
                    rows[from].last->nextC = NULL;
                    rows[from].nbElems    += 1;
                    rows[from].sum        += val;
                    /*Add the element in the column "to"*/
                    if(!cols[to].last){
                        cols[to].first  = rows[from].last;
                        cols[to].last   = rows[from].last;
                    }
                    else{
                        cols[to].last->nextL = rows[from].last;
                        cols[to].last = rows[from].last;
                    }
                    cols[to].nbElems += 1;
                    cols[to].sum     += val;
                    /*Update the number of non-nul elements*/
                    nnz += 1;
                } 
            }
            if(rows[from].nbElems > degStat->maxDeg)
                degStat->maxDeg = rows[from].nbElems;
            if(rows[from].nbElems < degStat->minDeg)
                degStat->minDeg = rows[from].nbElems;
            sumDeg += rows[from].nbElems;
        }
        degStat->meanDeg = (float)sumDeg/(float)m;
    }
    fclose(inFile);
    return nnz;
}      

/* Read the number of nodes from the adjacency matrix file */
int readNbNodes(char* fname){
    FILE*  inFile  = fopen(fname,"r");
    char   line[BUF_LINE_LEN];
    int    limit=BUF_LINE_LEN-1;
    char*  sym;
    int    j,nbNodes=0;
    
    if(!inFile){
        printf("File %s does not exist\n",fname);
        exit(EXIT_FAILURE);
    }    
    
    if(isSparseFile(fname)){
        fscanf(inFile,"%i",&nbNodes);
    }
    else{
        /*Read the first line*/
        fgets(line,BUF_LINE_LEN,inFile);
	check_fgets(line,limit,1,fname);	/* Check buffer overflow */
        sym = (char*)strtok(line," \n");
        if (sym != NULL){
            /*Scan all the nodes*/
            for(j=0;sym != NULL;++j){
                nbNodes++;
                sym = (char*)strtok(NULL, " \n");
            }
        }
    }

    fclose(inFile);
    return nbNodes;
}


/* Write a vector of double in a file*/
void writeVec(double* I,int m,char* fname){
    
    int i;
    FILE* outFile = fopen(fname,"w");
     
    for(i=0;i<m;i++){
        fprintf(outFile,"%e\n",I[i]);
    }
    fclose(outFile);    
}

/* Write a vector of double in a sparse file*/
void writeVec_sparse(double* I,int m,char* fname){
    
    int i;
    FILE* outFile = fopen(fname,"w");
     
    for(i=0;i<m;i++){
        if(I[i] > 0)
            fprintf(outFile,"%i:%e\n",i+1,I[i]);
    }
    fclose(outFile);    
}

/* Write a matrix of double in a file*/
void writeMat(double** P,int m,char* fname){
    
    int i,j; 
    FILE* outFile = fopen(fname,"w");
    
    for(i=0;i<m;i++){
        for(j=0;j<m;j++){
            fprintf(outFile,"%e ",P[i][j]);
        }
        fprintf(outFile,"\n");
    }
    fclose(outFile);    
}

/*Write the MEPT in a sparse file*/
void writeMEPT_sparse(char* fname,SparseDim* rows,int m){
    
    int i,add;
    FILE* outFile = fopen(fname,"w");
    SparseElem* elem;
    
    /*Write the number of lines (i.e. nodes)*/
    fprintf(outFile,"%i\n",m);
    
    /*Write the sparse matrix entries*/
    for(i=0;i<m;i++){
        /*Check if the row contains at least one non-zero entry*/
        for(add=0,elem=rows[i].first;!add&&elem;elem=elem->nextC){
            if(elem->ept>0){
                add=1;
            }
        }
        if(add){           
            /*Write the line number*/
            fprintf(outFile,"%i",i+1);
            elem = rows[i].first;
            while(elem){
                if(elem->ept > 0)
                    fprintf(outFile," %i:%e",(elem->c)+1,elem->ept);
                elem = elem->nextC;
            }
            fprintf(outFile,"\n");
        }
    }        
    fclose(outFile);    
}    
  
/*Write the absolute MEPT difference in a sparse file*/
void writeDiff_sparse(char* fname,SparseDim* rows,int m){
    
    int i,add;
    FILE* outFile = fopen(fname,"w");
    SparseElem* elem;
    
    /*Write the number of lines (i.e. nodes)*/
    fprintf(outFile,"%i\n",m);
    
    /*Write the sparse matrix entries*/
    for(i=0;i<m;i++){
        /*Check if the row contains at least one non-zero entry*/
        for(add=0,elem=rows[i].first;!add&&elem;elem=elem->nextC){
            if(elem->diff>0){
                add=1;
            }
        }
        if(add){           
            /*Write the line number*/
            fprintf(outFile,"%i",i+1);
            elem = rows[i].first;
            while(elem){
                if(elem->diff > 0)
                    fprintf(outFile," %i:%e",(elem->c)+1,elem->diff);
                elem = elem->nextC;
            }
            fprintf(outFile,"\n");
        }
    }        
    fclose(outFile);    
}  

/*Export a matrix of MEPN to the graphviz dot format*/
void writeDot(char* fname,SparseDim* rows,int m,int undirected){
    
    int i;
    SparseElem* elem;
    FILE* outFile = fopen(fname,"w");
    
    if(undirected){
        /*Write header*/
        fprintf(outFile,"graph MEPN {\n");
        fprintf(outFile,"rankdir=LR;\n");
        fprintf(outFile,"size=\"8,5\";\n");
        fprintf(outFile,"orientation=land;\n");
        fprintf(outFile,"node [shape = circle,fontcolor=black];\n");
        /*Define edges*/
        for(i=0;i<m;i++){
            elem = rows[i].first;
            while(elem && elem->c <= i){
                if(elem->diff > 0){
                    fprintf(outFile,"%i -- %i [label = \"%e\"];\n",i+1,(elem->c)+1,elem->diff);
                }
                elem = elem->nextC;
            }
        }
        /*Write trailer*/
        fprintf(outFile,"}\n");
    }
    else{
        /*Write header*/
        fprintf(outFile,"digraph MEPN {\n");
        fprintf(outFile,"rankdir=LR;\n");
        fprintf(outFile,"size=\"8,5\";\n");
        fprintf(outFile,"orientation=land;\n");
        fprintf(outFile,"node [shape = circle,fontcolor=black];\n");
        /*Define edges*/
        for(i=0;i<m;i++){
            elem = rows[i].first;
            while(elem){
                if(elem->ept > 0)
                    fprintf(outFile,"%i -> %i [label = \"%e\"];\n",i+1,(elem->c)+1,elem->ept);
                elem = elem->nextC;
            }
        }
        /*Write trailer*/
        fprintf(outFile,"}\n");
    }
    fclose(outFile);
}

/****************************Miscellaneous********************************/

/* Print command line help for the limited kwalks*/
void printHelpKwalk(){
    printf("Usage     : lkwalk -g graphFile -l walkLen -k relNodes -o outFile [options]\n\n");
    printf("Function  : Compute the limited K-walk Expected Passage Times\n");
    printf("Arguments : -g graphFile  the file containing the (sparse) adjacency matrix\n");
    printf("            -l walkLen    the length of the limited k-walk\n");
    printf("            -k relNodes   the list of nodes of interest\n");
    printf("                          groups are separated by # and nodes are separated by :\n");
    printf("            -o outFile    the output file containing the expected times\n");
    printf("Options   : -u            consider walks up to length walkLen\n");
    printf("            -p relProbas  initial probas for nodes of interest\n");
    printf("                          groups are separated by # and nodes are separated by :\n");
	printf("            -i nbInflat   the number of inflation iterations (default 0)\n");
    printf("            -d            generate a dot file (outFile.dot)\n");    
    printf("            -v verbose    verbosity level (default 1)\n");
    printf("            -h            display this help\n\n");
    printf("Version   : %1.2f\n",VERSION);
    printf("Author    : J. Callut and P. Dupont (2010)\n\n");
}

/* Display a vector of integers */
void show_int_vec(int* vec,int n){
    int i;
    for(i=0;i<n;++i){
        printf("%i ",vec[i]);
    }
    printf("\n");
}

/* Display a vector of doubles */
void show_vec(double* vec,int n){
    int i;
    for(i=0;i<n;++i){
        printf("%e ",vec[i]);
    }
    printf("\n");
}

/* Display a matrix of doubles */
void show_int_mat(int** mat,int l,int c){
    int i,j;
    for(i=0;i<l;++i){
        for(j=0;j<c;++j){
            printf("%i ",mat[i][j]);
        }
        printf("\n");
    }
}

/* Display a matrix of doubles */
void show_mat(double** mat,int l,int c){
    int i,j;
    for(i=0;i<l;++i){
        for(j=0;j<c;++j){
            printf("%e ",mat[i][j]);
        }
        printf("\n");
    }
}

/* Display an array of SparseDim */
void show_SparseDims(SparseDim* sdim,int m,char mode){
    int i;
    SparseElem* elem;
    printf("\n");
    for (i=0;i<m;i++){
        elem = sdim[i].first;
        while(elem){
            if(mode == 'l'){
                printf("%i:%e ",(elem->c)+1,elem->val);
                elem = elem->nextC;                
                }
            else if(mode == 'c'){
                printf("%i:%e ",(elem->l)+1,elem->val);
                elem = elem->nextL;
                }
        }
        printf("\n");
    }
}  

/*Tokenize a list of filenames separated by commas*/
int tokenizeFnames(char* list,char** fnames){
    
    int   nbFiles = 0;
    char* fname   = (char*)strtok(list,",");
    
    while (fname != NULL){
        fnames[nbFiles++] = fname;
        fname = (char*)strtok(NULL,",");
    }
    return nbFiles;
}

/*Add an extension to a file name*/
void addExtension(char* fname,char* extension,char* target){
    if(snprintf(target,FNAME_LEN,"%s.%s",fname,extension) >= FNAME_LEN){
        printf("Argument outFile is too long\n");
        exit(EXIT_FAILURE);
    }       
}    


/*Tokenize a list of integers separated by :*/
int tokenizeInt(char* ilist,int* olist){
    
    int   size = 0;
    char* elem = (char*)strtok(ilist,":");
    
    while (elem != NULL){
        olist[size++] = atoi(elem);
        elem = (char*)strtok(NULL,":");
    }
    return size;
}


/*Tokenize a list of doubles separated by :*/
int tokenizeDouble(char* ilist,double* olist){
    
    int   size = 0;
    char* elem = (char*)strtok(ilist,":");
    
    while (elem != NULL){
        olist[size++] = atof(elem);
        elem = (char*)strtok(NULL,":");
    }
    return size;
}

/*Return the number of groups and the max number of nodes per group*/
void getGroupsInfo(char* list,int *nbGroups, int *maxNPG){
 
    int   nbNodes=0;
    char  *elemG,*elemN;
    char  *ctxG=NULL;
    char  *ctxN;
    
    /*Initialize output arguments*/
    *nbGroups = 0;
    *maxNPG   = 0;
    
    /*Tokenize on the group delimiter # */
    elemG = (char*)strtok_r(list,"#",&ctxG); 

    while (elemG != NULL){
        /*Tokenize on the node delimiter : */
        elemN   = (char*)strtok_r(elemG,":",&ctxN);
        nbNodes = 0;
        while(elemN != NULL){
            elemN = (char*)strtok_r(NULL,":",&ctxN);
            nbNodes++;
        }
        /*Update the max number of nodes per group*/
        if(nbNodes > *maxNPG){
            *maxNPG = nbNodes;
        }
        (*nbGroups)++;
        elemG = (char*)strtok_r(NULL,"#",&ctxG);
    }       
}    

/*Tokenize groups of relevant nodes*/
int tokenizeKGroups(char* list,int** kgroups){
 
    int   nbGroups=0,nbNodes=0;
    char  *elemG,*elemN;
    char  *ctxG=NULL;
    char  *ctxN;
    
    /*Tokenize on the group delimiter # */

    elemG = (char*)strtok_r(list,"#",&ctxG);
    
    while (elemG != NULL){
        /*Tokenize on the node delimiter : */
        elemN   = (char*)strtok_r(elemG,":",&ctxN);
        nbNodes = 0;
        while(elemN != NULL){
            kgroups[nbGroups][++nbNodes] = atoi(elemN)-1;
            elemN = (char*)strtok_r(NULL,":",&ctxN);
        }
        kgroups[nbGroups++][0] = nbNodes;
        elemG = (char*)strtok_r(NULL,"#",&ctxG);
    }
    return nbGroups;
}    

/*Tokenize initial probabilities of groups*/
int tokenizeKGProbas(char* list,double** kgprobas){
 
    int   nbGroups=0,nbNodes=0;
    char  *elemG,*elemN;
    char  *ctxG=NULL;
    char  *ctxN;
    
    /*Tokenize on the group delimiter # */
    elemG = (char*)strtok_r(list,"#",&ctxG);
    
    while (elemG != NULL){
        /*Tokenize on the node delimiter : */
        elemN   = (char*)strtok_r(elemG,":",&ctxN);
        nbNodes = 0;
        while(elemN != NULL){
            kgprobas[nbGroups][++nbNodes] = atof(elemN);
            elemN = (char*)strtok_r(NULL,":",&ctxN);
        }
        kgprobas[nbGroups++][0] = nbNodes;
        elemG = (char*)strtok_r(NULL,"#",&ctxG);
    }
    return nbGroups;
}  

/* Check whether the maximal size of the reading buffer is reached */
void check_fgets(char *s, int limit, int iline, char *filename) {
  if(strlen(s)==limit) {
    fprintf(stderr,"Line buffer overflow at line %d of \"%s\"\n",iline,filename);
    fprintf(stderr,"You should probably increase BUF_LINE_LEN (=%d) declared in tools.h\n",BUF_LINE_LEN);
    exit(EXIT_FAILURE);
  }
}
/****************************Matrix operations****************************/


/* Set a vector with element i equal 1 and others equal 0 */
void setE_i(double* vec,int m,int i){
    int j;
    
    for(j=0;j<m;++j){
        vec[j] = 0.0;
    }
    vec[i] = 1.0;
}

/* Compute the sum of a vector of doubles*/
double vecSum(double* vec,int m){
    int i;
    double sum =0.0;
    
    for(i=0;i<m;i++){
        sum += vec[i];
    }
    return sum;    
}

/* Compute the sum of a vector of doubles*/
int int_vecSum(int* vec,int m){
    int i;
    int sum = 0;
    
    for(i=0;i<m;i++){
        sum += vec[i];
    }
    return sum;    
}

/*Divide all elements of a vector of doubles by a number a*/
void vecNormalize(double* vec,int m,double a){    
    int i;
    
    for(i=0;i<m;i++)
        vec[i] /= a;
}    

/*Divide all elements of a matrix of doubles by a number a*/
void matNormalize(double** P,int m,double a){
    int i,j;
    
    for(i=0;i<m;i++)
        for(j=0;j<m;j++)
            P[i][j] /= a;    
}    

/*Divide all elements of a sparse matrix by a number a*/
void matNormalize_sparse(SparseDim* rows,int m,double a){
    
    int i;
    SparseElem* elem;
      
    /*Normalize the entries*/
    for(i=0;i<m;i++){
        elem = rows[i].first;
        while(elem){
            elem->ept /= a;
            elem = elem->nextC;
        }
    }        
}  

/*Compute the position of the maximum element in a vector of doubles*/
int max_vec_pos(double* vec,int m){
    int i, maxPos=0;
    double maxVal=vec[0];
    
    for(i=0;i<m;++i){
        if(vec[i] > maxVal){
            maxVal = vec[i];
            maxPos = i;
        }    
    }
    return maxPos;
}

/*Enforce the stochasticity of a matrix of positive double */
void stochMat(double** P,int m){
    int    i,j;
    double sum;
    
    for(i=0;i<m;++i){
        sum = vecSum(P[i],m);
        if(sum > 0){
            for(j=0;j<m;++j){
                P[i][j] /= sum;
            }
        }
    }            
}

/*Enforce the stochasticity of a sparse matrix of positive double*/
void stochMat_sparse(SparseDim* rows,int m){
    int i;
    SparseElem* elem;
    
    for(i=0;i<m;i++){
        elem = rows[i].first;
        while(elem){
            elem->val /= rows[i].sum;
            elem = elem->nextC;
        }    
    }    
}    


/*Enforce the stochasticity of a vector of positive double */
double stochVector(double* vec,int m){
    int i;    
    double sum = vecSum(vec,m);
    
    if (sum>0){
        for(i=0;i<m;++i)
            vec[i] = vec[i]/sum;
    }
    return sum;
}    


/* Build a forward lattice up to wlen*/
void buildLatticeForward(double** lat,int wlen,int* kgroup,double* kgproba,double** P,int m){
    
    int i,j,k;
    
    /* Clean the lattice */
    init_mat(lat,m,wlen+1);
    
    /* Initialization for time 0 */
    for(i=1;i<=kgroup[0];i++){
        lat[kgroup[i]][0] = kgproba[i];
    }
    
    /* Propagate probabilities in the lattice */
    for (k=1;k<=wlen;++k){
        for (i=0;i<m;i++){
            for (j=0;j<m;j++){
                if(P[j][j] != 1 && P[j][i] > 0)
                    lat[i][k] += lat[j][k-1]*P[j][i];
            }
        }
    }
}

/* Build a forward lattice up to wlen*/
void buildLatticeForward_sparse(double** lat,int wlen,int* kgroup,double* kgproba,SparseDim* cols,int m,int verbose){
    
    int j,k;
    SparseElem* elem;
    
    /* Clean the lattice */
    init_mat(lat,m,wlen+1);
    
    /* Initialization for time 0 */
    for(j=1;j<=kgroup[0];j++){
        lat[kgroup[j]][0] = kgproba[j];
    }
    
    /* Propagate probabilities in the lattice */
    if(verbose >= 2)
        printf("Forward lattice  : 0");
    for (k=1;k<=wlen;++k){
        if(verbose >= 2 && k%10 == 0){
            printf("...%i",k);
            fflush(stdout);
        }
        for (j=0;j<m;j++){
            elem = cols[j].first;
            while(elem){
                if(!cols[elem->l].absorbing)
                    lat[j][k] += lat[elem->l][k-1]*(elem->val);
                elem = elem->nextL;
            }    
        }
    }
    if(verbose >= 2)
        printf("\n");
}

/* Build a backward lattice up to wlen*/
void buildLatticeBackward(double** lat,int wlen,double** P,int m){
    
    int i,j,k;
    
    /* Clean the lattice */
    init_mat(lat,m,wlen+1);
    
    
    /* Initialization for time wlen */
    for(i=0;i<m;i++){
        if(P[i][i] == 1)
            lat[i][wlen] = 1;
    }    
       
    /* Propagate probabilities in the lattice */
    for (k=wlen;k>=1;--k){ 
        for(j=0;j<m;j++){
            for (i=0;i<m;i++){
                if(P[i][i] != 1 && P[i][j] > 0)
                    lat[i][k-1] += P[i][j]*lat[j][k];
            }
        }
    }
}

/* Build a backward lattice up to wlen*/
void buildLatticeBackward_sparse(double** lat,int wlen,SparseDim* rows,int m,int verbose){
    
    int i,k;
    SparseElem* elem;
    
    /* Clean the lattice */
    init_mat(lat,m,wlen+1);
    
    /* Initialization for time wlen */
    for(i=0;i<m;i++){
        if(rows[i].absorbing)
            lat[i][wlen] = 1;
    }    
       
    /* Propagate probabilities in the lattice */
    if(verbose >= 2)
        printf("Backward lattice : ");
    for (k=wlen;k>=1;--k){
        if(verbose >= 2 && k%10 == 0){
            printf("%i...",k);
            fflush(stdout);
        }               
        for(i=0;i<m;i++){
            if(!rows[i].absorbing){
                elem = rows[i].first;
                while(elem){
                    lat[i][k-1] += lat[elem->c][k]*(elem->val);
                    elem = elem->nextC;
                }    
            }
        }
    }
    if(verbose >= 2)
        printf("0\n");
}

/*Compute the cumulated absorption probability*/
double cumulatedAbsorptionProba(double** latB,int wlen,int* kgroup,double* kgproba){
	
	int    j,k;
	double mass = 0.0;
	
	for (k=wlen-1;k>=0;--k){
		for(j=1;j<=kgroup[0];j++){
			mass += kgproba[j]*latB[kgroup[j]][k];
	    }
	}
	
	return mass;
}

/*Set all the groups being absorbing except the starting group*/
void setAbsStates(double** P,double** PStart,int m,int** kgroups,int nbGroups,int start){
    
    int i,j,k;
    int k_ij;
    
    /*Copy P -> PStart*/
    for (i=0;i<m;i++){
        for(j=0;j<m;j++){
            PStart[i][j] = P[i][j];
        }
    }
    
    /*Set all the groups being absorbing except the starting group*/
    for(i=0;i<nbGroups;i++){
        if(i != start){
            for(j=1;j<=kgroups[i][0];j++){
                k_ij = kgroups[i][j];
                for(k=0;k<m;k++){
                    PStart[k_ij][k] = (k == k_ij) ? 1 : 0;
                }
            }
        }
    }
}

/*Set all the groups being absorbing except the starting group*/
void setAbsStates_sparse(SparseDim* rows,SparseDim* cols,int m,int** kgroups,double** kgprobas,int nbGroups,int start){
    
    int i,j;
    SparseElem* selem;
    
    /*Reset the absorbing and starting states*/
    for(i=0;i<m;i++){
        rows[i].absorbing  = 0;
        cols[i].absorbing  = 0;
        rows[i].starting   = 0;
        cols[i].starting   = 0;
        rows[i].firstAbs   = NULL; /*Controlled memory leak*/
        rows[i].lastAbs    = NULL; /*Controlled memory leak*/       
    }
    
    /*Set all the groups being absorbing except the starting group*/
    /*Set the starting group*/
    for(i=0;i<nbGroups;i++){
        for(j=1;j<=kgroups[i][0];j++){
            if(i != start){
                rows[kgroups[i][j]].absorbing = 1;
                cols[kgroups[i][j]].absorbing = 1;
                selem = cols[kgroups[i][j]].first;
                while(selem){
                    if(selem->l != selem->c){
                        if(!rows[selem->l].firstAbs){
                            rows[selem->l].firstAbs = malloc(sizeof(Elem));
                            rows[selem->l].lastAbs  = rows[selem->l].firstAbs;
                        }
                        else{
                            rows[selem->l].lastAbs->next = malloc(sizeof(Elem));
                            rows[selem->l].lastAbs = rows[selem->l].lastAbs->next;
                        }
                        rows[selem->l].lastAbs->selem = selem;
                        rows[selem->l].lastAbs->next = NULL;
                    }
                    selem = selem->nextL;
                }                    
            }
            else{
                rows[kgroups[i][j]].starting = kgprobas[i][j];
                cols[kgroups[i][j]].starting = kgprobas[i][j];                
            }    
        }
    }
}

/*Compute the absolute EPT difference : |E_ij - E_ji|*/
void diffEPT(double** E,int m){
    
    int i,j;
    
    for(i=0;i<m;i++){
        E[i][i] = 0;
        for(j=0;j<i;j++){
            E[i][j] = ABS(E[i][j] - E[j][i]);
            E[j][i] = E[i][j];
        }    
    }        
}    

/*Compute the absolute EPT difference : |E_ij - E_ji|*/
void diffEPT_sparse(SparseDim* rows,SparseDim* cols,int m){
    
    int i;
    double diff;
    SparseElem* elemL;
    SparseElem* elemC;
    
    for(i=0;i<m;i++){        
        elemL = rows[i].first;
        elemC = cols[i].first;        
        /*Scan the two list in parallel*/
        while(elemL && elemC){
            if(elemL->c >= i){            
                diff = ABS((elemL->cur) - (elemC->cur));
                elemL->diff += diff;
                elemC->diff += diff;
                elemL->cur   = 0;
                elemC->cur   = 0;
            }
            elemL = elemL->nextC;
            elemC = elemC->nextL;            
        }
    }
}

/*Test if a matrix is symmetric*/
int isSymmetric(double** P,int m){
    int i,j;
    int sym = 1;
    
    for(i=0;i<m && sym ;i++){
        for(j=0;j<m && sym;j++){
            if (P[i][j] != P[j][i])
                sym = 0;
        }    
    }  
    return sym;
}

/*Test if a sparse matrix is symmetric*/
int isSymmetric_sparse(SparseDim* rows,SparseDim* cols,int m){
    int i,sym=1;
    SparseElem* elemL;
    SparseElem* elemC;
    
    for(i=0;sym && i<m;i++){
        if(rows[i].nbElems == cols[i].nbElems){
            elemL = rows[i].first;
            elemC = cols[i].first;           
            
            /*Scan the two list in parallel*/
            while(sym && elemL && elemC){                
                if(elemL->c != elemC->l){
                    sym = 0;
                    }
                else if(elemL->val != elemC->val){
                    sym = 0;
                    }
                else{
                    elemL = elemL->nextC;
                    elemC = elemC->nextL;
                }
            }
        }
        else
            sym = 0;
    }
    return sym;
}
  

/*Get the number of non-nul elements in a matrix*/
int matrixNNZ(double** P,int m){
    int i,j;
    int nnz = 0;

    for(i=0;i<m;i++){
        for(j=0;j<m;j++){
            if (P[i][j] != 0)
                nnz++;
        }    
    }
    return nnz;    
}    

/*Get the number of non-nul elements in a matrix*/
void matrixNNZ_sparse(SparseDim* rows,int m,int* nnz,int* nnr){
    int i;
    int *touched = malloc(m*sizeof(int));
    SparseElem* elem;
    *nnz=0,*nnr=0;
    
    init_int_vec(touched,m);    
    for(i=0;i<m;i++){
        elem = rows[i].first;
        while(elem){
            if(elem->ept > 0){
                if (!touched[i]){
                    touched[i]=1;
                    (*nnr)++;                    
                }    
                if (!touched[elem->c]){
                    touched[elem->c]=1;
                    (*nnr)++;
                }    
                (*nnz)++;
            }    
            elem = elem->nextC;
        }    
    }
    free(touched);
}    
 

/*Copy edge passage times to edge value for inflation*/
void copyForInflation(SparseDim* rows,int undirected,int m){
    
    int i;
    SparseElem* elem;

	for(i=0;i<m;i++){
		elem = rows[i].first;
		rows[i].sum = 0.0; /*The col sum should also be updated*/
		while(elem){
			/*Copy edge passage times to edge value*/
			if(undirected)
				elem->val = elem->diff;
			else
				elem->val = elem->ept;
			rows[i].sum += elem->val;
			elem->ept  = 0.0;
			elem->cur  = 0.0;
			elem->diff = 0.0;
	 		elem = elem->nextC;
		}	
	}
}

/*Copy a vetcor into another (must be allocated)*/
void copy_mat(double** src,double** dest,int l,int c){
	int i,j;
	if(src != NULL && dest != NULL){
		for(i=0;i<l;i++){
			for(j=0;j<c;j++)
				dest[i][j] = src[i][j]; 
		}
	}
}

/*************************************************************************/
