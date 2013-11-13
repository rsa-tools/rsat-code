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
/* Modified    : 09/10/2010 by P. Dupont                                 */
/*               - Cleanup code + fix in tools.c                         */
/* License     : GPL 3.0                                                 */ 
/*               GNU GENERAL PUBLIC LICENSE                              */
/*               Version 3, 29 June 2007                                 */
/*               For details, see the file gpl-3.0.txt in this bundle    */
/*************************************************************************/

#include "tools.h"

/*Function computing the limited k-walks*/
int lkwalk (int argc, char* argv[]);

int lkwalk (int argc, char* argv[]){
    
  /*getopt() variables */
    int opt;
     
  /*Command line options*/
    int       gflg=0,lflg=0,kflg=0,oflg=0;
    int       pflg=0,iflg=0,vflg=0;
    char*     optarg_dup;
    char*     graphFname=NULL;
    char*     outFname;
    char      outEPT[FNAME_LEN],outNPT[FNAME_LEN];
    char      outEPT_Norm[FNAME_LEN],outNPT_Norm[FNAME_LEN];
    char      outDIF[FNAME_LEN],outDOT[FNAME_LEN];
    int**     kgroups=NULL;
    double**  kgprobas=NULL;
    double**  kgprobas_stoch=NULL;
    int       maxNPG;
    int       nbGroups=0;
    int       nbNodes;
    int       nbEdges;
    int       nbKnodes;
    int       wlen=0;
    int       nbInflat = 0;
    int       upto     = 0;
    int       genDot   = 0;
    int       verbose  = 2;
    
  /*Working variables */
    int         i,j,k,t,l;
    int         undirected;
    int         start_cpu,stop_cpu;
    int         nnz,nnr;
	int			inflatIter;
    double      p_wlen;
    double      tmp,accF,accB;
    double      update;
    double      mass_i   = 0;
    double      mass_abs = 0;
    double      expWlen  = 0;
    double*     N;
    double**    latF;
    double**    latB;
    SparseDim*  rows;
    SparseDim*  cols;
    SparseElem* selem;
    Elem*       elem;
    DegStat     degStat;

  /*Scan the command line*/
    while ((opt = getopt(argc,argv,"g:l:k:o:up:i:dv:h")) != EOF){
        
        switch(opt){
            case 'g':
                gflg++;
                graphFname = optarg;
                break;
            case 'l':
                lflg++;
                wlen = atoi(optarg);
                break;
            case 'k':
                kflg++;
                optarg_dup = strdup(optarg);
                getGroupsInfo(optarg,&k,&maxNPG);
		if(pflg && k != nbGroups){
                    printf("Initial probabilities mismatch !\n");
                    return EXIT_FAILURE;
                }
                else{
                    nbGroups = k;
                }
                kgroups = int_mat_alloc(nbGroups,maxNPG+1);
                tokenizeKGroups(optarg_dup,kgroups);
                break;
            case 'o':
                oflg++;
                outFname = optarg;
                addExtension(outFname,"E",outEPT);
                addExtension(outFname,"N",outNPT);
                addExtension(outFname,"E.norm",outEPT_Norm);
                addExtension(outFname,"N.norm",outNPT_Norm);
                addExtension(outFname,"dif",outDIF);
                addExtension(outFname,"dot",outDOT);
                break;
            case 'u':
                upto = 1;
                break;
            case 'p':
                pflg++;
                optarg_dup = (char*)strdup(optarg);
                getGroupsInfo(optarg,&k,&maxNPG);
                if(kflg && k != nbGroups){
                    printf("Initial probabilities mismatch !\n");
                    return EXIT_FAILURE;
                }
                else{
                    nbGroups = k;
                }
                kgprobas = mat_alloc(nbGroups,maxNPG+1);
				kgprobas_stoch = mat_alloc(nbGroups,maxNPG+1);
                tokenizeKGProbas(optarg_dup,kgprobas);
                break;

			case 'i':
				iflg++;
				nbInflat = atoi(optarg);
				break;
            case 'd':
                genDot = 1;
                break;
            case 'v':
                vflg++;
                verbose = atoi(optarg);
                break;
            case 'h':
                printHelpKwalk();
                return EXIT_FAILURE;
                break;
        }
    }    
    
   
    /*Check mandatory arguments*/
    if(!gflg || !lflg || !kflg){
        printf("Mandatory argument(s) missing\n");
        printHelpKwalk();
        return EXIT_FAILURE;
    }
    
    /*If a dot output is asked an output file has to be provided */
    if(genDot && !oflg){
        printf("Dot file can not be generated if no outFile is provided\n");
        return EXIT_FAILURE;
    }    
    
    /*Check if there is at least 2 relevant groups*/
    if(nbGroups < 2){
        printf("There must be at least two relevant nodes (in distinct groups)\n");
        return EXIT_FAILURE;
    }
    
    /*Check for the number of starting probabilities*/
    if(pflg){
        for(i=0;i<nbGroups;i++){
            if(kgroups[i][0] != (int)kgprobas[i][0]){
                printf("Initial probabilities mismatch !\n");
                return EXIT_FAILURE;
            }
        }
    }
    
    /*Compute the number of relevant nodes (in all groups)*/
    for(nbKnodes=0,i=0;i<nbGroups;i++){
        nbKnodes += kgroups[i][0];
    }
    
    /*If no starting probabilities are provided, set a uniform distribution*/
    if(!pflg){
        kgprobas = mat_alloc(nbGroups,maxNPG+1);
 		kgprobas_stoch = mat_alloc(nbGroups,maxNPG+1);
        /*Put the uniform distribution*/
        for(i=0;i<nbGroups;i++){
            kgprobas[i][0] = kgroups[i][0];
            for(j=1;j<=kgroups[i][0];j++){
                kgprobas[i][j] = 1.0/(float)nbKnodes;
            }
        }
    }
    /*Check if the provided probabilities sum up to 1*/
    else{
        tmp=0.0;
        for(i=0;i<nbGroups;i++){
            for(j=1;j<=kgprobas[i][0];j++){
                tmp += kgprobas[i][j];
            }
        }
        if(ABS(tmp-1.0) > PROB_EPS)
            printf("WARNING: The initial probabilities does not sum up to 1\n");
    }
    
    /*Get the number of nodes in the graph*/
    nbNodes = readNbNodes(graphFname);
    
    /*Check if the relevant nodes are out of range*/
    for(i=0;i<nbGroups;i++){
        for(j=1;j<=kgroups[i][0];j++){
            if(kgroups[i][j] < 0 || kgroups[i][j] >= nbNodes){
                printf("A node of interrest is out of range\n");
                return EXIT_FAILURE;
            }
        }
    }

    /*Allocate data structure*/
    N    = vec_alloc(nbNodes);
    latF = mat_alloc(nbNodes,wlen+1);
    latB = mat_alloc(nbNodes,wlen+1);
    rows = malloc(sizeof(SparseDim)*nbNodes);
    cols = malloc(sizeof(SparseDim)*nbNodes);
    
    /*Initialize the sparse transition matrix*/
    init_SparseDims(rows,nbNodes);
    init_SparseDims(cols,nbNodes);
    
    /*Read the adjacency matrix*/
    nbEdges    = readMat_sparse(graphFname,rows,cols,nbNodes,&degStat);    

    undirected = isSymmetric_sparse(rows,cols,nbNodes);

    /*Make it stochastic*/
    stochMat_sparse(rows,nbNodes);
    
    /*Print statistics about the graph*/
    if (verbose > 0){
        printf("\nNumber of nodes : %i\n",nbNodes);
        printf("Number of edges : %i\n",nbEdges);
        printf("Directed        : %s\n",(undirected)?"no":"yes");
        printf("Mean degree     : %2.2f\n",degStat.meanDeg);
        printf("Min  degree     : %i\n",degStat.minDeg);
        printf("Max  degree     : %i\n",degStat.maxDeg);
        printf("-------------------------------------------------------\n");
    }    

    /*Initialize passage times in nodes*/
    init_vec(N,nbNodes);

    /*Start the CPU chronometer*/
    start_cpu = clock();
    
	/*Loop on the inflation iteration*/
	for(inflatIter=0;inflatIter<=nbInflat;inflatIter++){
		
		/*Display the inflation number*/
		if(verbose > 0 && nbInflat > 0 && inflatIter > 0){
			printf("INFLATION ITERATION %i\n",inflatIter);
			printf("-------------------------------------------------------\n");
		}	

		/*Reset and copy data for the inflation iteration*/
		if(inflatIter > 0){
			undirected = 0;
			/*Reset variables*/
		    mass_abs = 0.0;
		    expWlen  = 0.0;
			/*Reset passage times in nodes*/
			init_vec(N,nbNodes);
			/*Copy the expectation passage times on the edges on the edge values*/
			copyForInflation(rows,undirected,nbNodes);
			/*Stochastize the matrix*/
			stochMat_sparse(rows,nbNodes);
		}

	    /*Loop on the relevant nodes*/
	    for(i=0;i<nbGroups;i++){
        
	        if(verbose > 0){
	            if(kgroups[i][0] == 1)
	                printf("Computing E[N(i,j) | Start in node %i]\n",kgroups[i][1]+1);
	            else    
	                printf("Computing E[N(i,j) | Start in group %i]\n",i+1);            
	        }    
                
	        /*Normalize the initial distribution*/
			copy_mat(kgprobas,kgprobas_stoch,nbGroups,maxNPG+1); /*We take a copy of kgprobas otherwise it would change through inflation interations*/
	        mass_i = stochVector(kgprobas_stoch[i]+1,kgprobas_stoch[i][0]);
        
	        /*Build the transition matrix w.r.t group i and set the starting states*/
	        setAbsStates_sparse(rows,cols,nbNodes,kgroups,kgprobas_stoch,nbGroups,i);
        
	        /*Build the forward/backward lattices up to length wlen*/
	        buildLatticeForward_sparse(latF,wlen,kgroups[i],kgprobas_stoch[i],cols,nbNodes,verbose);
	 		buildLatticeBackward_sparse(latB,wlen,rows,nbNodes,verbose); 
        
	        /*LINEAR UPTO MODE*/
	        if (upto){
		
				/*Compute the cumulated absorption mass */
				p_wlen = cumulatedAbsorptionProba(latB,wlen,kgroups[i],kgprobas_stoch[i]);
				mass_abs += mass_i*p_wlen;
		
	            for(k=0;k<nbNodes;k++){
	                accB = 0;
	                accF = 0;
	                for(l=wlen-1;l>=0;l--){
	                    if(!rows[k].absorbing){
	                        accB += latB[k][l];
	                        accF += latF[k][l];
	                        if(l > 0){
                            
	                           /************************/
	                           /*Transient to Transient*/
	                           /************************/
                            
	                            selem = cols[k].first;
	                            while (selem){
	                                if(!cols[selem->l].absorbing){
	                                    update = (mass_i*(selem->val)*latF[selem->l][l-1]*accB)/p_wlen;
	                                    selem->cur   += update;
	                                    selem->ept   += update;
	                                    expWlen      += update;
	                                    N[selem->l]  += update;
	                                }
	                                selem = selem->nextL;
	                            }
	                        }
	                        else{
                            
	                           /**************************/
	                           /*Absorbing from Transient*/
	                           /**************************/
                            
	                            elem = rows[k].firstAbs;
	                            while(elem){
	                                update = (mass_i*accF*(elem->selem->val))/p_wlen;
	                                elem->selem->cur  += update;
	                                elem->selem->ept  += update;
	                                expWlen           += update;
	                                N[k]              += update;
	                                elem = elem->next;
	                            }                            
	                        }
	                    }
	                }
	            }
	        }
	        /*EXACT LENGTH MODE*/
	        else{
	            /*Compute P[wlen]*/
	            for(p_wlen=0,l=1;l<=kgroups[i][0];l++){
	                p_wlen += kgprobas_stoch[i][l]*latB[kgroups[i][l]][0];
	            }
	            mass_abs += mass_i*p_wlen;
            
	            if(verbose > 1)
	                printf("P[Absorption | Start in group %i] = %e\n",i+1,p_wlen);
            
	            /*If there is a strictly positive proba, compute E[E(l,c)|wlen]*/
	            if(p_wlen > 0){
	                /*Iterate on "from" nodes*/
	                for(l=0;l<nbNodes;l++){
	                    /*If l is not absorbing, compute E[N(l)|wlen]**/
	                    if(!rows[l].absorbing){
	                        /*Iterate on the "to" nodes*/
	                        selem = rows[l].first;
	                        while(selem){
	                            if(!rows[selem->c].absorbing){
                                
	                                /*****************************/
	                                /*Transient destination nodes*/
	                                /*****************************/
                                
	                                for(t=0,tmp=0;t<wlen-1;t++){
	                                    tmp += latF[l][t]*(selem->val)*latB[selem->c][t+1];
	                                }
                                
	                                /*Update the expectation matrix*/
	                                update      = mass_i*tmp/p_wlen;
	                                selem->cur += update;
	                                selem->ept += update;
	                                N[l]       += update;
	                                expWlen    += update;
	                            }
	                            else{
                                
	                                /*****************************/
	                                /*Absorbing destination nodes*/
	                                /*****************************/
                                
	                                tmp = latF[l][wlen-1]*(selem->val);
                                
	                                /*Update the expectation matrix*/
	                                update      = mass_i*tmp/p_wlen;
	                                selem->cur += update;
	                                selem->ept += update,
	                                N[l]       += update;
	                                expWlen    += update;
	                            }
	                            selem = selem->nextC;
	                        }/*End while(selem)*/
	                    }/*End if(non-absorbing from node)*/
	                }/*End for(from node)*/
	            }/*End if(P[wlen] > 0)*/
	        }/*End if(upto)*/
	        /*If the graph is undirected compute |E_ij - E_ji|*/
	        if(undirected) diffEPT_sparse(rows,cols,nbNodes);
	        if(verbose > 0)
	            printf("-------------------------------------------------------\n");        
	    }/*End for(relevant nodes)*/
	}/*End for(inflation)*/	

    /*Stop the CPU chronometer*/
    stop_cpu = clock();
    
    /*Print informations*/
    if (verbose > 0){
        printf("CPU Time (sec)                  : %e\n",((double)(stop_cpu - start_cpu)/CLOCKS_PER_SEC));
        printf("Absorption probability ratio    : %e\n",mass_abs);
        printf("Expected walk length            : %e\n",expWlen);       
    }
    
    /*Write the EPT to output files*/
    if(oflg){
        writeVec_sparse(N,nbNodes,outNPT);
        writeMEPT_sparse(outEPT,rows,nbNodes);
    }
    
    /*Normalize the EPT by expWlen*/
    vecNormalize(N,nbNodes,expWlen);
    matNormalize_sparse(rows,nbNodes,expWlen);
    
    /*Write to output files*/
    if(oflg){
        writeVec_sparse(N,nbNodes,outNPT_Norm);
        writeMEPT_sparse(outEPT_Norm,rows,nbNodes);
    }            
    
    /*If the graph is undirected write |E_ij - E_ji|*/
    if(undirected && oflg){
        writeDiff_sparse(outDIF,rows,nbNodes);
    }    
    
    /*If required write a graphviz dot file*/
    if(genDot){
        writeDot(outDOT,rows,nbNodes,undirected);
    }    
    
    /*Display the relative size of the extracted subgraph*/
    if(verbose > 0){
        matrixNNZ_sparse(rows,nbNodes,&nnz,&nnr);
        printf("Subgraph size (#edges)          : %i\n",nnz);
        printf("Subgraph size (#nodes)          : %i\n",nnr);
        printf("Subgraph relative size (#edges) : %e\n",(float)nnz/(float)nbEdges);        
        printf("-------------------------------------------------------\n");
        printf("Prototype by J.Callut, INGI dept. 2008 (c)\n");
	printf("Last update by P. Dupont, UCL/ICTEAM/INGI 2010\n\n");
    }    
    
    /*Exit successfully :-)*/
    return EXIT_SUCCESS;
}

/*Main function for the stand-alone command line tool*/
int main (int argc, char* argv[]){
  
    /*Call the lkwalk function with command line arguments*/
    int exit_value = lkwalk (argc,argv);
    
    /*Exit with the value return by lkwalk()*/
    exit(exit_value);
}

