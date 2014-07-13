#!/usr/bin/python
#-*-coding: utf-8-*-

### imports ###
#Used to use print 3.x function
from __future__ import print_function
import os
import sys
import json
import platform
import time #To get the current date


from workflows_functions import *


#if using python 3.x
try:
    from urllib.request import urlopen
#else
except ImportError:
    from urllib import urlopen


python_version = platform.python_version_tuple()

### functions needed for retrieve Species###


def connection_server_species(v):

    """
    This function will send a request to a server and then return the answer
    :param v: verbosity index
    :return: A Json object. The response from the server to the request
    """

    url = "http://mgalileo.genoscope.cns.fr/api/microscope/organisms/list"

    if v >= 3:
        print("Sending a request to the http server...")
    try:
        response = urlopen(url)
        data = response.read()
        decoded_data = data.decode('utf-8')
        json_dict = json.loads(decoded_data)

        if v >= 3:
            print("Server sent a response")
    except:
        print("Can't reach server at : %s" % url)
        print(sys.exc_info()[1])
        if os.path.exists('org'):
            os.remove('org')
        sys.exit()

    return json_dict


    
    
    
def retrieve_microscope_species(outfile):
  
     
    print( "; Job started", time.strftime("%Y-%m-%d %H:%M:%S") + '\n')
     
    #Send the request to the Genoscope server and get a Json object
    org_supported_json = RequestGenoscope("/organisms/list")

    #Get a Python object from a JSON object : here an array containing dictionaries with each organism's informations 
    org_supported = json.loads(org_supported_json)         
    
    
    ########According to the user choice, the program runs regarding to a few columns or all columns 
    out_columns = ['ospeciesCode','odirname','ogram','oid','osynonym','otaxon','oreplicNb','okingdom','oname','ostrain']
    
    #Variable containing all strings matching to each value of the dictionary for a given species separated by tabulations
    file_content = ""
    
    for i in range(len(org_supported)):
    #Dictionary containing informations about a given species
            current_species_dict = org_supported[i] 
                   
    #Feeding of file_content variable 
            for j in range(len(out_columns)):
                if out_columns[j] in current_species_dict.keys():
                    file_content += unicode(current_species_dict[out_columns[j]]) + '\t'
            file_content += '\n'
    
    if outfile:
      #Creation of the tabulated file
      Save_Tabulated_File(out_columns,outfile,file_content)
    else:
      Display_Tabulated_Results(out_columns,file_content)
    
    print( "; Job done", time.strftime("%Y-%m-%d %H:%M:%S") ) 

  
  
 ##################################################################
 ## REQUEST 2 :
 
def genomicObjects_by_orgname(organism, outfile):
 
    ## Job_started variable is created here to get the real moment of the beginning of the program. It will be used in the output file
    Job_started = time.strftime("%Y-%m-%d %H:%M:%S")

    print( "; Job started", Job_started)

    # Send the request to the Genoscope server and get a Json object
    print( 'Sending the request to the server and checking the request status')
    supported_org_json = RequestGenoscope("/organisms/list")

    # From the JSON object we get a Python object : a list containing dictionaries containing dictionaries
    print( 'Treatment of the server response')
    supported_org = json.loads(supported_org_json)

    # This variable is initialized in order to get the organism Genoscope ID (oid)
    org_id = ''

    # Identify the query organism in the list of supported organisms and collect some of its Genoscope ID (oid)
    for i in range (len(supported_org)):
        current_species_dict = supported_org[i]
        if organism in current_species_dict[u'odirname']:
            org_id = str(current_species_dict[u'oid'])

    # Check that the query organism is found in the list of supported
    # organisms : if it isn't listed, an error message is displayed on
    # the STDOUT and we quit the program
    if org_id == '' :
        print( "Organism" + ' ' + organism + ' ' + "is not supported at Genoscope server. \nTo get the list of supported organisms use : 'python microscope_get.py organisms'")
        sys.exit()

    ## New request to get informations about organism genes and get a new Json object
    org_genes_json = RequestGenoscope("/genomes/genomic-objects/list-from-organism/" + org_id)

    ## Get a Python object from a JSON object : here it is an array/lis containing dictionaries
    org_genes = json.loads(org_genes_json)

    ## Variable containing all strings matching to interesting values from dict_gene_val for a given gene separated by tabulations
    file_content = ''
    print( 'Obtention of the response data')

    ## Error 1 counters (list) ---> Error = there is a reaction but there is not an EC number
    print( 'Errors search')
    errors_Reaction_No_Ec = []
    reactionInerror = []
    ## Error 2 counter (list) ---> Error = there is an ID but no gene name
    errors_ID_No_GeneName = []


    ## Column headers list
    #list_columnHeaders = ['refSeq','uniprot','trembl','goOriId','evidence','protId']
    list_columnHeaders = ['refSeq','uniprot','trembl','goOriIdE','goTypeE','goFrameE','goBeginE','goEndE','goStatusE','goMutationE','goLabelE','goGeneNameE','goSynonymsE','goProductE','goProductTypeE','goFunctionE','goLocalizationE','goClassE','goLengthE','goEcE','goEvidenceE','goCreationDateE','goAmigeneStatusE','goPmidE','goProcessE','protId','goOriId']
    ## From the Python object (array) we get the dictionary containing informations about a given gene
    print( 'Data analysis...')

    ## This dictionary contains column headers as keys and (in this
    ## worder) the list of data matching with the column header, the
    ## counter of Non Null data, the counter of Null data and the
    ## counter of Unique data as values
    dict_counters = {'refSeq': [[],0,0,0],'uniprot':[[],0,0,0] ,'trembl':[[],0,0,0] ,'goOriId':[[],0,0,0],'protId':[[],0,0,0]}
    for i in range(len(org_genes)):

        ## dict_gene is a dictionary containing informations about a
        ## given gene. dict_gene key is u'genomicObject' and its value
        ## is a dictionary containing the given gene characteristics
        dict_gene = org_genes[i]

        ## Feeding of file_content variable
#        file_content += unicode(dict_gene_val[u'goGeneName']) + '\t' + unicode(dict_gene_val[u'goSynonyms']) + '\t' + unicode(dict_gene_val[u'goEc']) + '\t' + unicode(dict_gene_val[u'goOriId']) + '\t' + unicode(dict_gene_val[u'goReaction']) + '\t' + unicode(dict_gene_val[u'goId']) + '\t' + unicode(dict_gene_val[u'goFunction']) + '\t' + unicode(dict_gene_val[u'goProductType']) + '\t' + unicode(dict_gene_val[u'goProcess']) + '\t' + unicode(dict_gene_val[u'goProduct']) + '\n'
        file_content += '\t'.join([unicode(dict_gene[u'refSeq']),
                                   unicode(dict_gene[u'uniprot']),
                                   unicode(dict_gene[u'trembl']),
                                   unicode(dict_gene[u'goOriIdE']),
                                   unicode(dict_gene[u'goTypeE']),
                                   unicode(dict_gene[u'goFrameE']),
                                   unicode(dict_gene[u'goBeginE']),
                                   unicode(dict_gene[u'goEndE']),
                                   unicode(dict_gene[u'goStatusE']),
                                   unicode(dict_gene[u'goMutationE']),
                                   unicode(dict_gene[u'goLabelE']),
                                   unicode(dict_gene[u'goGeneNameE']),
                                   unicode(dict_gene[u'goSynonymsE']),
                                   unicode(dict_gene[u'goProductE']),
                                   unicode(dict_gene[u'goProductTypeE']),
                                   unicode(dict_gene[u'goFunctionE']),
                                   unicode(dict_gene[u'goLocalizationE']),
                                   unicode(dict_gene[u'goClassE']),
                                   unicode(dict_gene[u'goLengthE']),
                                   unicode(dict_gene[u'goEcE']),
                                   unicode(dict_gene[u'goEvidenceE']),
                                   unicode(dict_gene[u'goCreationDateE']),
                                   unicode(dict_gene[u'goAmigeneStatusE']),
                                   unicode(dict_gene[u'goPmidE']),
                                   unicode(dict_gene[u'goProcessE']),
                                   unicode(dict_gene[u'protId']),
                                   unicode(dict_gene[u'goOriId'])])
        file_content += '\n'

    ## Column headers dictionary
    dict_columnHeaders = dict([(1,"refSeq\trefSeq ID"),
                               (2,"uniprot\tuniprot ID"),
                               (3,"trembl\ttrembl"),
                               (4,"goOriIdE\tgene Genoscope ID"),
                               (5,"goTypeE"),
                               (6,"goFrameE"),
                               (7,"goBeginE"),
                               (8,"goEndE"),
                               (9,"goStatusE"),
                               (10,"goMutationE"),
                               (11,"goLabelE"),
                               (12,"goGeneNameE"),
                               (13,"goSynonymsE"),
                               (14,"goProductE"),
                               (15,"goProductTypeE"),
                               (16,"goFunctionE"),
                               (17,"goLocalizationE"),
                               (18,"goClassE"),
                               (19,"goLengthE"),
                               (20,"goEcE"),
                               (21,"goEvidenceE"),
                               (22,"goCreationDateE"),
                               (23,"goAmigeneStatusE"),
                               (24,"goPmidE"),
                               (25,"goProcessE"),
                               (26,"protId"),
                               (27,"goOriId")])

    
    
    if outfile:
      #Creation of the tabulated file
      Save_Tabulated_File(list_columnHeaders,outfile,file_content)
    else:
      Display_Tabulated_Results(list_columnHeaders,file_content)
     

    
    print( "\n; Job done", time.strftime("%Y-%m-%d %H:%M:%S"))

  
  
def get_MetaCycReaction(reaction):
    ## Job_started variable is created here to get the real moment of the beginning of the program. It will be used in the output file
    Job_started = time.strftime("%Y-%m-%d %H:%M:%S")

    ## The user has to enter an organism name to launch the program else we quit it
    if reaction == None:
        print( "Error: reaction name is mandatory. \nUse the option:\n\t-mrID my_Reaction_Name (ex: '6-ACETYLGLUCOSE-DEACETYLASE-RXN')")
        sys.exit()


    ## New request to get informations about organism genes and get a new Json object
    mrId_genes_json = RequestGenoscope("/networks/microcyc/reactions/get/" + reaction)

    ## Get a Python object from a JSON object : here it is an array/lis containing dictionaries
    org_genes = json.loads(mrId_genes_json)

    ## Variable containing all strings matching to interesting values from dict_gene_val for a given gene separated by tabulations
    file_content = ''
    print( 'Obtention of the response data')

   

    ## Column headers list
    #list_columnHeaders = ['refSeq','uniprot','trembl','goOriId','evidence','protId']
    list_columnHeaders = ['mrIdE','mrNameE','mrEcE','mrOfficialEcE','mrSpontaneousE','coefficientE','natureE','compartmentE','cidE']
    ## From the Python object (array) we get the dictionary containing informations about a given gene
    print( 'Data analysis...')

    ## This dictionary contains column headers as keys and (in this
    ## worder) the list of data matching with the column header, the
    ## counter of Non Null data, the counter of Null data and the
    ## counter of Unique data as values
#    for i in range(len(org_genes)):

        ## dict_gene is a dictionary containing informations about a
        ## given gene. dict_gene key is u'genomicObject' and its value
        ## is a dictionary containing the given gene characteristics
    dict_gene = org_genes
    dict_gene_val = dict_gene[u'mcrEs']

    for j in range(len(dict_gene_val)):
        dict_gene_valj = dict_gene_val[j]
            ## Feeding of file_content variable
        file_content += '\t'.join([unicode(dict_gene[u'mrIdE']),
                                       unicode(dict_gene[u'mrNameE']),
                                       unicode(dict_gene[u'mrEcE']),
                                       unicode(dict_gene[u'mrOfficialEcE']),
                                       unicode(dict_gene[u'mrSpontaneousE']),
                                       unicode(dict_gene_valj[u'coefficientE']),
                                       unicode(dict_gene_valj[u'natureE']),
                                       unicode(dict_gene_valj[u'compartmentE']),
                                       unicode(dict_gene_valj[u'cidE'])])
        file_content += '\n'

    ## Column headers dictionary
    dict_columnHeaders = dict([(1,"mrIdE"),
                               (2,"mrNameE"),
                               (3,"mrEcE"),
                               (4,"mrOfficialEcE"),
                               (5,"mrSpontaneaousE"),
                               (6,"coefficientE"),
                               (7,"natureE"),
                               (8,"compartmentE"),
                               (9,"cidE")])

    if outfile:
      #Creation of the tabulated file
      Save_Tabulated_File(out_columns,outfile,file_content)
    else:
      Display_Tabulated_Results(list_columnHeaders,file_content)
 
    print( "\n; Job done", time.strftime("%Y-%m-%d %H:%M:%S"))

 
 
  
def get_MetaCycReactionFromOrg(organism):
      ## Job_started variable is created here to get the real moment of the beginning of the program. It will be used in the output file
    Job_started = time.strftime("%Y-%m-%d %H:%M:%S")

    ## The user has to enter an organism name to launch the program else we quit it
    if organism == None:
        print( "Error: organism name is mandatory. \nUse the option:\n\t-org My_Organism_Name")
        sys.exit()

    ## The job can start here
    else :
        print( "; Job started", Job_started)

        # Send the request to the Genoscope server and get a Json object
        print( 'Sending the request to the server and checking the request status')
        supported_org_json = RequestGenoscope("/organisms/list")

        # From the JSON object we get a Python object : a list containing dictionaries containing dictionaries
        print( 'Treatment of the server response')
        supported_org = json.loads(supported_org_json)

    # This variable is initialized in order to get the organism Genoscope ID (oid)
    org_id = ''

    # Identify the query organism in the list of supported organisms and collect some of its Genoscope ID (oid)
    for i in range (len(supported_org)):
        current_species_dict = supported_org[i]
        if organism in current_species_dict[u'odirname']:
            org_id = str(current_species_dict[u'oid'])

    # Check that the query organism is found in the list of supported
    # organisms : if it isn't listed, an error message is displayed on
    # the STDOUT and we quit the program
    if org_id == '' :
        print( "Organism" + ' ' + organism + ' ' + "is not supported at Genoscope server. \nTo get the list of supported organisms, use the command : microscope_get organisms")
        sys.exit()

    ## New request to get informations about organism genes and get a new Json object
    org_genes_json = RequestGenoscope("/networks/microcyc/reactions/list/" + org_id)

    ## Get a Python object from a JSON object : here it is an array/lis containing dictionaries
    org_genes = json.loads(org_genes_json)

    ## Variable containing all strings matching to interesting values from dict_gene_val for a given gene separated by tabulations
    file_content = ''
    print( 'Obtention of the response data')


    ## Column headers list
    list_columnHeaders = ['reactionId','genes']

    ## From the Python object (array) we get the dictionary containing informations about a given gene
    print( 'Data analysis...')

    
    for i in range(len(org_genes)):

        ## dict_gene is a dictionary containing informations about a
        ## given gene. dict_gene key is u'genomicObject' and its value
        ## is a dictionary containing the given gene characteristics
        dict_gene = org_genes[i]
        dict_gene_val = dict_gene[u'genes']
        for j in range(len(dict_gene_val)):
        #    dict_gene_val = dict_gene[u'genes']
        ## Feeding of file_content variable
        
            file_content += '\t'.join([unicode(dict_gene[u'reactionId']),
                                       unicode(dict_gene_val[j])])
            file_content += '\n'

    ## Column headers dictionary
    dict_columnHeaders = dict([(1,"reaction\treaction ID"),
                               (2,"gene\tgene ID")])

    if outfile:
      #Creation of the tabulated file
      Save_Tabulated_File(out_columns,outfile,file_content)
    else:
      Display_Tabulated_Results(list_columnHeaders,file_content)
 
    print( "\n; Job done", time.strftime("%Y-%m-%d %H:%M:%S"))
    
    
    
def allReactionList():
      ## Job_started variable is created here to get the real moment of the beginning of the program. It will be used in the output file
    Job_started = time.strftime("%Y-%m-%d %H:%M:%S")


    ## New request to get informations about organism genes and get a new Json object
    mrId_genes_json = RequestGenoscope("/networks/microcyc/reactions/list/")

    ## Get a Python object from a JSON object : here it is an array/lis containing dictionaries
    org_genes = json.loads(mrId_genes_json)

    ## Variable containing all strings matching to interesting values from dict_gene_val for a given gene separated by tabulations
    file_content = ''
    print( 'Obtention of the response data')

   

    ## Column headers list
    #list_columnHeaders = ['refSeq','uniprot','trembl','goOriId','evidence','protId']
    list_columnHeaders = ['mrIdE','mrNameE','mrEcE','mrOfficialEcE','mrSpontaneousE','coefficientE','natureE','compartmentE','cidE']
    ## From the Python object (array) we get the dictionary containing informations about a given gene
    print( 'Data analysis...')

    ## This dictionary contains column headers as keys and (in this
    ## worder) the list of data matching with the column header, the
    ## counter of Non Null data, the counter of Null data and the
    ## counter of Unique data as values
#    for i in range(len(org_genes)):

        ## dict_gene is a dictionary containing informations about a
        ## given gene. dict_gene key is u'genomicObject' and its value
        ## is a dictionary containing the given gene characteristics
    
    for i in range(len(org_genes)):
      dict_gene = org_genes[i]
      dict_gene_val = dict_gene[u'mcrEs']

      for j in range(len(dict_gene_val)):
	dict_gene_valj = dict_gene_val[j]
	# Feeding of file_content variable
	file_content += '\t'.join([unicode(dict_gene[u'mrIdE']),
                                  unicode(dict_gene[u'mrNameE']),
                                  unicode(dict_gene[u'mrEcE']),
                                  unicode(dict_gene[u'mrOfficialEcE']),
                                  unicode(dict_gene[u'mrSpontaneousE']),
                                  unicode(dict_gene_valj[u'coefficientE']),
                                  unicode(dict_gene_valj[u'natureE']),
                                  unicode(dict_gene_valj[u'compartmentE']),
                                  unicode(dict_gene_valj[u'cidE'])])
        file_content += '\n'

    ## Column headers dictionary
    dict_columnHeaders = dict([(1,"mrIdE"),
                               (2,"mrNameE"),
                               (3,"mrEcE"),
                               (4,"mrOfficialEcE"),
                               (5,"mrSpontaneaousE"),
                               (6,"coefficientE"),
                               (7,"natureE"),
                               (8,"compartmentE"),
                               (9,"cidE")])

    if outfile:
      #Creation of the tabulated file
      Save_Tabulated_File(out_columns,outfile,file_content)
    else:
      Display_Tabulated_Results(list_columnHeaders,file_content)
 
    print( "\n; Job done", time.strftime("%Y-%m-%d %H:%M:%S"))
    
