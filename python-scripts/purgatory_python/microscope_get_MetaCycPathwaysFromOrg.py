#!/usr/bin/python
#-*-coding: utf-8-*-

"""Version : v0.2

Description :
    This program was developed in order to get informations on all studied species'genes available on the Genoscope server.
    It runs only in Python 2.7 and with the workflows_functions.py file.

Authors :
    - Alexandra BOMANE    alexandra.bomane@laposte.net
    - Melanie COLOMBIER    mel.colombier@free.fr
    - Aurélie BERGON       aurelie.bergon@inserm.fr
Usage :
     microscope_get_microscope_get_MetaCycReactions.py -org organism_name  [-o outfile_name] [--output_file outfile_name]

Options :
    -h | -help
        Display full help message.
    -org organism_name
        Selected organism. Mandatory argument.
        To get available organisms, refer to supported_organisms_genoscope.py
    -o | --output_file  outfile_name
        If no output file is specified, the result is printed on the standard output.
        Else, 2 files are returned : a tabulated file presenting the server response and an error report.

Execution :
    The program launch a first request to the Genoscope server.

    This request is the same as in supported_organisms_genoscope.py (http://mgalileo.genoscope.cns.fr/api/microscope/organisms/list).
    It returns all public organisms of MicroScope database in json format.

    The Json object is parsed in order to get the Genoscope Id of the chosen organism.

    This Id allows the program to send a second request to the server.
    This request was established thanks the Genoscope API "/genomes/genomic-objects/list/Organism'sGenoscopeId" that returns a Json object too.
    The second request response is parsed to make the result more readable.

    Eventually, according to the use of -o | --output_file option, we have 2 possibilities :
        - The result is displayed on the standard output and the error detail is displayed in the standard error.
        - The result is written in a tabulated file and a precise error report file is returned too.
"""

import json #To parse the server response
import time #To add the current date to the standard file name
import argparse #To add options
import sys  #To print in STDERR and to get the command line
import socket   #To get host name
import os   #To get the path of the output file directory
from collections import Counter #To do the statistic analysis
from workflows_functions import RequestGenoscope    #To send our request to the Genoscope server
from workflows_functions import Display_Tabulated_Results    #To display the result on the standard output


if __name__=='__main__':

    ## Job_started variable is created here to get the real moment of the beginning of the program. It will be used in the output file
    Job_started = time.strftime("%Y-%m-%d %H:%M:%S")

    ## Declaration of arguments to argparse
    parser = argparse.ArgumentParser(add_help = True)
    parser.add_argument('-org', action = 'store' , dest = 'organismID', default = None, help = 'Select a reaction (mandatory; ex: "Escherichia_coli_K-12_DH10B")')
    parser.add_argument('-o','--output_file', action = 'store', dest = 'outfilen', default = False, help = 'The user can choose the outfile name (optional)')
    args = parser.parse_args()

    ## Variable containing organism chosen by the user : it's the mandatory argument
    organism = args.organismID

    ## Variable containing the file name chosen by the user : if not specified, the output will be displayed on the STDOUT
    outfilen = args.outfilen

    ## To get the line command
    command_args = ""
    for arg in sys.argv:
        command_args = command_args + arg + ' '

    ## The user has to enter an organism name to launch the program else we quit it
    if organism == None:
        print "Error: organism name is mandatory. \nUse the option:\n\t-org My_Organism_Name"
        sys.exit()

    ## The job can start here
    else :
        print "; Job started", Job_started

        # Send the request to the Genoscope server and get a Json object
        print '; Sending the request to the server and checking the request status'
        supported_org_json = RequestGenoscope("/organisms/list")

        # From the JSON object we get a Python object : a list containing dictionaries containing dictionaries
        print '; Treatment of the server response'
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
        print "Organism" + ' ' + organism + ' ' + "is not supported at Genoscope server. \nTo get the list of supported organisms, use the command : supported_organisms_genoscope.py"
        sys.exit()
    
    print("organism:",org_id,organism)

    ## New request to get informations about organism genes and get a new Json object
    mpId_genes_json = RequestGenoscope("/networks/microcyc/pathways/list/" + org_id)
    
    ## Get a Python object from a JSON object : here it is an array/lis containing dictionaries
    org_genes = json.loads(mpId_genes_json)

    ## Variable containing all strings matching to interesting values from dict_gene_val for a given gene separated by tabulations
    file_content = ''
    print 'Obtention of the response data'

   

    ## Column headers list
    #list_columnHeaders = ['refSeq','uniprot','trembl','goOriId','evidence','protId']
    list_columnHeaders = ['mpIdE','mpNameE','mpTypeE','mpHasReactionE','mrIdE','mrNameE','mrEcE','mrOfficialEcE','mrSpontaneousE','coefficientE','natureE','compartmentE','cidE']
    ## From the Python object (array) we get the dictionary containing informations about a given gene
    print 'Data analysis...'

    ## This dictionary contains column headers as keys and (in this
    ## worder) the list of data matching with the column header, the
    ## counter of Non Null data, the counter of Null data and the
    ## counter of Unique data as values
#    for i in range(len(org_genes)):

        ## dict_gene is a dictionary containing informations about a
        ## given gene. dict_gene key is u'genomicObject' and its value
        ## is a dictionary containing the given gene characteristics
    

    for p in range(len(org_genes)):
        dict_pways = org_genes[p]
        for r in range(len(dict_pways)):
            dict_react=dict_pways[u'metaCycReactionExtensions']
            for e in range(len(dict_react)):
                dict_react_val=dict_react[e]
                dict_equat = dict_react_val[u'mcrEs']
                for j in range(len(dict_equat)):
                    dict_equat_val = dict_equat[j]
                    ## Feeding of file_content variable
                    file_content += '\t'.join([unicode(dict_pways[u'mpIdE']),
                                           unicode(dict_pways[u'mpNameE']),
                                           unicode(dict_pways[u'mpTypeE']),
                                           unicode(dict_pways[u'mpHasReactionE']),
                                           unicode(dict_react_val[u'mrIdE']),
                                           unicode(dict_react_val[u'mrNameE']),
                                           unicode(dict_react_val[u'mrEcE']),
                                           unicode(dict_react_val[u'mrOfficialEcE']),
                                           unicode(dict_react_val[u'mrSpontaneousE']),
                                           unicode(dict_equat_val[u'coefficientE']),
                                           unicode(dict_equat_val[u'natureE']),
                                           unicode(dict_equat_val[u'compartmentE']),
                                           unicode(dict_equat_val[u'cidE'])])
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


    # The user uses the -o | --output_file option
    if outfilen:
        # File name chosen by the user
        print 'Creation of the tabulated files'
        outfile = outfilen

        # Job_done variable is created here to get the real moment of the end of the program. It will be used in the output file
        Job_done = time.strftime("%Y-%m-%d %H:%M:%S")


        # Tabulated file written
        f = open(outfile,'w')
        f.write(file_content.encode('utf-8'))
        f.close()

        print "; Job done", time.strftime("%Y-%m-%d %H:%M:%S")

    # The user doesn't use the -o | --output_file option : result is displayed on the standard output
    else:
        Display_Tabulated_Results(list_columnHeaders,file_content)

        print "\n; Job done", time.strftime("%Y-%m-%d %H:%M:%S")
