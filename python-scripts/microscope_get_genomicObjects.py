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
     microscope_get_genomicObjects.py -org organism_name  [-o outfile_name] [--output_file outfile_name]

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
    parser.add_argument('-org', action = 'store' , dest = 'organism', default = None, help = 'Selected organism (mandatory)')
    parser.add_argument('-o','--output_file', action = 'store', dest = 'outfilen', default = False, help = 'The user can choose the outfile name (optional)')
    args = parser.parse_args()

    ## Variable containing organism chosen by the user : it's the mandatory argument
    organism = args.organism

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
        print 'Sending the request to the server and checking the request status'
        supported_org_json = RequestGenoscope("/organisms/list")

        # From the JSON object we get a Python object : a list containing dictionaries containing dictionaries
        print 'Treatment of the server response'
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

    ## New request to get informations about organism genes and get a new Json object
    org_genes_json = RequestGenoscope("/genomes/genomic-objects/list-from-organism/" + org_id)

    ## Get a Python object from a JSON object : here it is an array/lis containing dictionaries
    org_genes = json.loads(org_genes_json)

    ## Variable containing all strings matching to interesting values from dict_gene_val for a given gene separated by tabulations
    file_content = ''
    print 'Obtention of the response data'

    ## Error 1 counters (list) ---> Error = there is a reaction but there is not an EC number
    print 'Errors search'
    errors_Reaction_No_Ec = []
    reactionInerror = []
    ## Error 2 counter (list) ---> Error = there is an ID but no gene name
    errors_ID_No_GeneName = []


    ## Column headers list
    #list_columnHeaders = ['refSeq','uniprot','trembl','goOriId','evidence','protId']
    list_columnHeaders = ['refSeq','uniprot','trembl','goOriIdE','goTypeE','goFrameE','goBeginE','goEndE','goStatusE','goMutationE','goLabelE','goGeneNameE','goSynonymsE','goProductE','goProductTypeE','goFunctionE','goLocalizationE','goClassE','goLengthE','goEcE','goEvidenceE','goCreationDateE','goAmigeneStatusE','goPmidE','goProcessE','protId','goOriId']
    ## From the Python object (array) we get the dictionary containing informations about a given gene
    print 'Data analysis...'

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

        ## The user uses the -o | --output_file option : calculation
        ## of statistics to write them in an error report
        if outfilen:
            # dict_counters updated according to dict_gene_val data
            for col_name in list_columnHeaders:
                dict_counters[col_name][0].append(dict_gene["u'"+col_name+"'"])
                if dict_gene["u'"+col_name+"'"] != None:
                    dict_counters[col_name][1] += 1 #Non Null counter
                else:
                    dict_counters[col_name][2] += 1 #Null counter
            for col_name in list_columnHeaders:
                dict_counters[col_name][3] = dict(Counter(dict_counters[col_name][0])).values().count(1) #Unique counter

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


    # The user uses the -o | --output_file option
    if outfilen:
        # File name chosen by the user
        print 'Creation of the tabulated files'
        outfile = outfilen

        # Job_done variable is created here to get the real moment of the end of the program. It will be used in the output file
        Job_done = time.strftime("%Y-%m-%d %H:%M:%S")


        # Tabulated file written
        f = open(outfile,'w')
        f.write(command_args + '\n' + ';\tOutput file\t' + os.path.abspath(outfile) + '\n')
        f.write(';\tErrors report\t' + os.path.abspath(error_file_name) + '\n')
        f.write(";\tOutput format\ttab\n")
        f.write(";\tErrors report format\ttxt\n")
        f.write(";\tOrganism\t" + organism + '\n;\n')

        ## Write the column content (number of values for each field)
        f.write(";\tColumn contents\n")
        f.write('\t'.join([';','col_nb','col_head','nonull','uniq','null']) + '\n')
        for header in list_columnHeaders:
            f.write(';\t' 
                    + '\t'.join([str(list_columnHeaders.index(header) + 1), 
                                 "%-14s" % header, str(dict_counters[header][1]), 
                                 str(dict_counters[header][3]), 
                                 str(dict_counters[header][2])]) 
                    + '\n')

        ## Write column headers
        f.write(";\t\n;\tColumn headers\n")
        for key in sorted(dict_columnHeaders.iterkeys()):
            f.write(";\t%d\t%s" % (key, dict_columnHeaders[key]) + '\n')
        f.write('#' + list_columnHeaders[0] + '\t' + '\t'.join(list_columnHeaders[1:]) + '\n')
        f.write(file_content.encode('utf-8'))
        f.write(';\tHost name\t' + socket.gethostname() + '\n' + ';\tJob started\t' + Job_started + '\n;\tJob done\t' + Job_done)
        f.close()

        print "; Output names :", os.path.abspath(outfile) + '\n' + os.path.abspath(error_file_name)
        print "; Job done", time.strftime("%Y-%m-%d %H:%M:%S")

    # The user doesn't use the -o | --output_file option : result is displayed on the standard output
    else:
        Display_Tabulated_Results(list_columnHeaders,file_content)

        print "\n; Job done", time.strftime("%Y-%m-%d %H:%M:%S")
