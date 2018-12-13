#!/usr/bin/python
#-*-coding: utf-8-*-

"""Name : supported_organisms_ensembl
Version: v0.1
Description: this program was developed in order to get informations 
on all the studied species available on the MICROME server. It runs only in Python 2.7.
Authors: alexandra.bomane@laposte.net, mel.colombier@free.fr
Usage: python supported_organisms_ensembl.py [-o outfile_name] [--output_file outfile_name]
Options: 
    -h  Display full help message
    -help
        Same as -h
    -o outfile_name
        If no output file is specified, the standard output is used.
    --output_file outfile_name
        Same as -o
"""

import httplib2 #To send the request using http protocol
import sys
import time #To add the current date to the standard file name
import json #For parsing the server response
import argparse #To add options

if __name__=='__main__':
 
    print "; Job started", time.strftime("%Y-%m-%d_%Hh%M")
    
    #.cache is a file containing in memory all requests parameters such as the status, the url, the returned JSON object...
    http = httplib2.Http(".cache")
    
    #Specify the url of the server 
    resp, content = http.request("http://www.microme.eu/reaction_matrix/api/genomes.json", method="GET", headers={"Content-Type":"application/json"})
    
    #Check the status of the request: 
    if not resp.status == 200:
        print "Invalid response: ", resp.status
        sys.exit()
    
    #Response return as a JSON object : here an array containing dictionaries with each organism's informations 
    decoded = json.loads(content)         
    
    #Current date
    date = time.strftime("%Y_%m_%d")
    
    #Creation of arguments
    parser = argparse.ArgumentParser(add_help = True)
    parser.add_argument('-o','--output_file', action = 'store', dest = 'outfile_n', default = False, help = 'The user can choose the outfile name')
    args = parser.parse_args()
    
    outfile_n, outfileN = args.outfile_n, args.outfileN
    
    # Choose of the outfile name
    if outfile_n:
        outfile = outfile_n  #File name chosen by the user
    else:
        outfile = "supported_organisms_microme_" + date + ".tab" #Standard file name according to the date
    
    #Creation of the tabulated file
    f = open(outfile, "w")
    
    f.write("Assembly Name" + '\t' + "Name" + '\t' + "Supperregnum" + '\t' + "Taxonomy" + '\t' + "Uniprot Code" + '\t' + "Assembly Ac" + '\n')
    
    for i in range(len(decoded)):
        dict_species = decoded[i]     
    
    #file_content variable contents all strings matching to each value of the dictionary for a given species separated by tabulations    
        file_content = str(dict_species['assembly_name']) + '\t' + str(dict_species['name']) + '\t' + str(dict_species['superregnum']) + '\t' + str(dict_species['taxonomy_id']) + '\t' + str(dict_species['uniprot_oscode']) + '\t' + str(dict_species['assembly_ac'])
        
        f.write(file_content + '\n')
    
    f.close()
    
    print "; Job done, Job name :", outfile, time.strftime("%Y-%m-%d_%Hh%M")


