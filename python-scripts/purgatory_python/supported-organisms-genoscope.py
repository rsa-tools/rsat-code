#!/usr/bin/python
# -*- coding: utf-8 -*-

######################################################
### -------------------- HELP -------------------- ###

'''
NAME
  %(progname)s

VERSION
  %(version)s

CATEGORY
  Metabolism

AUTHORS
  Quentin DA COSTA <quentin.da-costa@etumel.univmed.fr>
  Jonathan VERNEAU <jonathan.verneau@etumel.univmed.fr>

DESCRIPTION
  Get the list of supported organisms from the MicroScope database
(http://www.genoscope.cns.fr/agc/microscope/).

OUTPUT FORMAT
  The program returns a tab-delimited text file with one row per
organism.

Columns:

1. ms_Organism_id
   Identifier of the organism in the MicroScope database (inherited
   from MetaCyc).

2. ms_Number
   Number of the organism in the list (the order reflects the import
   history).

3. Organism_name
   Scientific name of the organism.

4. Taxon_id
   Taxonomic ID (TAXID) as defined in the Taxonomy database at NCBI
   (http://www.ncbi.nlm.nih.gov/taxonomy/).
'''

##########################################################
### -------------------- PACKAGES -------------------- ###

import suds
import time
import argparse
import sys
import os

###########################################################
### -------------------- FUNCTIONS -------------------- ###

## Open a connection to the GenoScope Web services
def Connection():
    url = 'http://www.genoscope.cns.fr/microcyc-webservice/services/fr.genoscope.agc.microcyc.business.MicrocycService.wsdl'
    return suds.client.Client(url)
    


def Create_file(fileName):
    client=Connection()
    ## Open and write a new file
    f=open(fileName,"w")

    ## Print generic header (program name, execution date)
    f.write('; ' + str.join(" ", sys.argv)+'\n')
    f.write(';\n')
    f.write('; Time: '+time.strftime('%Y-%m-%d %H:%M\n',time.localtime()))
    f.write(';\n')

    ## Print column content
    header="""; Column content
;\t1\tms_Organism_id\tIdentifier of the organism in the MicroScope database (inherited from MetaCyc)
;\t2\tms_Number\tNumber of entry in the database (reflects history of import in MicroScope)
;\t3\tOrganism_name\tGenus + species name
;\t4\tTaxon_id\tIAXID as specified in the NCBI taxonomy database
;
"""
    f.write(header)

    ## Print header line
    i=0
    f.write('#ms_Organism_id\tms_Number\tOrganism_name\tTaxon_id\n')

    print "Loading the list of organisms id and name..."
    ## Send the resquest to the MicroScope Web service
    getallPgdb=client.service.getAllPgdbs().pgdbVerySimpleVO
    print "Organisms loaded"

    ## Print the organism information in the output file
    for org in getallPgdb:
        organism_name=""
        taxon=""
        if "commonName" in getallPgdb[i].__keylist__:
            organism_name=str(getallPgdb[i].commonName)
        if "taxon" in getallPgdb[i].__keylist__:
            taxon=str(getallPgdb[i].taxon)
        line=str(getallPgdb[i].frameVersion)+'\t'+str(getallPgdb[i].id)+'\t'+organism_name+'\t'+taxon
        f.write(line+"\n")

        ## Display results if parser argument -p is specified
        if args.prt:
            print line

        i=i+1

    ## Close output file
    f.close()



##############################################################
### -------------------- MAIN PROGRAM -------------------- ###

if __name__ == '__main__':
    ## Read command line arguments
    parser = argparse.ArgumentParser(version='12.03',add_help=0)
    parser.add_argument("-h", "--help", action="store_true", dest="help", help='Show this help message and exit')
    parser.add_argument('-o', '--outputfile', action= "store", dest="file", help="Output file path+name",metavar='#', default='supported-organisms-genoscope.tab')
    parser.add_argument('-p', '--print', action= "store_true", dest="prt", help="Echo the result on the screen (in addition to the output file)")
    args = parser.parse_args()                           

    ## Function help in arguments parser
    if args.help:
        doc =  globals()['__doc__'] % {'version' : parser.version, 'progname' : parser.prog}
        print doc
        parser.print_help()
        sys.exit(0)

    ## Start calcul execution of programm  
    time1=time.time()

	## Create_file function call
    Create_file(args.file)

    ## Record end-time execution of programm
    time2=time.time()
    print "List of organisms exported in file", args.file, '\nin',"%.2f" % (time2-time1),'seconds'
    print "File created in", os.getcwd().replace('\\','/')

    ## This line serves only to let the DOS window open in Windows OS after execution of the program double clicking
#   raw_input("Enter to exit... ")
