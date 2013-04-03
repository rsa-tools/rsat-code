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
  Create a file that contains data returned by getGPR method: Genes list of
organism (different types of gene IDs), EC numbers of GPR in witch genes are
involved, taxonomy_id and species, from the frameVersion (ms_Organism_id).

OUTPUT FORMAT
  The program returns a tab-delimited text file with one row per
organism.

Columns:

1. Query
   Name or ID of genes.

2. (Reaction
   If -r argument is specified; RXN ID of reactions from GPR in
   witch genes are involved.)

3. EC
   EC numbers of GPR in witch genes are involved.

4. Qualifier
   Type of gene ID.

5. Name
   Names of genes.

6. Taxon_id
   Taxonomic ID (TAXID) as defined in the Taxonomy database at NCBI
   (http://www.ncbi.nlm.nih.gov/taxonomy/).

7. Organism_name
   Scientific name of the organism.

usage: gene-ec-genoscope.py -id # [-v] [-h] [-o #] [-r]

required argument:
  -id #, --ms_org_id #	Microscope organism ID e.g. 'ECOLI149-150'

optional arguments:
  -v, --version         show program's version number and exit
  -h, --help		Show this help message and exit
  -o #, --outputfile #	Name of the output file
  -r, --reactions	Add a column reaction with reactions of GPR in witch
			genes are involved
'''

##########################################################
### -------------------- PACKAGES -------------------- ###

import suds
import time
import argparse
import sys
import os
import re

###########################################################
### -------------------- FUNCTIONS -------------------- ###

def Connection():
	## Open a connection to the GenoScope Web services
	url = 'http://www.genoscope.cns.fr/microcyc-webservice/services/fr.genoscope.agc.microcyc.business.MicrocycService.wsdl'
	client = suds.client.Client(url, timeout=3600)
	return client



def getGPR(ms_Organism_id):
	client=Connection()

	print "Running getGPR on GenoScope Web service..."
	## Running getGPR for the microscope organism ID specified on GenoScope Web service
	try:
		getGPR=client.service.getGPR(ms_Organism_id)

	## Intercept suds.WebFault error if you enter a ID that is in the list of microscope organism IDs but is'nt available (only two are in this case) and return the list of available IDs
	except suds.WebFault:
		l=[]
		getallPgdb=client.service.getAllPgdbs().pgdbVerySimpleVO
		## Creation of the list of the available microscope organism ID
		for i in range(len(getallPgdb)):
			if "commonName" in getallPgdb[i].__keylist__ and "taxon" in getallPgdb[i].__keylist__:
				l.append(str(getallPgdb[i].frameVersion))
##		print "List of the available attributes:\n", l
		raise AttributeError, ms_Organism_id
	return getGPR



def FileName(getGPR):
	client=Connection()

	try:
		## File name built from the species name and the taxon
		fileName=str(getGPR.commonName)+"-"+str(getGPR.taxon)
		## Replace all is not an alphanumeric character by a "_" then two or more consecutive "_" by one
		regex = re.compile('\W')
		fileName=regex.sub("_",fileName)
		regex2=re.compile("_{2,}")
		fileName=regex2.sub("_",fileName)
		fileName=fileName+"-gene-ec-genoscope.tab"
		return fileName

	## Intercept attribute error if you enter an unavailable attribute and return the list of available attributes
	except AttributeError:
		l=[]
		getallPgdb=client.service.getAllPgdbs().pgdbVerySimpleVO
		## Creation of the list of the available microscope organism ID
		for i in range(len(getallPgdb)):
			if "commonName" in getallPgdb[i].__keylist__ and "taxon" in getallPgdb[i].__keylist__:
				l.append(str(getallPgdb[i].frameVersion))
		print "List of the available attributes:\n", l
		raise AttributeError, args.id



def Create_file(fileName, getGPR):
	if args.file!=None:
		fileName=args.file
	## List of the different types of gene id
	geneKeys=["accession1", "frame", "id", "proteinId", "refseqId"]

	## Print species name
	print "For", str(getGPR.commonName)+":"

	print "Writing data in file", fileName+"..."

	## Move in '../../../data/GPR' directory
#	os.chdir(r'../../data/GPR')

	if args.file==None:
		regex = re.compile('\W')
		## Directory name built from the species name and the taxon
		dirName=str(getGPR.commonName)+"-"+str(getGPR.taxon)
		dirName=regex.sub("_",dirName)
		regex2=re.compile("_{2,}")
		dirName=regex2.sub("_",dirName)
		## Create the dirName directory only if it does'nt exists
		if dirName not in os.listdir(os.getcwd()):
			os.mkdir(dirName)
		## Move in 'dirName' directory
		os.chdir(dirName)

	## Open the file "fileName" in write mode
	File=open(fileName, "w")

	File.write('; ' + str.join(" ", sys.argv)+'\n')
	File.write(';\n')
	File.write('; Time: '+time.strftime('%Y-%m-%d %H:%M\n',time.localtime()))
	File.write(';\n')


	## If -r argument is specified
	if args.reactions:
		## Print detail of column content. Not very useful for the time being since we only have two output fields, but can be useful in general when some files wil have many columns
		header="""; Column content
;\t1\tQuery
;\t2\tReaction
;\t3\tEC number
;\t4\tQualifier
;\t5\tName
;\t6\tTaxonomy_id
;\t7\tSpecies
;
"""
		File.write(header)

		## Write header line
		File.write("#Query\tReaction\tEC\tQualifier\tName\tTaxonomy_id\tSpecies\n")
		## Genes list and GPRs list
		genes=getGPR.genes.geneSimpleVO2
		gprs=getGPR.gprs.gprSimpleVO
		#print dir(genes[0])

		## len(list) = length of the "list" list = numbre of elements; range() fonction creates a list from an integer: range(10)=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
		## For each integer in list created
		for i in range(len(genes)):
			name=str(genes[i].shortName)
			if name==":error":
				name=""
			## Creation of a empty dictionary
			ReactionECDictionary={}
			## For each integer in list created
			for j in range(len(gprs)):
				## If "gpr" field (key) (gpr = list of frames of genes associated with GPR) of GPR number j is not set to None... 
				if gprs[j].gpr!=None:
					## Retrieve the frames list converted in string and split at each " and "
					frameList=str(gprs[j].gpr).split(" and ")
					## For each frame in frames list...
					for f in frameList:
						## If this frame is the same as the frame of gene number i
						if f==genes[i].frame:
							## Add Reaction (key) and EC (value) of the GPR number j to the Reaction-EC dictionary
							ReactionECDictionary[str(gprs[j].frameReaction)]=str(gprs[j].ecNumber)

			## If the Reaction-EC dictionary isn't left empty...
			if ReactionECDictionary!={}:
				itemList=ReactionECDictionary.items()
				## For each Reaction-EC tuple (key, value), write line in the file
				for frameReaction_EC in itemList:
					line=name+"\t"+frameReaction_EC[0]+"\t"+frameReaction_EC[1]+"\t"+"name"+"\t"+name+"\t"+str(getGPR.taxon)+"\t"+str(getGPR.commonName)+"\n"
					File.write(line)
					## For each type of gene ID...
					for k in geneKeys:
						## If this type is part of the keys list (print dir(genes[i]) to see different methods)
						if k in genes[i].__keylist__:
							kValue=str(genes[i][k])
							if kValue!="0":
								## One line by EC number and by type of gene ID: if there are 3 EC and 3 types, it creates 9 lines for one gene
								line=kValue+"\t"+frameReaction_EC[0]+"\t"+frameReaction_EC[1]+"\t"+k+"\t"+name+"\t"+str(getGPR.taxon)+"\t"+str(getGPR.commonName)+"\n"
								File.write(line)
			else:
				line=name+"\t\t\t"+"name"+"\t"+name+"\t"+str(getGPR.taxon)+"\t"+str(getGPR.commonName)+"\n"
				File.write(line)
				for k in geneKeys:
					if k in genes[i].__keylist__:
						kValue=str(genes[i][k])
						if kValue!="0":
							line=kValue+"\t\t\t"+k+"\t"+name+"\t"+str(getGPR.taxon)+"\t"+str(getGPR.commonName)+"\n"
							File.write(line)


	## If -r argument is not specified
	else:
		header="""; Column content
;\t1\tQuery
;\t2\tEC number
;\t3\tQualifier
;\t4\tName
;\t5\tTaxonomy_id
;\t6\tSpecies
;
"""
		File.write(header)

		## Write header line
		File.write("#Query\tEC\tQualifier\tName\tTaxonomy_id\tSpecies\n")
		## Genes list and GPRs list
		genes=getGPR.genes.geneSimpleVO2
		gprs=getGPR.gprs.gprSimpleVO
		#print dir(genes[0])

		## len(list) = length of the "list" list = numbre of elements; range() fonction creates a list from an integer: range(10)=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
		## For each integer in list created
		for i in range(len(genes)):
			name=str(genes[i].shortName)
			if name==":error":
				name=""

			## Creation of an empty list
			ECList=[]
			## For each integer in list created
			for j in range(len(gprs)):
				## If "gpr" field (key) (gpr = list of frames of genes associated with GPR) of GPR number j is not set to None... 
				if gprs[j].gpr!=None:
					## Retrieve the frames list converted in string and split at each " and "
					frameList=str(gprs[j].gpr).split(" and ")
					## For each frame in frames list...
					for f in frameList:
						## If this frame is the same as the frame of gene number i
						if f==genes[i].frame:
							ecNumber=str(gprs[j].ecNumber)
							## Prevent the storage of redundant EC number in the EC list
							if ecNumber not in ECList:
								## Add EC number of the GPR number j to the EC list
								ECList.append(ecNumber)

			## If the EC list isn't left empty...
			if ECList!=[]:
				## For each EC number, write line in the file
				for EC in ECList:
					line=name+"\t"+EC+"\t"+"name"+"\t"+name+"\t"+str(getGPR.taxon)+"\t"+str(getGPR.commonName)+"\n"
					File.write(line)
					## For each type of gene ID...
					for k in geneKeys:
						## If this type is part of the keys list (dir(genes[i]) to see different methods of genes[i])
						if k in genes[i].__keylist__:
							kValue=str(genes[i][k])
							if kValue!="0":
								## One line by EC number and one line by type of gene ID: if there are 3 EC and 3 types, it creates 9 lines for one gene
								line=kValue+"\t"+EC+"\t"+k+"\t"+name+"\t"+str(getGPR.taxon)+"\t"+str(getGPR.commonName)+"\n"
								File.write(line)
			else:
				line=name+"\t\t"+"name"+"\t"+name+"\t"+str(getGPR.taxon)+"\t"+str(getGPR.commonName)+"\n"
				File.write(line)
				for k in geneKeys:
					if k in genes[i].__keylist__:
						kValue=str(genes[i][k])
						if kValue!="0":
							line=kValue+"\t\t"+k+"\t"+name+"\t"+str(getGPR.taxon)+"\t"+str(getGPR.commonName)+"\n"
							File.write(line)
	## Close file
	File.close



##############################################################
### -------------------- MAIN PROGRAM -------------------- ###

if __name__=="__main__":

	## Read command line arguments
	parser = argparse.ArgumentParser(version='12.03',add_help=0)
	parser.add_argument('-id', '--ms_org_id', action= "store", dest="id")
	parser.add_argument('-o', '--outputfile', action= "store", dest="file")
	parser.add_argument("-h", "--help", action="store_true", dest="help")
	parser.add_argument("-r", "--reactions", action="store_true", dest="reactions")
	args = parser.parse_args()

	## Required id argument
	if args.id==None:
		doc =  globals()['__doc__'] % {'version' : parser.version, 'progname' : parser.prog}
		print doc
		sys.exit(0)

	## Function help in arguments parser
	if args.help:
		doc =  globals()['__doc__'] % {'version' : parser.version, 'progname' : parser.prog}
		print doc
		sys.exit(0)

	## Start time calcul execution of programm
	time1=time.time()

	## Necessary functions called
	getGPR=getGPR(args.id)
	fileName=FileName(getGPR)
	Create_file(fileName, getGPR)

	## Record end-time execution of programm
	time2=time.time()

	print 'Writing completed in',"%.2f" % (time2-time1),'seconds'
	print "File created in", os.getcwd().replace('\\','/')

