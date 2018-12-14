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
  Create a file that contains data returned by getAllReactions method:
reactions and compounds list. Determines also the list of arcs.

OUTPUT FORMAT
  The program returns a tab-delimited text file with one row per
node.

Columns:

1. NODES
   Normal/reverse reactions (RXN ID) and compounds names.

2. Object type
   Node type: reaction or compound.

3. Exclusion attribute
   Reaction or compound excluded during the reconstruction of metabolic
   network.

4. Label
   Compounds.

8. FORMULA

+ List of arcs: substrate - (normal or reverse) reaction
                (normal or reverse) reaction - product

If -n argument is specified:
Create a file that contains the list of nodes: reactions (RXN IDs, EC numbers),
compounds (scientist names, common names, CPD IDs).

1. query
   Reactions (RXN IDs, EC numbers).
   Compounds (scientist names, common names, CPD IDs).

2. id
   RXN IDs, compounds names.

3. qualifier
   Type of node.

4. name
   Reactions, compounds.
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

def reactions():
	print "Running getAllReactions on GenoScope Web service..."
	## Open a connection to the GenoScope Web services
	url='http://www.genoscope.cns.fr/microcyc-webservice/services/fr.genoscope.agc.microcyc.business.MicrocycService.wsdl'
	client=suds.client.Client(url, timeout=3600)
	## Running getAllReactions on GenoScope Web service
	reactions=client.service.getAllReactions().reactionSimpleVO
	return reactions



def Create_file(fileName, reactions):
	print "Running getCompoundSimple on GWs for each compound (this may take a long time)\nand writing data in file", fileName+"..."

	## Open the file "fileName" in write mode
	File=open(fileName, "w")

	File.write(';' + str.join(" ", sys.argv)+'\n')
	File.write(';\n')
	File.write(';Time: '+time.strftime('%Y-%m-%d %H:%M\n',time.localtime()))
	File.write(';\n')


	## If -n argument is specified
	if args.nodes:
		## Print detail of column content. Not very useful for the time being since we only have two output fields, but can be useful in general when some files wil have many columns
		header=""";Column content
;\t1\tquery
;\t2\tid
;\t3\tqualifier
;\t4\tname
;
"""
		File.write(header)

		## Write header line
		File.write(";query\tid\tqualifier\tname\n")

		url='http://www.genoscope.cns.fr/microcyc-webservice/services/fr.genoscope.agc.microcyc.business.MicrocycService.wsdl'
		client=suds.client.Client(url, timeout=3600)

		## Regular expression matching html tag
		htmlTag = re.compile('<[^>]*>|&|;')

		## Creation of an empty list
		compoundsList=[]
		## Create a list from an integer (number of reactions): range(10)=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
		for i in range(len(reactions)):
			## Retrieve frame, EC number and common name of reaction number i
			frame=str(reactions[i].frame)
			ecNumber=str(reactions[i].ecNumber)
			commonName=str(reactions[i].commonName)
			## If commonName and ecNumber are worth "NIL", write one line for frame
			if commonName=="NIL" and ecNumber=="NIL":
				line=frame+"\t"+frame+"\treactions\tNA\n"
			## If commonName is worth "NIL", write one line for frame and one for EC number
			elif commonName=="NIL":
				line=frame+"\t"+frame+"\treactions\tNA\n"+ecNumber+"\t"+frame+"\tEC\tNA\n"
			## If ecNumber is worth "NIL", write one line for frame and one for common name
			elif ecNumber=="NIL":
				line=frame+"\t"+frame+"\treactions\t"+commonName+"\n"+commonName+"\t"+frame+"\treactions\t"+commonName+"\n"
			## Else write one line for frame, one for common name and one for EC number
			else:
				commonName=htmlTag.sub("", commonName)
				line=frame+"\t"+frame+"\treactions\t"+commonName+"\n"+commonName+"\t"+frame+"\treactions\t"+commonName+"\n"+ecNumber+"\t"+frame+"\tEC\t"+commonName+"\n"
			File.write(line)

			## If there are left compounds...
			if reactions[i].leftCompounds!=None:
				leftCompounds=reactions[i].leftCompounds.reactionParticipantSimpleVO
				for lC in leftCompounds:
					## Retrieve compoundFrame of lC
					compoundFrame=str(lC.compoundFrame)
					## Remove "|" character if it is in compoundFrame
					compoundFrame=compoundFrame.replace("|", "")
					## If compoundFrame is not in compoundsList...
					if compoundFrame not in compoundsList:
						## Add the compound to the list of compounds
						compoundsList.append(compoundFrame)
						## Running getCompoundSimple for the lC compoundId on GenoScope Web service
						CompoundSimple=client.service.getCompoundSimple(lC.compoundId)
						## If "commonName" attribute is part of the keys list (print dir(CompoundSimple) to see different methods of CompoundSimple)
						if "commonName" in CompoundSimple.__keylist__:
							compoundName=str(CompoundSimple.commonName)
							## Replace htmlTag by nothing
							compoundName=htmlTag.sub("", compoundName)
							## If compoundFrame is different from compoundName, write two lines, else write one (prevent the repetition of two identical rows)
							if compoundFrame!=compoundName:
								line=compoundFrame+"\t"+compoundFrame+"\tcompounds\t"+compoundName+"\n"+compoundName+"\t"+compoundFrame+"\tcompounds\t"+compoundName+"\n"
							else:
								line=compoundFrame+"\t"+compoundFrame+"\tcompounds\t"+compoundFrame+"\n"
						## Else write one line without common name
						else:
							line=compoundFrame+"\t"+compoundFrame+"\tcompounds\tNA\n"
						File.write(line)

			## If there are right compounds
			if reactions[i].rightCompounds!=None:
				rightCompounds=reactions[i].rightCompounds.reactionParticipantSimpleVO
				for rC in rightCompounds:
					compoundFrame=str(rC.compoundFrame)
					compoundFrame=compoundFrame.replace("|", "")
					if compoundFrame not in compoundsList:
						compoundsList.append(compoundFrame)
						CompoundSimple=client.service.getCompoundSimple(rC.compoundId)
						if "commonName" in CompoundSimple.__keylist__:
							compoundName=str(CompoundSimple.commonName)
							compoundName=htmlTag.sub("", compoundName)
							if compoundName!=compoundFrame:
								line=compoundFrame+"\t"+compoundFrame+"\tcompounds\t"+compoundName+"\n"+compoundName+"\t"+compoundFrame+"\tcompounds\t"+compoundName+"\n"
							else:
								line=compoundFrame+"\t"+compoundFrame+"\tcompounds\t"+compoundFrame+"\n"
						else:
							line=compoundFrame+"\t"+compoundFrame+"\tcompounds\tNA\n"
						File.write(line)


	## If -n argument is not specified
	else:
		header=""";Column content
;\t1\tNODES
;\t2\tObject type
;\t3\tExclusion attribute
;\t4\tLabel
;\t5\tPubChem
;\t6\tCAS
;\t7\tKEGG_identifier
;\t8\tFORMULA
;
"""
		File.write(header)

		## Write header line
		File.write(";NODES\tObject type\tExclusion attribute\tLabel\tPubChem\tCAS\tKEGG_identifier\tFORMULA\n")

		url='http://www.genoscope.cns.fr/microcyc-webservice/services/fr.genoscope.agc.microcyc.business.MicrocycService.wsdl'
		client=suds.client.Client(url)

		htmlTag = re.compile('<[^>]*>|&|;')

		compoundsList=[]
		arcs=";ARCS\n"
		## Create a list from an integer: range(10)=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
		for i in range(len(reactions)):
			reactionFrame=str(reactions[i].frame)
			## Write two lines by reaction (normal and reverse)
			line=reactionFrame+">\tReaction\t"+reactionFrame+"\tNA\tNA\tNA\tNA\tNA\n"+reactionFrame+"<\tReaction\t"+reactionFrame+"\tNA\tNA\tNA\tNA\tNA\n"
			File.write(line)

			if reactions[i].leftCompounds!=None:
				leftCompounds=reactions[i].leftCompounds.reactionParticipantSimpleVO
				for lC in leftCompounds:
					compoundFrame=str(lC.compoundFrame)
					compoundFrame=compoundFrame.replace("|", "")
					## Two arcs by reaction and compound (normal and reverse)
					arc=compoundFrame+"\t"+reactionFrame+">\n"+reactionFrame+"<\t"+compoundFrame+"\n"
					arcs+=arc
					if compoundFrame not in compoundsList:
						compoundsList.append(compoundFrame)
						CompoundSimple=client.service.getCompoundSimple(lC.compoundId)
						if "commonName" in CompoundSimple.__keylist__:
							labelCompound=str(CompoundSimple.commonName)
							labelCompound=htmlTag.sub("", labelCompound)
						else:
							labelCompound="NA"
						## If "inchi" attribute is part of the keys list
						if "inchi" in CompoundSimple.__keylist__:
							inchi=str(CompoundSimple.inchi)
							if inchi=="NIL":
								inchi="NA"
						else:
							inchi="NA"
						line=compoundFrame+"\tCompound\t"+compoundFrame+"\t"+labelCompound+"\tNA\tNA\tNA\t"+inchi+"\n"
						File.write(line)

			if reactions[i].rightCompounds!=None:
				rightCompounds=reactions[i].rightCompounds.reactionParticipantSimpleVO
				for rC in rightCompounds:
					compoundFrame=str(rC.compoundFrame)
					compoundFrame=compoundFrame.replace("|", "")
					arc=reactionFrame+">\t"+compoundFrame+"\n"+compoundFrame+"\t"+reactionFrame+"<\n"
					arcs+=arc
					if compoundFrame not in compoundsList:
						compoundsList.append(compoundFrame)
						CompoundSimple=client.service.getCompoundSimple(rC.compoundId)
						if "commonName" in CompoundSimple.__keylist__:
							labelCompound=str(CompoundSimple.commonName)
							labelCompound=htmlTag.sub("", labelCompound)
						else:
							labelCompound="NA"
						if "inchi" in CompoundSimple.__keylist__:
							inchi=str(CompoundSimple.inchi)
							if inchi=="NIL":
								inchi="NA"
						else:
							inchi="NA"
						line=compoundFrame+"\tCompound\t"+compoundFrame+"\t"+labelCompound+"\tNA\tNA\tNA\t"+inchi+"\n"
						File.write(line)
#		print compoundsList
		## Write arcs at the end of file 
		File.write(arcs)
	## Close file
	File.close



##############################################################
### -------------------- MAIN PROGRAM -------------------- ###

if __name__=="__main__":
	## Read command line arguments
	parser = argparse.ArgumentParser(version='12.04',add_help=0)
	parser.add_argument("-h", "--help", action="store_true", dest="help", help="Show this help message and exit")
	parser.add_argument('-o', '--outputfile', action= "store", dest="file", metavar='#', default="MICROME_v15_directed-metab-network.tab", help="Name of the output file")
	parser.add_argument("-n", "--nodes", action="store_true", dest="nodes", help="Create a file that contain nodes (reactions and compounds with their ID: EC number, RXN and CPD)")
	args = parser.parse_args()

	## Function help in arguments parser
	if args.help:
		doc =  globals()['__doc__'] % {'version' : parser.version, 'progname' : parser.prog}
		print doc
		parser.print_help()
		sys.exit(0)

	time1=time.time()

	## Necessary functions called
	reactions=reactions()
	Create_file(args.file, reactions)

	## Record end-time execution of programm
	time2=time.time()

	print 'Writing completed in',"%.2f" % (time2-time1),'seconds'
	print "File created in", os.getcwd().replace('\\','/')

