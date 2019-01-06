#!/usr/bin/python
#-*-coding: utf-8-*-

"""This program contains shared functions for a series of REST client
programs enabling to collect data from various bioinformatics
resources (Ensembl, Ensembl Genomes, Genoscope and MICROME

Currently, the concerned programs using workflows_functions are :
 - supported_organisms_ensembl.py
 - supported_organisms_ensemblgenomes.py
 - supported_organisms_genoscope.py
 - supported_organisms_microme.py
 - gene2reaction_from_genoscope.py

"""

import httplib2 #To send the request using the http protocol
import sys #To quit the program or display a error message on the standard error

def RequestGeneric(url_root,url_extension):
	"""The function RequestGeneric() sends a request to a server. 
	
	@type	url_root : string
	@param	url_root : server URL
	@type	url_extension : string
	@param	url_extension : extension of server URL, specifying the request
	@returns : this function returns the JSON object obtained from a server
	"""
	
	## Note by JvH (2014-09-31): this cache directory should be
	## cleaned after download, because it can take a large space
	## (I realized it occupied 3.6Gb after having downloaded some
	## bacteria).
	url_root = url_root
	http = httplib2.Http(".cache")
	
	## Send the request to the server
	resp, response = http.request(url_root + url_extension, method = "GET", headers = {"Content-Type" : "application/json"})
	
	## Catch errors
	if not resp.status == 200:
		print( "Invalid response: ", resp.status)
		sys.exit()
	
	## Return the result : here a Json object
	return response



def RequestGenoscope(url_extension):
	"""The function RequestGenoscope() sends a request to the Genoscope server. 
	
	@type	url_extension : string
	@param	url_extension : extension of server URL, specifying the request
	@returns : this function returns the JSON object obtained from the Genoscope REST project
	"""
	
	return RequestGeneric("http://mgalileo.genoscope.cns.fr/api/microscope",url_extension)
	


def RequestMicrome(url_extension):
	"""The function RequestMicrome() sends a request to the MICROME server. 
	
	@type	url_extension : string
	@param	url_extension : extension of server URL, specifying the request
	@returns : this function returns the JSON object obtained from the MICROME REST project
	"""
	
	return RequestGeneric("http://www.microme.eu/reaction_matrix/api",url_extension)
	

def RequestEnsembl(url_extension):
	"""The function RequestEnsembl sends a request to the Ensembl server. 
	
	@type	url_extension : string
	@param	url_extension : extension of server URL, specifying the request
	@returns : this function returns the JSON object obtained from the Ensembl REST project
	"""
	
	return RequestGeneric("http://beta.rest.ensembl.org",url_extension)
	

def RequestEnsemblGenomes(url_extension):
	"""This function send a request to the Ensembl Genomes server. 
	
	@type	url_extension : string
	@param	url_extension : extension of server URL, specifying the request
	@returns : this function returns the JSON object obtained from the Ensembl Genomes REST project
	"""
	
	return RequestGeneric("http://beta.rest.ensemblgenomes.org",url_extension)
	

def Test_Of_Chosen_Out_Columns(arg_list,dict_keys):
	"""This function test the existence of a column argument chosen by the user.
	
	@type	arg_list : list
	@param	arg_list : list of columns arguments available 
	@type	dict_keys : list
	@param	dict_keys : list of keys contained in the dictionary giving informations to complete the output content 
	@returns : this function displays a message error on the standard error and the program is killed if a column argument doesn't exist
	"""
	
	#This list contains not available arguments
	uncorrect_argument = []
	for arg in arg_list:
		if arg not in dict_keys:
			uncorrect_argument.append(arg)
	if len(uncorrect_argument) != 0:
		for uncorrect in uncorrect_argument:
			print >> sys.stderr, 'This chosen argument is not available :', uncorrect
		print >> sys.stderr, '\nEnd of the operation'
		sys.exit() 
		

def Creation_Of_The_Tabulated_File(arg_list,file_name,file_content):
	"""This function creates a tabulated file containing informations got
        from the used server.
	
	@type	arg_list : list
	@param	arg_list : list of columns arguments available 
	@type	file_name : string
	@param	file_name : name given to the tabulated file by the user
	@type	file_content : string
	@param	file_content : informations presented as a tabulated text. Refer to Complete_Output_Content
	@returns : this function creates the tabulated file

        """
	
	f = open(file_name,'w')
	f.write('#' + arg_list[0] + '\t' + '\t'.join(arg_list[1:]) + '\n')
	f.write(file_content.encode('utf-8'))
	f.close()
	
def Display_Tabulated_Results(arg_list,output):
	"""This function displays tabulated results containing informations got from the used server
	
	@type	arg_list : list
	@param	arg_list : list of columns arguments available 
	@type	output : string
	@param	output : informations presented as a tabulated text.
	@returns : this function displays the result on the standard output 
	"""
	
	print( '#' + arg_list[0] + '\t' + '\t'.join(arg_list[1:]))
	print( output)
	

def Save_Tabulated_File(list_columnHeaders,outfile,file_content):
	"""This function creates a tabulated file containing informations got from the used server
	
	@type	arg_list : list
	@param	arg_list : list of columns arguments available 
	@type	file_name : string
	@param	file_content : informations presented as a tabulated text. Refer to Complete_Output_Content
	@returns : this function creates the tabulated file
	"""

	f = open(outfile,'w')
	f.write('#' + list_columnHeaders[0] + '\t' + '\t'.join(list_columnHeaders[1:]) + '\n')
	f.write(file_content.encode('utf-8'))
	f.close()

      
