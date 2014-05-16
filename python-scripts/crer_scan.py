#!/usr/bin/python3
# -*- coding: utf-8 -*-

#===============================================================================
# Documentation
#===============================================================================

__doc__ ="""

Scanning of predicted sites on sequence.
And Detection of putative cis-regulatory enriched regions (CRERs).

usage: 
	   crer_scan.py [-h] 
					[-i FILE]
					[-o OUTFILE]
					[-v VERBOSE]
					[-in_format FORMAT]
					[-s]
					[-return_limits]
					[-return_limits_filtered] 
					[-uth_site_pval UTH_SITE_PVAL] 
					[-lth_score LTH_SITE_SCORE]
					[-uth_score UTH_SITE_SCORE]
					[-lth_crer_size LTH_CRER_SIZE]
					[-uth_crer_size UTH_CRER_SIZE]
					[-lth_crer_sites LTH_CRER_SITES]
					[-uth_crer_sites UTH_CRER_SITES]
					[-lth_crer_sites_distance LTH_CRER_SITES_DISTANCE]
					[-uth_crer_sites_distance UTH_CRER_SITES_DISTANCE]
					[-uth_crer_pval UTH_CRER_PVAL]
					[-uth_crer_eval UTH_CRER_EVAL]
					[-lth_crer_sig LTH_CRER_SIGNIFICANCE]
					[-uth_overlap UTH_OVERLAP]
					[-nopval]
					[-pre_table]

Optional arguments:

  -h, --help
		show this help message and exit
		
  -i FILE
		If not specified, input is read from STDIN
		
  -in_format FORMAT
		input_format. 
		Default: ft (produced by RSAT matrix-scan and dna-pattern).
		Supported: ft, bed
		
  -o OUTFILE
		  Output file in ft format
		  
  -v VERBOSE
		  level of verbose.
		  Messages are wrote on standard error.
		  Supported: Integer = 1,2,3
		  By default : 1 = No message
		  Level 2 : moderately density of messages
		  Level 3 : High density
		  
  -s
		Sort the list of sites.
		Very recommended.
		The sites are sorted by center position

  -return_limits
		  return every limits of sequences
		  By default : no return any limits
		  
  -return_limits_filtered
		return the limits filtered of the sequence. 
		Only the sequence limits of CRERs.
		By default : no return any limits
		
  -uth_site_pval UTH_SITE_PVAL
		maximal p-value of sites to be considered
		Recommended to be the higher site p-value considered.
		Default = 1e-4
		
  -lth_score LTH_SITE_SCORE
		minimal site score to be considered
		
  -uth_score UTH_SITE_SCORE
		maximal site score to be considered

  -lth_crer_size LTH_CRER_SIZE
		minimal size of the enriched region (in bp). 
		Default: minimal site size = 30bp
		
  -uth_crer_size UTH_CRER_SIZE
		maximal size of the enriched region (in bp). 
		Default: maximal site site = 500bp
		
  -lth_crer_sites LTH_CRER_SITES
		minimal number of sites covered by the enriched region. 
		Default: minimal number of sites = 2
		
  -uth_crer_sites UTH_CRER_SITES
		maximal number of sites covered by the enriched region

  -lth_crer_sites_distance LTH_CRER_SITES_DISTANCE
		distance between successive sites to be considered. 
		A minimal inter-site distance can be used to prevent overlap between redundant matrices. 
		Default = minimal distance = 1 bp
		
  -uth_crer_sites_distance UTH_CRER_SITES_DISTANCE
		distance between successive sites to be considered. 
		A maximal inter-site distance can be used to prevent merging distinct modules into a single one. 
		Note: the maximal inter-site distance is one of the most influential parameters in cluster-buster. 
		Default = maximal distance = 35 bp
		
  -uth_crer_pval UTH_CRER_PVAL
		maximal binomial p-value
		Default: 1e-4
		
  -uth_crer_eval UTH_CRER_EVAL
		maximal e-value
		Default: 1e-4
		
  -lth_crer_sig LTH_CRER_SIGNIFICANCE
		minimal binomial significance
		Default: 2
		
  -uth_overlap UTH_OVERLAP
		maximal overlap to define two distinct sites
		
  -nopval
		  compute crer without p value
		  
  -pre_table
		  compute a table where is all possible p_value
		  Useful where there is a huge number of sites to scan.
"""

__author__ = """
Marie Artufel, Lucie Khamvongsa-Charbonnier *


* in alphabetic order
"""

__date__ = "18 April 2014"

__credits__ = """
Mr Van-Helden for his support
Mr Gonzalez for his advices
Mr Tichit for his psychological support at the breakfast time
And Genopole class room for its liberality
"""

"""
Us to be the Best
Mr Jack for his support 
Mr Gonzalez for his legendary A+ in the bus
Genopole class room for its liberality
and last not least Mr Tichit for nothing
"""

__version__ = 'v.1.00 '


#===============================================================================
# Modules
#===============================================================================

import argparse
import sys
import tempfile
import re

import time
import datetime

import scipy.stats
import math


#===============================================================================
# Classes
#===============================================================================

# Creation of site's class
class Site(object):
	"""
	This class is intended for sites of transcription factors identified by the matrix-scan program of the Regulatory Sequence Analysis Tools (RSAT).
	
	The attributes of this class are : 
		seq_id -> the identity of the sequence
		ft_name -> the name/identity of the transcription factor
		strand -> the strand of DNA where is the site of transcription factor
		start -> the beginning of the site given by the number of the first nucleotide
		center -> the center of the site given by the number of the first nucleotide
		end -> the end of the given site by the number of the last nucleotide
		weight -> the weight or score of the given site
	"""
	
	
	def __init__(self,seq_id,ft_name,strand,start,end,weight):
		"""
		Class initiation 
		
		The arguments are:
			seq_id -> the identity of the sequence. Type = string
			ft_name -> the name/identity of the transcription factor. Type = string
			strand -> the strand of DNA where is the site of transcription factor. Type = string
			start -> the beginning of the site given by the number of the first nucleotide. Type = string, float,int
			center -> the center of the site given by the number of the first nucleotide. Type = string, float, int
			end -> the end of the given site by the number of the last nucleotide. Type = string, float, int
			weight -> the weight or score of the given site. Type = string, float, int
		"""
			
		self.__seq_id = seq_id
		self.__ft_name = ft_name
		self.__strand = strand
		self.__start = float(start)
		self.__center = int((float(end)+float(start))/2)
		self.__end = float(end)
		self.__weight = float(weight)
	
	
	def seq_id(self):
		""" Return the identity of the sequence and do not take arguments"""
		
		return self.__seq_id
	
	
	def ft_name(self):
		"""Return the name of the transcription factor and do not take arguments"""
		
		return self.__ft_name
	
	
	def strand(self):
		"""Return the name of the strand where the site is and do not take arguments"""
		
		return self.__strand
	
	
	def start(self):
		"""Return site's start and do not take arguments"""
		
		return self.__start
	
	
	def end(self):
		"""Return site's end and do not take arguments"""
		
		return self.__end
	
	
	def center(self):
		"""Return site's center and do not take arguments"""
		
		return self.__center
	
	
	def weight(self):
		"""Return site's weight and do not take arguments"""
		
		return self.__weight


# Creation of subclass to treat features file		
class Siteft(Site):
	"""
	This subclass is intended for sites of transcription factors identified by the matrix scan software from feature's files
	because depending on the format of the input file we don't have the same information
	the bed format contains less information than the format feature (. ft.)
	
	The attributes of this class are those of the class site and the following attributes : 
		seq_id -> the identity of the sequence
		ft_name -> the name/identity of the transcription factor
		strand -> the strand of DNA where is the site of transcription factor
		start -> the beginning of the site given by the number of the first nucleotide
		center -> the center of the site given by the number of the first nucleotide
		end -> the end of the given site by the number of the last nucleotide
		weight -> the weight or score of the given site
		sequence -> the sequence of the site of the transcription factor 
		Pval -> the P-value of this site identified by the matrix scan software
		ln_Pval -> the natural logarithm of the p-value 
		sig -> the binomial significance of the site
	"""
	
	
	def __init__(self,seq_id,ft_name,strand,start,end,sequence,weight,Pval,ln_Pval,sig):
		""" Class initiation.
		
		The arguments are:
		Class initiation 
		
		The arguments are:
			seq_id -> the identity of the sequence. Type = string
			ft_name -> the name/identity of the transcription factor. Type = string
			strand -> the strand of DNA where is the site of transcription factor. Type = string
			start -> the beginning of the site given by the number of the first nucleotide. Type = string, float,int
			center -> the center of the site given by the number of the first nucleotide. Type = string, float, int
			end -> the end of the given site by the number of the last nucleotide. Type = string, float, int
			weight -> the weight or score of the given site. Type = string, float, int
			sequence -> the sequence of the site of the transcription factor. Type = str
			Pval -> the P-value of this site identified by the matrix scan software. Type = string, float, int
			ln_Pval -> the natural logarithm of the p-value. Type = string, float, int
			sig -> the binomial significance of the site. Type = string, float, int
		"""
		
		Site.__init__(self,seq_id,ft_name,strand,start,end,weight)
		self.__sequence = sequence
		self.__Pval = float(Pval)
		self.__ln_Pval = float(ln_Pval)
		self.__sig = float(sig)
	
	
	def sequence(self):
		"""Return site's sequence and do not take arguments"""
		
		return self.__sequence
	
	
	def Pval(self):
		"""Return site's P-value and do not take arguments"""
		
		return self.__Pval
	
	
	def ln_Pval(self):
		"""Return logarithm site's P-value and do not take arguments"""
		
		return self.__ln_Pval
	
	
	def sig(self):
		"""Return site's significance and do not take arguments"""
		
		return self.__sig
	
	
# Creation of class to treat sites on the same sequence
class SiteSetOnSameSeq(object):
	"""
	This class relates to a site collection on the same sequence as you can analyze several sites from different sequences at the same time
	
	The attributes of this class are : 
		seq_id -> the identity of the sequence
		sites_list -> the list of sites on the same sequence
		sorted -> this argument is false if the list of sites on a same sequence is not sorted and if it is this argument is true
	"""
	
	
	def __init__(self,seq_id, sites_list, sorted):
		"""
		Class initiation
		
		The arguments are:
			seq_id -> the identity of the sequence. Type = string
			sites_list -> the list of sites on the same sequence. Type = list of Site or Siteft object
			sorted -> this argument is false if the list of sites on a same sequence is not sorted and if it is this argument is true. Type = Boolean
		"""
		
		self.__seq_id = seq_id
		self.__sites_list = sites_list
		self.__sorted = sorted
	
	
	def seq_id(self):
		"""Return the identity of the sequence and do not take arguments"""
		
		return self.__seq_id
	
	
	def get_sites(self):
		"""Return a list of sites who are in the same sequence and do not take arguments"""
		
		return self.__sites_list
	
	
	def sorted(self):
		"""Return if the list is sorted or not and do not take arguments"""
		
		return self.__sorted
	
	
	def sort(self):
		"""
		This method sorts the sites according to their attributes 'start' thanks to a local function which is designated 'lambda'
		once the list sorted, 'sorted' attribute is changed. 
		Do not take arguments.
		"""
		
		self.__sites_list = sorted(self.__sites_list,key=lambda x: x.center()) #this function sorts by center's position
		self.__sorted = True
	
	
	def add_site_same_seq(self,site):
		"""	print()
		This method adds a site to the list of site located on the same sequence. 
		For this, we verify that the site is well to add the sequence concerned.
		
		Take one argument:
			site -> Site object to add at the SiteSetOnSameSeq object. Type = Site or Siteft object
		"""
		
		# Verification of sequence identity and add the site to the list
		if self.seq_id() == site.seq_id():
			
			self.__sites_list.append(site)
	
	
	def get_sorted_sites(self):
		"""
		This method allows us to obtain the list of sites sorts.
		Do not take arguments.
		"""
		
		# Use the sort method defined earlier and return the sorted list
		self.sort()
		return self.__sites_list
	
	
# Creation of class to treat a list of SiteSetOnSameSeq
class SiteSet(object):
	"""
	This class relates to a collection of several collections of sites on the same sequence.
	Example, if we consider three sequences at the same time we have a collection containing three collections
	that corresponds to each of the sequences.
	
	The attributes of this class are : 
		site_set_list -> This is a list containing the different lists where each corresponds to a sequence
		sorted -> this argument is false if the list of sites on a same sequence is not sorted and if it is this argument is true
	"""
	
	def __init__(self,site_set_list, sorted):
		"""
		Class initiation
		
		The arguments are:
			site_set_list -> This is a list containing the different lists where each corresponds to a sequence. Type = list of SiteSetOnSameSeq object
			sorted -> this argument is false if the list of sites on a same sequence is not sorted and if it is this argument is true. Type: Boolean
		"""
		
		self.__site_set_list = site_set_list
		self.__sorted = sorted
	
	
	def site_set_list(self):
		"""Return  the list of site set. Do not take arguments."""
		
		return self.__site_set_list
	
	
	def sorted(self):
		"""Return if the collection is sorted or not. Do not take arguments."""
		
		return self.__sorted
	
	
	def sort(self):
		"""This method sorts the sites by the sort method of SiteSetOnSameSeq class. Do not take arguments."""
		
		# For every SiteSetOnSameSeq object in the list, use the method sort defined in that class.
		for i in self.__site_set_list:
			i.sort()
		self.__sorted = True
	
	
	def add_site(self,site):
		"""
		This method allows us to add a site in the collection of site collections.
		This is done in several stages.
		At first we check if the collection is empty collections. In that case we add the site object by instantiating the class SiteSetOnSameSeq.
		If the collection of collections is not empty, we check if there is a corresponding SiteSetOnSameSeq object 
		then the site is added to the right SiteSetOnSameSeq object.
		If the good object is not found we create a new object of the class SiteSetOnSameSeq for it.
		
		The Argument is:
			site -> Site object to add. Type = Site or Siteft object
		"""
		
		# If the object contain an empty site list, the site is transform as an SiteSetOnSameSeq object and added to the liste
		if self.__site_set_list == []:
			self.__site_set_list.append(SiteSetOnSameSeq(site.seq_id(),[site],False))
		
		else:
			
			add = False
			# Identify if there is always an SiteSetOnSame object in list which correspond to the site and add this.
			for y in self.site_set_list():
			
				if y.seq_id() == site.seq_id():
					y.add_site_same_seq(site)
					add = True
			
			# If the site do not added create a new SiteSetOnSameSeq object with this site and this to the list (site_set_list)
			if not add : 
				self.__site_set_list.append(SiteSetOnSameSeq(site.seq_id(),[site],False))
	
	
	def get_sites(self):
		"""Return a list of all the sites studied. Do not take arguments."""
		
		# Creation of a list(get_sites) with every sites contained in site_set_list
		get_sites=[]
		
		for SiteSetOnSameSeq in self.site_set_list():
			list_siteset=SiteSetOnSameSeq.get_sites()
			
			for sites in list_siteset:
				get_sites.append(sites)
		
		return get_sites
	
	
	def get_sorted_sites(self):
		"""Return a sorted list of all the sites studied. Do not take arguments."""
		
		# Using get_sites method to obtain every sites and sort this sites by center position
		list_sites=self.get_sites()
		get_sorted_sites=sorted(list_sites,key=lambda x: x.center())
		
		return get_sorted_sites


#===============================================================================
# Functions
#===============================================================================

def read_features(file_name,pval_threshold,lth_score_threshold, uth_score_threshold, verbose):
	""" 
	Read a matrix scan's output file.
	Return respectively,
		- a list of site objects 
		- a list of sequence limits which are list with the sequence,start position and end position of the sequence
		- the number of matrix used to predict binding sites
	
	This function will take as arguments:
		file_name -> name of file in feature format. Type = string
		pval_threshold -> maximal p-value for read sites. Type = float
		lth_score_threshold -> minimal score for read sites. Type = float
		uth_score_threshold -> maximal score for read sites. Type = float
		verbose -> level of verbose. Type = int
	"""
	
	# Verbose message
	if verbose >= 3 :
		time_warn("Starting to parse")
	
	# opening and reading of file
	filehandle = open(file_name,'r')
	read = filehandle.readlines() # list of lines
	tab_list=list(set(read)) # create a list without duplicates

	# Initiation of variables
	sites_list =[]
	i = 1
	list_limit =[]
	j = 0
	object_list=[]
	max_site_pvalue = 0
	
	# treatment line by line and extraction of interesting information in a list
	for line in tab_list :
		
		# Extraction of tab information
		if line[0] != ";" and line[0]!="#":
			table = line.split()
			
			# Extraction of limit sequence informations in a list (init_table) 
			# which is added to another list (list_limit)
			if table[1]== 'limit':
				
				init_table = []
				init_table.append(table[0])
				init_table.append(table[4])
				init_table.append(table[5])
				list_limit.append(init_table)
			
			# Extraction of sites informations
			else:
				
				del table[1] # Remove the second element of the table ('site')
				table.insert(0,i) # Insertion of an index
				sites_list.append(table) # Add this list of table elements (site's informations) to the list of sites
				i = i+1
	
	# Extraction of other information like the number of matrix
	for line in read :
		
		# When there is the word "Matches" the followed informations are extracted
		if line[2:9]=="Matches":
			matrix = read[j+2:] 
			matrix_list = []
			
			# while the word "TOTAL" doesn't found the number find (second element of the row) are added in a list (matrix_list)
			for row in matrix:
				
				if row[3:8] == "TOTAL" :
					break
				
				else:
					num = row.split()
					matrix_list.append(int(num[1]))
					
		j= j+1
	
	# Number of matrix is the maximum of the matrix_list
	number_of_matrix = max(matrix_list)
	
	
	if verbose >=3 :
		time_warn("	Information extracted from the input")
	
	
	# every site is converted in SiteFt class object
	for j in sites_list :
		
		# Replace "." row in 0
		if j[7] == ".":
			j[7] = 0
		
		# If the site's p-value and score are between the thresholds
		# Sites are instantiated and added to the list of Site object
		if pval_threshold :
			if pval_threshold >= eval(j[8]):
				lth_pval_site_threshold = True
			else:
				lth_pval_site_threshold = False 
		
		if not pval_threshold:
			lth_pval_site_threshold = True
			
		if lth_score_threshold :
			if lth_score_threshold <= float(j[7]):
				lth_threshold = True
			else:
				lth_threshold = False 
		
		if not lth_score_threshold:
			lth_threshold = True
		
		if uth_score_threshold :
			if uth_score_threshold >= float(j[7]):
				uth_threshold = True
			else:
				uth_threshold = False 
		
		if not uth_score_threshold:
			uth_threshold = True
		
		if lth_threshold and uth_threshold and lth_pval_site_threshold:
			
			max_site_pvalue = max(max_site_pvalue,eval(j[8]))
			
			j[0] = Siteft(j[1],j[2],j[3],j[4],j[5],j[6],j[7],j[8],j[9],j[10])
			object_list.append(j[0])
	
	# Verbose message
	if verbose >=3 :
		time_warn("	Site object instantiated")
	
	
	filehandle.close()
	
	return object_list,list_limit,number_of_matrix,max_site_pvalue
	
	
def read_stdin(stdin,format,pval_threshold,lth_score_threshold, uth_score_threshold,verbose):
	""" 
	Read a list of sites enter on standard input
	Return respectively,
		- a list of site objects 
		- a list of sequence limits which are list with the sequence,start position and end position of the sequence
		- the number of matrix used to predict binding sites
	
	This function will take as arguments:
		stdin -> standard input name. Type = string
		format -> format of sites, it can be bed format "bed" or features format "ft". Type = string
		pval_threshold -> maximal p-value for read sites. Type = float
		lth_score_threshold -> minimal score for read sites. Type = float
		uth_score_threshold -> maximal score for read sites. Type = float
		verbose -> level of verbose. Type = int
	"""
	
	# Verbose message
	# Print information on standard error
	if verbose >= 3 :
		time_warn("Starting to parse")
	
	# Read sites in ft format
	if format == 'ft':
		
		# Extraction of text on standard input
		read = stdin
		read = read.strip()
		read = read.split("\n")
		tab_list = list(set(read)) # create a list without duplicates
		
		# Initiation of variables
		sites_list=[]
		i = 1
		list_limit=[]
		j=0
		max_site_pvalue = 0
		
		# treatment line by line and extraction of interesting information in a list
		for line in tab_list :
			
			# Extraction of tab information
			if line[0] != ";" and line[0]!="#":
				line = line.split("\t")
				
				# Extraction of limit sequence informations in a list (init_table) 
				# which is added to another list (list_limit)
				if line[1]== 'limit':
					init_table = []
					init_table.append(line[0])
					init_table.append(line[4])
					init_table.append(line[5])
					list_limit.append(init_table)
				
				# Extraction of sites informations	
				else :
					del line[1] # Remove the second element of the table ('site')
					line.insert(0,i) # Insertion of an index
					sites_list.append(line) # Add this list of table elements (site's informations) to the list of sites
					i = i+1
		
		# Extraction of other information like the number of matrix
		for line in read :
			
			# When there is the word "Matches" the followed informations are extracted
			if line[2:9]=="Matches":
				matrix=read[j+2:]
				matrix_list=[]
				
				# while the word "TOTAL" doesn't found the number find (second element of the row) are added in a list (matrix_list)
				for row in matrix :
					
					if row[3:8] == "TOTAL" :
						break
					else:
						num = row.split()
						matrix_list.append(int(num[1]))
			j= j+1
		
		# Number of matrix is the maximum of the matrix_list
		number_of_matrix = max(matrix_list)
		
		
		object_list=[]
		
		if verbose >=3 :
			time_warn("	Information extracted from the input")
			
		# every site is converted in siteft class object
		# If the site's p-value and score are between the thresholds
		# Sites are instantiated and added to the list of Site object
		for j in sites_list :
			
			# Replace incompatible element with a zero
			if j[7] == ".":
				j[7]=0
			
			# If the site's p-value and score are between the thresholds
			# Sites are instantiated and added to the list of Site object
			if pval_threshold :
				if pval_threshold >= eval(j[8]):
					lth_pval_site_threshold = True
				else:
					lth_pval_site_threshold = False 
			
			if not pval_threshold:
				lth_pval_site_threshold = True
				
			if lth_score_threshold :
				if lth_score_threshold <= float(j[7]):
					lth_threshold = True
				else:
					lth_threshold = False 
			
			if not lth_score_threshold:
				lth_threshold = True
			
			if uth_score_threshold :
				if uth_score_threshold >= float(j[7]):
					uth_threshold = True
				else:
					uth_threshold = False 
			
			if not uth_score_threshold:
				uth_threshold = True
			
			if lth_threshold and uth_threshold and lth_pval_site_threshold:
				
				max_site_pvalue = max(max_site_pvalue,eval(j[8]))
				j[0] = Siteft(j[1],j[2],j[3],j[4],j[5],j[6],j[7],j[8],j[9],j[10])
				
				object_list.append(j[0])
			
		# Verbose message
		if verbose >=3 :
			time_warn("	Site object instantiated")
		
		
		return object_list, list_limit, number_of_matrix, max_site_pvalue

	# Read sites in bed format
	if format == 'bed':
		
		# Extraction of text on standard input
		read= stdin
		read=read.strip()
		read=read.split("\n")
		
		
		# Initiation of variables
		sites_list=[]
		i = 1
		j = 0
		list_limit=[]
		max_site_pvalue = False
		
		# treatment line by line and extraction of interesting information in a list
		# identification of table with sites
		for line in read :
			j=j+1
			if line[0:13]=="browser dense":
				line_table=read[j-1:]
		
		# Extraction of table only
		del line_table[0:2]
		
		# Extraction of element to do the limits of the sequence
		for line in line_table:
			if line[0:7] == "browser":
				init_table = []
				# Extraction of sequence name
				first_init_table = line.split()
				first_init_table = first_init_table[-1].split(":")
				init_table.append(first_init_table[0])
				
				# Extraction of sequence positions
				first_init_table = first_init_table[-1].split("-")
				init_table.append(first_init_table[0])
				init_table.append(first_init_table[1])
				
				#add to the list with every limits of every sequence(list_limit)
				list_limit.append(init_table)
				
			else :
				# Treatment of sites table
				table = line.split()
				table.insert(0,i) # Insertion of an index
				sites_list.append(table) # Add this list of table elements (site's informations) to the list of sites
				i = i+1
		
		
		if verbose >=3 :
			time_warn("	Information extracted from the input")
		
		object_list=[]
		matrix = []
		
		# every site is converted in site class object
		for site in sites_list :
			
			# Replace incompatible element with a zero
			if site[5] == ".":
				site[5]=0
			
			# If the site's p-value and score are between the thresholds
			# Sites are instantiated and added to the list of Site object
			if lth_score_threshold :
				if lth_score_threshold<= float(site[5]):
					lth_threshold = True
				else:
					lth_threshold = False 
			
			if not lth_score_threshold:
				lth_threshold = True
			
			if uth_score_threshold :
				if uth_score_threshold >= float(site[5]):
					uth_threshold = True
				else:
					uth_threshold = False 
			
			if not uth_score_threshold:
				uth_threshold = True
			
			if lth_threshold and uth_threshold:
				
				# list of every transcription factors
				matrix.append(site[4])
				site[0] = Site(site[1],site[4],site[6],site[2],site[3],site[5])
				object_list.append(site[0])
		
		# We consider that there is only one matrix per sites (not True every time)
		matrix=list(set(matrix)) # Remove duplicates
		number_of_matrix = len(matrix) 
		
		if verbose >=3 :
			time_warn("	Site object instantiated")
		
		
		return object_list, list_limit , number_of_matrix, max_site_pvalue



def read_bed(file_name,lth_score_threshold,uth_score_threshold,verbose):
	"""
	Read bed file.
	Return respectively,
		- a list of site objects 
		- a list of sequence limits which are list with the sequence,start position and end position of the sequence
		- the number of matrix used to predict binding sites
	
	this function will take as arguments :
		file_name -> filehandle of read file in format bed. Type = string
		lth_score_threshold -> minimal score for read sites. Type = float
		uth_score_threshold -> maximal score for read sites. Type = float
		verbose -> level of verbose. Type = int
	"""
	
	if verbose >= 3 :
		time_warn("Starting to parse")
	
	# opening and reading of file
	filehandle = open(file_name,'r')
	read = filehandle.readlines()
	
	# Initiation of variables
	sites_list=[]
	i = 1
	j = 0
	list_limit=[]
	
	# treatment line by line and extraction of interesting information in a list
	# identification of table with sites
	for line in read :
		j=j+1
		if line[0:13]=="browser dense":
			line_table=read[j-1:]
	
	# Extraction of table only
	del line_table[0:2]
	
	# Extraction of element to do the limits of the sequence
	for line in line_table:
		if line[0:7] == "browser":
			init_table = []
			# Extraction of sequence name
			first_init_table = line.split()
			first_init_table = first_init_table[-1].split(":")
			init_table.append(first_init_table[0])
			
			# Extraction of sequence positions
			first_init_table = first_init_table[-1].split("-")
			init_table.append(first_init_table[0])
			init_table.append(first_init_table[1])
			
			#add to the list with every limits of every sequence(list_limit)
			list_limit.append(init_table)
			
		else :
			# Treatment of sites table
			table = line.split()
			table.insert(0,i) # Insertion of an index
			sites_list.append(table) # Add this list of table elements (site's informations) to the list of sites
			i = i+1
	
	
	if verbose >=3 :
		time_warn("	Information extracted from the input")
	
	object_list=[]
	matrix = []
	
	# every site is converted in site class object
	for site in sites_list :
		
		# Replace incompatible element with a zero
		if site[5] == ".":
			site[5]=0
		
		# If the site's p-value and score are between the thresholds
		# Sites are instantiated and added to the list of Site object
		if lth_score_threshold :
			if lth_score_threshold<= float(site[5]):
				lth_threshold = True
			else:
				lth_threshold = False 
		
		if not lth_score_threshold:
			lth_threshold = True
		
		if uth_score_threshold :
			if uth_score_threshold >= float(site[5]):
				uth_threshold = True
			else:
				uth_threshold = False 
		
		if not uth_score_threshold:
			uth_threshold = True
		
		if lth_threshold and uth_threshold:
			
			# list of every transcription factors
			matrix.append(site[4])
			site[0] = Site(site[1],site[4],site[6],site[2],site[3],site[5])
			object_list.append(site[0])
	
	# We consider that there is only one matrix per sites (not True every time)
	matrix=list(set(matrix)) # Remove duplicates
	number_of_matrix = len(matrix) 
	
	if verbose >=3 :
		time_warn("	Site object instantiated")
		
	filehandle.close()
	
	return object_list, list_limit, number_of_matrix



def write_crer(crer,filehandle):
	""" 
	This function write a file CRER by CRER in feature format
	There is no return
	
	it takes two arguments: 
		crer -> list of crer to write. Type = list of list with crer's informations. 
		filehandle -> standard output or the name of the output file. Type = string
	"""
	
	# Write in standard output or output file crer line, element correspond to the crer objects attributes
	filehandle.write("%s	%s	%s	%s	%d	%d	%d	%.2g	%.2e	%.2g	%.2g	%.2g	%d	%d\n" % (crer[0],"CRER",crer[1],crer[2],crer[3],crer[4],crer[5],crer[6],crer[7],crer[8],crer[9],crer[10],crer[11],crer[12]))
	
	# Close the output file if it has one
	if filehandle != sys.stdout :
		filehandle.close()


def time_warn(step):
	""" 
	This function write in standard err the time and the step associated to follow the script execution 
	Return the time and the step (in argument)
	
	it take one argument
		step -> name of the step. Type = string
	"""
	
	# Date extraction 
	T = str(datetime.datetime.now())
	
	# Write this information in standard err
	sys.stderr.write("%s	%s\n"% (T,step))



#===============================================================================
# Main of crer_scan
#===============================================================================

if __name__=='__main__':
	
	# To ignore the warnings
	sys.warnoptions = "ignore"
	
	# Start time
	time_start= time.strftime("%Y-%d-%m.%H:%M:%S")
	seconds_start= time.time()
	
	#===========================================================================
	# Arguments definition 
	# (thanks to argparse module)
	#===========================================================================
	
	parse = argparse.ArgumentParser()
	
	parse.add_argument('-i', action='store', dest='file',default=False, help = "If not specified, input is read from STDIN ")
	parse.add_argument('-in_format',action='store',dest='format',default='ft', help = "input_format. Default: ft (produced by RSAT matrix-scan and dna-pattern). Supported: ft, bed")
	parse.add_argument('-o',action = 'store', dest='outfile',help ="Output file in ft format")
	parse.add_argument('-v',action = 'store', dest = 'verbose',default=1,type = int ,help = "level of verbose. Messages are wrote on standard error. Supported: Integer = 1,2,3. By default : 1 = No message. Level 2 : moderately density of messages. Level 3 : High density")
	
	parse.add_argument('-s', action='store_true', dest='sort', default=None, help ="sort the list of sites. Very recommended. The sites are sorted by center position")
	parse.add_argument('-return_limits',action='store_true',dest='return_limits',default= None, help = "return every limits of sequences. By default : no return any limits")
	parse.add_argument('-return_limits_filtered',action='store_true',dest='return_limits_filtered',default= None, help = "return the limits filtered of the sequence. Only the sequence limits of CRERs. By default : no return any limits")
	parse.add_argument('-lth_score',action='store',dest='lth_site_score',default= None ,type= float, help = "minimal site score to be considered")
	parse.add_argument('-uth_score',action='store',dest='uth_site_score',default= None, type= float, help = "maximal site score to be considered")
	parse.add_argument('-uth_site_pval', action='store', dest='uth_site_pval', default= 1e-4, type = float, help='maximal p_value of sites to be considered. recommended to be the higher site p_value considered')
	
	# Arguments for return CRERs
	parse.add_argument('-lth_crer_size',action='store',dest='lth_crer_size',default=None,type= float, help = "minimal size of the enriched region (in bp). Default: minimal size = 30bp")
	parse.add_argument('-uth_crer_size',action='store',dest='uth_crer_size',default=None, type= float, help = "maximal size of the enriched region (in bp). Default:  maximal site = 500bp")
	
	parse.add_argument('-lth_crer_sites',action='store',dest='lth_crer_sites', default=None, type = float, help="minimal number of sites covered by the enriched region. Default: minimal number of sites = 2")
	parse.add_argument('-uth_crer_sites',action='store',dest='uth_crer_sites',default=None, type = float, help ="maximal number of sites covered by the enriched region. ")
	
	parse.add_argument('-lth_crer_sites_distance',action='store',dest='lth_crer_sites_distance',default=None,type = float,  help = "distance between successive sites to be considered. A minimal inter-site distance can be used to prevent overlap between redundant matrices. Default = minimal distance = 1")
	parse.add_argument('-uth_crer_sites_distance',action='store',dest='uth_crer_sites_distance',default=None,type = float,  help ="distance between successive sites to be considered. A maximal inter-site distance can be used to prevent merging distinct modules into a single one. Note: the maximal inter-site distance is one of the most influential parameters in cluster-buster. Default: maximal distance = 35")
	
	parse.add_argument('-uth_crer_pval',action='store',dest='uth_crer_pval',default = None, type= float, help = "maximal binomial p-value. Default: 1e-4")
	
	parse.add_argument('-uth_crer_eval',action='store',dest='uth_crer_eval',default = None, type= float, help = "maximal e-value. Default: 1e-4")
	
	parse.add_argument('-lth_crer_sig',action='store',dest='lth_crer_significance', type = float, default= None,help = "minimal binomial significance. Default: 2")

	parse.add_argument('-uth_overlap',action='store',dest='uth_overlap',default= None, type= int, help = "maximal overlap to define two distinct sites")

	parse.add_argument('-nopval',action='store_true',dest='nopval',default= False, help = "compute CRER without p value")
	parse.add_argument('-pre_table',action = 'store_true', dest = 'pre_table',default= False,help = "compute a table where is all possible p_value. Useful where there is a huge number of sites to scan.")
	
	# Arguments extraction
	args = parse.parse_args()
	
	# Variables for arguments extraction
	
	file,sort,format,outfile = args.file,args.sort,args.format,args.outfile
	verbose = args.verbose
	
	return_limits= args.return_limits
	return_limits_filtered = args.return_limits_filtered
	
	pval_threshold= args.uth_site_pval
	lth_score_threshold,uth_score_threshold = args.lth_site_score, args.uth_site_score
	
	lth_size, lth_sites, lth_dist, lth_sig = args.lth_crer_size, args.lth_crer_sites, args.lth_crer_sites_distance, args.lth_crer_significance
	uth_size, uth_sites, uth_dist, uth_pval,= args.uth_crer_size, args.uth_crer_sites, args.uth_crer_sites_distance, args.uth_crer_pval,
	uth_eval = args.uth_crer_eval
	uth_overlap = args.uth_overlap
	nopval = args.nopval
	pre_table=args.pre_table
	
	
	#===========================================================================
	# Parsing
	#===========================================================================
	
	# Do not accept '(,)' characters in output filename
	if outfile:
		expression_reg=re.compile("[\(\) ]")
		corresp=expression_reg.search(outfile)
		
		if corresp:
			   raise IOError ("Forbidden character in the outfile's name")
	
	# Verbose message
	# Print the start of the program in standard err
	if verbose >= 2:
		time_warn("Starting")
	
	# If there is not a file, take sites on standard input
	if file == False :
		read = sys.stdin.read()
		
		try:
			list_site, init_table,number_of_matrix,max_site_pvalue = read_stdin(read,format,pval_threshold, lth_score_threshold, uth_score_threshold,verbose)
			
		except (UnboundLocalError,IndexError,NameError,ValueError):
			sys.stderr.write("Parsing error :\n Sites in STDIN not in right format. /n Check the format: Default = ft. \n. Supported = ft or bed")
			sys.exit()
		
	else:
		
		# Test file format 
		if format=="ft": #test of the input file format 
			
			# Verbose message
			if verbose >= 3:
				time_warn("	FT format recognized")
			
			try :
				list_site,init_table,number_of_matrix,max_site_pvalue = read_features(file,pval_threshold, lth_score_threshold, uth_score_threshold,verbose)
			
			# Error messages
			except (UnboundLocalError,IndexError,NameError,ValueError):
				sys.stderr.write("Parsing error :\nFile input not in ft format")
				sys.exit()
				
			except ( IOError):
				sys.stderr.write("Parsing error :\nFile input not found ")
				sys.exit()
				
		elif format == "bed":
			
			# Verbose message
			if verbose >= 3:
				time_warn("	bed format recognized")
			
			try :
				list_site,init_table,number_of_matrix = read_bed(file,lth_score_threshold, uth_score_threshold,verbose)
				
			# Error messages
			except (UnboundLocalError,IndexError,NameError,ValueError):
				sys.stderr.write("Parsing error :\nFile input not in bed format")
				sys.exit()
			
			except ( IOError):
				sys.stderr.write("Parsing error :\nFile input not found ")
				sys.exit()
		
		else :
			
			raise IOError ("Unknowing format")
	
	# Verbose message
	if verbose >= 2:
		
		time_warn("Parsing Finished")
		
		
		
	
	#==========================================================================
	# Writing beginning of the file:
	#	- Command
	#	- Input file information
	#	- Description of the output file
	#==========================================================================
	
	# Verbose message
	if verbose >=3 :
		
		time_warn("Starting to write beginning of the file")
	
	# Definition of the output file name
	# By default the result is write on standard output
	if outfile :
		
		name = outfile
		
		try:
			filehandle = open(name,'w') # writing the beginning of the outfile

		except IOError:
			sys.stderr.write("Invalid path for Outfile: No such file or directory\n\n")
			filehandle = sys.stdout
			outfile = False
	
	else:
		filehandle = sys.stdout
	
	# writing the complete command with all arguments
	
	filehandle.write("; python3 crer_scan.py -i %s " % args.file)
	filehandle.write("-in_format %s " % args.format)
	if outfile:
		filehandle.write("-o %s " % args.outfile)
	filehandle.write("-v %d " % verbose)
	 
	if sort:
		filehandle.write("-s ")
	 
	if return_limits:
		filehandle.write("-return_limits ")
	 
	if return_limits_filtered:
		filehandle.write("-return_limits_filtered ")
	 
	if pval_threshold:
		filehandle.write("-uth_site_pval %g " % args.uth_site_pval)
	
	if lth_score_threshold:
		filehandle.write("-lth_score %g " % args.lth_site_score)
	 
	if uth_score_threshold :
		filehandle.write("-uth_score %g " % args.uth_site_score)
	
	if lth_size: 
		filehandle.write("-lth_crer_size %g " % args.lth_crer_size)
	
	if uth_size:
		filehandle.write("-uth_crer_size %g " % args.uth_crer_size)
	
	if lth_sites:
		filehandle.write("-lth_crer_sites %g " % args.lth_crer_sites)
	 
	if uth_sites:
		filehandle.write("-uth_crer_sites %g " % args.uth_crer_sites)
	
	if lth_dist:	 
		filehandle.write("-lth_crer_sites_distance %g " % args.lth_crer_sites_distance)
	
	if uth_dist: 
		filehandle.write("-uth_crer_sites_distance %g " % args.uth_crer_sites_distance)
	
	if uth_eval:
		filehandle.write("-uth_crer_eval %g " % args.uth_crer_eval)
	
	if uth_pval:
		filehandle.write("-uth_crer_pval %g " % args.uth_crer_pval)
	
	if lth_sig:
		filehandle.write("-lth_crer_sig %s " % args.lth_crer_significance)
	
	if uth_overlap:
		filehandle.write("-uth_overlap %g \n"% args.uth_overlap)
	
	if nopval :
		filehandle.write("-nopval " )
		
	if pre_table:
		filehandle.write("-pre_table ")
	

	# Information about the input file
	filehandle.write("; Input file\n;	input	%s\n" % file)
	filehandle.write("; File format	%s\n" % format )
	
	# Description of thresholds
	filehandle.write(";Thresholds	lower	upper\n")
	filehandle.write(";	crer_size	%s	%s\n" % (lth_size, uth_size))
	filehandle.write(";	crer_sites	%s	%s\n" % (lth_sites, uth_sites))
	filehandle.write(";	crer_sites_distance	%s	%s\n" % (lth_dist, uth_dist))
	filehandle.write(";	crer_pval	None	%s\n" % uth_pval)
	filehandle.write(";	crer_sig	%s\n	None" % lth_sig)
	filehandle.write(";	crer_eval	None	%s\n" % uth_eval)
	filehandle.write(";	site_pval	None	%s\n" % (pval_threshold))
	filehandle.write(";	site_score	%s	%s\n" % (lth_score_threshold, uth_score_threshold))
	if format == 'ft':
		filehandle.write(";	maximal site p-value	%g \n" % max_site_pvalue)
	
	# Description about the output table columns
	filehandle.write("; Output columns\n")
	filehandle.write(";	1	seq_id\n")
	filehandle.write(";	2	ft_type\n")
	filehandle.write(";	3	ft_name\n")
	filehandle.write(";	4	strand\n")
	filehandle.write(";	5	start = starting position of the crer, which corresponds to the start position of the leftmost site\n")
	filehandle.write(";	6	end = end position of the crer, which corresponds to the end position of the rightmost site\n")
	filehandle.write(";	7	hit_sum\n")
	filehandle.write(";	8	crer_sig\n")
	filehandle.write(";	9	crer_eval\n")
	filehandle.write(";	10	crer_pval\n")
	filehandle.write(";	11	hit_pval_product\n")
	filehandle.write(";	12	weight_sum\n")
	filehandle.write(";	13	crer_size\n")
	filehandle.write(";	14	crer_center_size = size taken in consideration for the p-value computation -> from center of the leftmost site to the center of the rightmost site (this size thus differs from the difference between crer_start and crer_end).\n")
	
	
	# Write the columns
	filehandle.write("#seq_id	ft_type	ft_name	strand	start	end	hit_sum	crer_sig	crer_eval	crer_pval	hit_pval_product	weight_sum	crer_size	crer_center_size\n")
	
	
	if outfile:
		filehandle.close()
	
	# Verbose message
	if verbose >=3 :
		time_warn("End writing beginning of the file")
	
	#===========================================================================
	#	 Computation of CRERs
	#	Creation and writing of CRERs
	#===========================================================================
	
	# Verbose message
	# Print the start of the crer computation in standard err
	if verbose >= 2:
		time_warn("Computation of crer")
	
	# Variables initiation
	sum_lengths = 0
	nb_crer_expect= 0
	nb_binom_computation = 0
	nb_sig_computation = 0
	nb_crer = 0
	
	# Do a temporary file 
	# Used to realized a pre-scan to compute e-value
	tmp = tempfile.TemporaryFile()
	
	# Create a SiteSet object with an empty list to put our site objects 
	SiteSet_object = SiteSet([], False)
	
	# Number of sites extracted
	length_list_site = len(list_site)
	
	# Computation of the binomial function
	# As the R function : pbinom(q,size,prob,lower.tail=FALSE)
	# p is the probability to have at least one site with the number of matrix.
	# With the maximum p-value of sites
	
	if format == 'ft':
		p = scipy.stats.binom.sf(0,number_of_matrix,max_site_pvalue)
	
	# With the p-value threshold
	if format == 'bed':
		p = scipy.stats.binom.sf(0,number_of_matrix,pval_threshold)
		
	# Computation of pre p-value in function of all possible size for crer
	if pre_table:
		if verbose >= 2:
			time_warn("Computation of table of all possible p_value")
			
		val_pval = {}
		for nb_expect_size in range(1,uth_size+1):
			N=nb_expect_size*2
			
			for nb_expect_sites in range(1,nb_expect_size+1):
				
				# Compute p-value for the number of positions with at least one site
				# pvalue positions = crer pvalue
				crer_pval = scipy.stats.binom.sf(nb_expect_sites-1,N,p)
				val_pval[nb_expect_size,nb_expect_sites]= crer_pval

	# Add site objects in SiteSet object
	for site in list_site:
		SiteSet_object.add_site(site)
	
	# Verbose message
	if verbose >=3 :
		time_warn("	Instanciated Site set on same seq and Site set classes")
	
	# Treat each SiteSetOnSameSeq object
	for SiteSetOnSameSeq in SiteSet_object.site_set_list():
		
		seq_id = SiteSetOnSameSeq.seq_id()
		
		
		
		#=======================================================================
		# Return limits of sequences
		#=======================================================================
		for seq in init_table :
			
			if seq[0]== seq_id:
				
				# Return all sequence limits
				if return_limits :
					# Write in temporary file
					if not nopval :
						data = "%s	limit	START_END	D	%s	%s	.	0	0	0	0	0	0\n" % (seq[0],seq[1], seq[2])
						data = data.encode('UTF-8')
						tmp.write(data)
					
					# Write in output file
					if nopval:
						if outfile:
							filehandle = open(name,'a') # opening of the file as alterable
						else :
							filehandle = sys.stdout
						
						filehandle.write("%s	limit	START_END	D	%s	%s	.	0	0	0	0	0	0\n" % (seq[0],seq[1], seq[2]))
				
				# Summary of all sequence's size
				sum_lengths = sum_lengths + (int(seq[2]) - int(seq[1]) + 1)
				
				
		#=======================================================================
		# Extract and sort sites in SiteSetOnSameSeq object by center positions
		#=======================================================================
		if sort :
			
			list_site_one_seq= SiteSetOnSameSeq.get_sorted_sites()
			
			if verbose >=3 :
				time_warn("	Sort sites")
				
		else:
			list_site_one_seq= SiteSetOnSameSeq.get_sites()
		
		#=======================================================================
		# Treat each site object (by its number i) in the site list 
		# as the first position of the crer for one sequence 
		#=======================================================================
		
		for i in range(0,len(list_site_one_seq)-1):
			
			# Define the beginning characteristics of the crer 
			
			# The number of sites is one (there is only the first site)
			# The score sum of the crer is the score of the first site
			# The first center position is the center position of the first site
			
			nb_sites= 1
			
			weight_crer = list_site_one_seq[i].weight()
			
			crer_start_center = list_site_one_seq[i].center()
			
			# Compute p-values product
			if format == 'ft':
				pval_prod = list_site_one_seq[i].Pval()
			
			# there is no site's p-value for bed files
			else:
				pval_prod = 0
			
			#=======================================================================
			# Treat for each ith site, each jth site.
			# which represent the following sites could be keep or not in the crer
			#=======================================================================
			
			for j in range(i+1,len(list_site_one_seq)):
				
				# Compute the start of the crer useful because sites are sorted by center
				crer_start = min(list_site_one_seq[i].start(),list_site_one_seq[j].start())
				
				# Compute overlap between two sites as a percent of overlap based on minimal site size. 
				site_size_1 = list_site_one_seq[j-1].end()-list_site_one_seq[j-1].start()+1
				site_size_2 = list_site_one_seq[j].end()-list_site_one_seq[j].start()+1
				minimal_site_size = min (site_size_1,site_size_2)
				overlap =  (min(list_site_one_seq[j].end(),list_site_one_seq[j-1].end()) - max(list_site_one_seq[j].start(),list_site_one_seq[j-1].start())) / minimal_site_size
				
				#=======================================================================
				# Overlap threshold
				#=======================================================================
				if uth_overlap:
					if overlap <= uth_overlap:
						overlap_threshold=True
					else:
						overlap_threshold=False
				if not uth_overlap:
					overlap_threshold=True
				
				if overlap_threshold : 
					# Verbose message
					if verbose >=3 :
							time_warn("	overlap is good")
					# Compute the size of the crer center to center 
					# with the max of the end position between the ith site end position and the jth site end position to avoid the inclusion of jth site in ith site
					size = max(list_site_one_seq[j].center(),list_site_one_seq[j-1].center()) - crer_start_center +1
					
					#=======================================================================
					# Size threshold
					# Verification of the size of the crer
					#=======================================================================
					if lth_size:
						if lth_size <= size:
							lth_size_threshold = True
						else :
							lth_size_threshold = False
								
					if not lth_size:
						lth_size_threshold = True
					   
					if uth_size:
						if uth_size >= size:
							uth_size_threshold = True
						else:
							uth_size_threshold = False
					if not uth_size:
						uth_size_threshold = True
					
					
					if lth_size_threshold and uth_size_threshold :
						
						# Verbose message
						if verbose >=3 :
							time_warn("	size is good")
						
						# Compute the distance between following sites in the crer
						# Define by the start position of the jth site and the previous site, jth-1 site.
						# ( We could try to estimate the overlapping)
						dist = list_site_one_seq[j].center() - list_site_one_seq[j-1].center() 
						
						#=======================================================================
						# Distance threshold
						# Verification of the distance between following sites
						# If the distance is not correct the jth site is not keeping in the crer 
						#=======================================================================
						if lth_dist:
							if lth_dist <= dist:
								lth_dist_threshold = True
							else :
								lth_dist_threshold = False
								
						if not lth_dist:
							lth_dist_threshold = True
						
						if uth_dist:
							if uth_dist >= dist:
								uth_dist_threshold = True
							else:
								uth_dist_threshold = False
								
						if not uth_dist:
							uth_dist_threshold = True
						
						if lth_dist_threshold and uth_dist_threshold : 
							
							# Verbose message
							if verbose >=3 :
								time_warn("	distance is good")
								
							# Increment the number of sites
							nb_sites = nb_sites + 1 
							
							# Compute of the weight sum with addition of the jth site score
							weight_crer = weight_crer + list_site_one_seq[j].weight()
							
							# Compute pval_product as multiplication of sites p-values in the crer
							# Only possible for features files because in bed files we don't have site p-values.
							if format == 'ft':
								pval_prod = pval_prod * list_site_one_seq[j].Pval()
							
							#=======================================================================
							# Sites number threshold
							# Verification of the number of sites 
							# Default minimal value is 2
							# if it is not correct we don't keep the jth site in the crer
							#=======================================================================
							if lth_sites:
								if lth_sites <= nb_sites:
									lth_sites_threshold = True
								else :
									lth_sites_threshold = False
								
							if not lth_sites:
								lth_sites_threshold = True
						
							if uth_sites:
								if uth_sites >= nb_sites:
									uth_sites_threshold = True
								else:
									uth_sites_threshold = False
								
							if not uth_sites:
								uth_sites_threshold = True
							
							if lth_sites_threshold and uth_sites_threshold :
								
								# Verbose message
								if verbose >=3 :
										time_warn("	Number of sites is good")
								
								# Compute of the density of the crer
								density_crer = nb_sites/size
								
								# Compute the expected number of sites, random expectation of sites
								exp_sites = pval_prod * size
								
								#=======================================================================
								# Expected sites number threshold
								# keep crers with more sites than random expectation
								#=======================================================================
								
								if nb_sites >= exp_sites:
									
									#=======================================================================
									# With p-value computation
									#=======================================================================
									if not nopval :
									
										# Computation for two strands
										N = size*2
										
										try:
											# Computation is already done in the pre-table
											if pre_table:
												crer_pval = val_pval[size,nb_sites]
											
											else:
												# Compute p-value for the number of positions with the number of sites 
												# the probability is probability to observe a hit for at least one matrix at a given position
												# pvalue for all positions = crer p-value
												crer_pval = scipy.stats.binom.sf(nb_sites-1,N,p)
												nb_binom_computation += 1
												
												# Verbose message
												if verbose >=3 :
													time_warn("		%d binomial computation" % nb_binom_computation)
											
										except (RuntimeWarning,ValueError):
											# Verbose message
											if verbose >= 2 :
												sys.stderr.write("Error in binomial computation\n")
											crer_sig = 0
											crer_pval = 0
										
										#=======================================================================
										# P-value threshold
										# E-value threshold
										#=======================================================================
										if uth_pval:
											if uth_pval >= crer_pval:
												uth_pval_threshold = True
											else :
												uth_pval_threshold = False
								
										if not uth_pval:
											uth_pval_threshold = True
						
										if uth_eval:
											if uth_eval >= crer_pval:
												uth_eval_threshold = True
											else:
												uth_eval_threshold = False
								
										if not uth_eval:
											uth_eval_threshold = True
											
										if uth_pval_threshold and uth_eval_threshold:
											
											crer_end = max(list_site_one_seq[j].end(),list_site_one_seq[i].end())
											crer_end_center = max(list_site_one_seq[j].center(),list_site_one_seq[i].center())
											
											crer_size = crer_end-crer_start+1
										
											nb_crer_expect += 1
											data = "%s	%s	%s	%d	%d	%d	%.2g	%.2g	%.2g	%d	%d	%d\n" % (list_site_one_seq[i].seq_id(),"crer","DR", crer_start,crer_end,nb_sites,0,crer_pval,pval_prod,weight_crer,crer_size,size)
										
											# Writing of this crer in the temporary file
											data = data.encode('UTF-8')	
											tmp.write(data)
											
											# Verbose message
											if verbose >=3 :
												time_warn("	Printing %d crer in tmp file" % nb_crer_expect)
									
									# Without p-value computation
									if nopval :
										
										# initiate variables not defined
										crer_sig = 0
										crer_pval = 0
										nb_binom_computation = 0
										e_val = 0
										
										nb_crer += 1
										crer_end = max(list_site_one_seq[j].end(),list_site_one_seq[i].end())
										crer_end_center = max(list_site_one_seq[j].center(),list_site_one_seq[i].center())
										crer_size = crer_end-crer_start+1
										
										# Put crer description in a list
										crer = [list_site_one_seq[i].seq_id(),"crer","DR", crer_start,crer_end,nb_sites,crer_sig,e_val,crer_pval,pval_prod,weight_crer,crer_size,size]
										
										# Open right file
										if outfile:
											#opening of the file as alterable
											filehandle = open(name,'a')
										else :
											filehandle = sys.stdout
										
										# Return the limit of sequence where there is the crer
										if return_limits_filtered:
											for seq in init_table :
												if crer[0]== seq[0]:
													filehandle.write("%s	limit	START_END	D	%s	%s	.	0	0	0	0	0	0	0\n" % (seq[0],seq[1], seq[2]))
													init_table.remove(seq)
										
										# Write the crer on file (or standard output)
										write_crer(crer,filehandle)
										
										# Print a verbose message on standard error every five crer 
										if verbose == 2 and nb_crer%5 == 0:
											time_warn("	Printing %d crer" % nb_crer)
										
										# Print a verbose message on standard error for every crer 
										if verbose >=3 :
											time_warn("	Printing %d crer" % nb_crer)
	
	# Read the temporary file
	if not nopval:
		
		tmp.seek(0)
		read = tmp.read()
		read = read.decode()
		read = read.strip().split('\n')
		nb_crer = 0
		
		if not nopval :
			for ligne in read :
				crer_tmp = ligne.split('\t')
				# Extract limits of sequences, only here if there is retur_limits option
				# Write limits on file
				if crer_tmp[1] == 'limit':
					if outfile:
						filehandle = open(name,'a')#opening of the file as alterable
					else :
						filehandle = sys.stdout
					
					filehandle.write("%s	limit	START_END	D	%s	%s		0	0	0	0	0	0	0\n" % (crer_tmp[0],crer_tmp[4], crer_tmp[5]))
			
				# Compute e-value and significance
				# E-value depend on number of expected crer 
				else:
					e_val = float(crer_tmp[7])* nb_crer_expect
					
					# For this moment we define sig_base as one and after we can compute this.
					# Computation of crer significance
					
					crer_sig = -math.log10(e_val)
					
					# Verbose message for significance
					if verbose >=3 :
						nb_sig_computation += 1
						time_warn("		%d significance computation" % nb_sig_computation)
					
					#=======================================================================
					# Significance threshold
					#=======================================================================
					if lth_sig:
						if lth_sig <= crer_sig:
							lth_sig_threshold = True
						else :
							lth_sig_threshold = False
								
					if not lth_sig:
						lth_sig_threshold = True
					
					if lth_sig_threshold :
						
						crer = [crer_tmp[0],crer_tmp[1],crer_tmp[2],int(crer_tmp[3]),int(crer_tmp[4]),int(crer_tmp[5]),crer_sig,e_val,float(crer_tmp[7]),float(crer_tmp[8]),int(crer_tmp[9]),int(crer_tmp[10]),int(crer_tmp[11])]
						
						nb_crer += 1
					
						# Write the crer
						
						if outfile:
							#opening of the file as alterable
							filehandle = open(name,'a')
						else :
							filehandle = sys.stdout
						
						# Write the limit of sequence where there is the crer
						if return_limits_filtered:
							for seq in init_table :
								if seq[0]== crer_tmp[0]:
									filehandle.write("%s	limit	START_END	D	%s	%s	.	0	0	0	0	0	0	0\n" % (seq[0],seq[1], seq[2]))
									init_table.remove(seq)
						
						# Write the crer on file (or standard output)
						write_crer(crer,filehandle)
					
						# Print a verbose message on standard error every five crer 
						if verbose == 2 and nb_crer%5 == 0:
							time_warn("	Printing %d crer" % nb_crer)
						# Print a verbose message on standard error for every crer 
						if verbose >=3 :
							time_warn("	Printing %d crer" % nb_crer)
		
		# Write an error if there is no crer in temporary file
	   # except (IndexError):
		#	sys.stderr.write("No Crer found\n")
			
	if verbose >=3 :
		time_warn("	Computation of crer finished")
	
	#=======================================================================
	# Write the end of output file
	#=======================================================================
	
	if outfile:
		filehandle = open(name,'a')#opening of the file as alterable
	else :
		filehandle = sys.stdout
	
	
	if verbose >=3 :
		time_warn("Writing the end of output file")
	
	# Statistics
	filehandle.write("; Number of binomial	%d\n"% nb_binom_computation)
	filehandle.write("; Number of sites scanned	%d\n" % length_list_site)
	filehandle.write("; Sum of sequence lengths	%d\n" % sum_lengths)
	filehandle.write("; Number of Crer	%d\n" % nb_crer)
	
	# finished times
	time_end= time.strftime("%Y-%d-%m.%H:%M:%S")
	
	# Write start time and finished time
	filehandle.write("; Job started	%s\n" % time_start)
	filehandle.write("; Job done	%s\n" % time_end)
	
	# finished time in second
	seconds_end=time.time()
	
	# Write time of the scan
	filehandle.write("; Seconds	%f\n" % (seconds_end - seconds_start))

	if outfile:
		filehandle.close()
	
	# Print time when the script is finished on standard error
	if verbose >= 2:
		time_warn("Finished")

