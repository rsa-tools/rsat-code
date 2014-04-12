#!/usr/bin/python3
# -*- coding: utf-8 -*-

import argparse
import pydoc
import sys
import crer_scan
import time
import datetime
import scipy.stats
import math
import tempfile
import re

# Creation of site's class
class Site(object):
	"""This class is intended for sites of transcription factors identified by the matrix-scan program of the Regulatory Sequence Analysis Tools (RSAT). 
	The attributes of this class are : seq_id -> the identity of the sequence
										ft_name -> the name/identity of the transcription factor
										strand -> the strand of DNA where is the site of transcription factor
										start -> the beginning of the site given by the number of the first nucleotide
										end -> the end of the given site by the number of the last nucleotide
										weight -> the weight or score of the given site
										"""
										
	
	def __init__(self,seq_id,ft_name,strand,start,end,weight):
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
	
	
#Creation of subclass to treat features file		
class Siteft(Site):
	"""This subclass is intended for sites of transcription factors identified by the matrix scan software from feature's files
	because depending on the format of the input file we don't have the same information
	the bed format contains less information than the format feature (. ft.)
	The attributes of this class are those of the class site and  the following attributes : sequence -> the sequence of the site of the transcription factor 
																							 Pval -> the P-value of this site identified by the matrix scan software
																							 ln_Pval -> the natural logarithm of the p-value 
																							 sig -> the binomial significance of the site"""
	
	def __init__(self,seq_id,ft_name,strand,start,end,sequence,weight,Pval,ln_Pval,sig):
		Site.__init__(self,seq_id,ft_name,strand,start,end,weight)
		self.__sequence = sequence
		self.__Pval = float(Pval)
		self.__ln_Pval = float(ln_Pval)
		self.__sig = float(sig)
		
	def sequence(self):
		"""Return site's sequence and do not take arguments"""
		return self.__sequence
	
	def Pval(self):
		"""Return site's P-value v"""
		return self.__Pval
	
	def ln_Pval(self):
		"""Return logarithm site's P-value and do not take arguments"""
		return self.__ln_Pval
	
	def sig(self):
		"""Return site's significance and do not take arguments"""
		return self.__sig
	
	
class Siteset_on_same_seq(object):
	"""this class relates to a site collection on the same sequence as you can analyze several sites from different sequences at the same time
	The attributes of this class are : seq_id -> the identity of the sequence
									list_sites -> the list of sites on the same sequence
									sorted -> this argument is false if the list of sites on a same sequence is not sorted and if it is this argument is true"""
	
	def __init__(self,seq_id, list_sites, sorted):
		self.__seq_id = seq_id
		self.__list_sites = list_sites
		self.__sorted = sorted
	
	def seq_id(self):
		"""Return the identity of the sequence"""
		return self.__seq_id
	
	def get_sites(self):
		"""Return a list of sites who are in the same sequence"""
		return self.__list_sites
	
	def sorted(self):
		"""Return if the list is sorted or not"""
		return self.__sorted
	
	def sort(self):
		"""This method sorts the sites according to their attributes 'start' thanks to a local function which is designated 'lambda'
		once the list sorted, 'sorted' attribute is changed"""
		self.__list_sites = sorted(self.__list_sites,key=lambda x: x.center()) #this function sorts by start's position
		self.__sorted = True
		
	def add_site_same_seq(self,site):
		"""This method adds a site to the list of site located on the same sequence. For this, we verify that the site is well to add the sequence concerned"""
		if self.seq_id() == site.seq_id():
			
			self.__list_sites.append(site)
			
	def get_sorted_sites(self):
		"""this method allows us to obtain the list of sites sorts"""
		self.sort()
		return self.__list_sites
	
	
class Site_set(object):
	"""this class relates to a collection of several collections of sites on the same sequence.
	Example, if we consider three sequences at the same time we have a collection containing three collections
	 that corresponds to each of the sequences.
	 The attributes of this class are : liste_site_set -> This is a list containing the different lists where each corresponds to a sequence
	 									sorted -> this argument is false if the list of sites on a same sequence is not sorted and if it is this argument is true"""
	
	def __init__(self,liste_site_set, sorted):
		self.__liste_site_set = liste_site_set
		self.__sorted = sorted
		
	def liste_site_set(self):
		"""Return  the list of site set"""
		return self.__liste_site_set
	
	def sorted(self):
		"""Return if the collection is sorted or not"""
		return self.__sorted
	
	def sort(self):
		"""This method sorts the sites by the sort method of Siteset_on_same_seq class """
		for i in self.__liste_site_set:
			i.sort()
		self.__sorted = True
		
	def add_site(self,site):
		"""This method allows us to add a site in the collection of site collections.
		This is done in several stages.
		At first we check if the collection is empty collections whether we add the site object by instantiating the class Siteset_on_same_seq.
		If the collection of collections is not empty then we check if the site is added to a sequence investigate 
		whether it is added to the good list if not we create a new object of the class Siteset_on_same_seq for him. """
		  
		if self.__liste_site_set == []:
			self.__liste_site_set.append(Siteset_on_same_seq(site.seq_id(),[site],False))
			
		else:
			
			l = ''
			for y in self.liste_site_set():
			
				if y.seq_id() == site.seq_id():
					y.add_site_same_seq(site)
					l = 'ok'
			
			if l != 'ok': 
				w=Siteset_on_same_seq(site.seq_id(),[site],False)
				self.__liste_site_set.append(w)
	
	def get_sites(self):
		"""Returns a list of all the sites studied"""
		get_sites=[]
		for siteset_on_same_seq in self.liste_site_set():
			list_siteset=siteset_on_same_seq.get_sites()
			for sites in list_siteset:
				get_sites.append(sites)
		return get_sites
	
	def get_sorted_sites(self):
		"""Returns a sorted list of all the sites studied"""
		list_sites=self.get_sites()
		get_sorted_sites=sorted(list_sites,key=lambda x: x.center())
		return get_sorted_sites
	
	
#Creation of CRER class
class CRER(object):
	"""This class relates to a CRERs (cis regulation enriched Region  )
	The attributes of this class are : 
						seq_id -> the identity of the sequence
						ft_name -> the name/identity of the transcription factor
						strand -> the strand of DNA where is the site of transcription factor
						start -> the beginning of the site given by the number of the first nucleotide
						end -> the end of the given site by the number of the last nucleotide
						hit_sum -> 
						crer_sig -> the significance of the CRER
						crer_eval -> the e-value of the CRER
						crer_pval -> the P-value of the CRER
						hit_pval_product -> hit_pval_product
						weight_sum -> the sum of weight
						crer_size -> the CRER'size"""
	
	def __init__(self,seq_id,ft_name,strand,start,end,hit_sum,crer_sig,crer_eval,crer_pval,hit_pval_product, weight_sum,crer_size,crer_center_size):
		self.__seq_id = seq_id
		self.__ft_name = ft_name
		self.__strand = strand
		self.__start = int(start)
		self.__end = int(end)
		self.__hit_sum = int(hit_sum)
		self.__crer_sig = float(crer_sig)
		self.__crer_eval = float(crer_eval)
		self.__crer_pval= float(crer_pval)
		self.__hit_pval_product = float(hit_pval_product)
		self.__weight_sum = float(weight_sum)
		self.__crer_size = int(crer_size)
		self.__crer_center_size = int(crer_center_size)
		
	def seq_id(self):
		"""Return crer's identity"""
		return self.__seq_id
	def ft_name(self):
		"""Return transcription factor's name"""
		return self.__ft_name
	def strand(self):
		"""Return strand's name"""
		return self.__strand
	def start(self):
		"""Return start position"""
		return self.__start
	def end(self):
		"""Return end position"""
		return self.__end
	def hit_sum(self):
		"""Return sum of hits"""
		return self.__hit_sum
	def crer_sig(self):
		"""Return crer's significance"""
		return self.__crer_sig
	def crer_eval(self):
		"""Return crer's E-value"""
		return self.__crer_eval
	def crer_pval(self):
		"""Return crer's P-value"""
		return self.__crer_pval
	def hit_pval_product(self):
		"""Return p-value of hits producted"""
		return self.__hit_pval_product
	def weight_sum(self):
		"""Return weight's sum"""
		return self.__weight_sum
	def crer_size(self):
		"""Return crer's size"""
		return self.__crer_size
	def crer_center_size(self):
		"""Return crer's size center to center"""
		return self.__crer_center_size

def read_features(file_name,pval_threshold,lth_score_threshold, uth_score_threshold, verbose):
	""" Read a matrix scan's output file.
	Return a list of site objects.
	This function will take as arguments:
	
		file_name = filehandle of read file in feature format
		pval_threshold = maximal p-value for read sites
		lth_score_threshold = minimal score for read sites
		uth_score_threshold = maximal score for read sites
	"""
	if verbose >= 3 :
		TimeWarn("Starting to parse")
	
	# opening and reading of file
	lecture=open(file_name,'r')
	lect=lecture.readlines()
	l_tab=list(set(lect))
	liste=[]
	i = 1
	list_limit=[]
	j = 0
	# treatement line by line and extraction of interesting information in a list
	for ligne in l_tab :
		if ligne[0] != ";" and ligne[0]!="#":
			lis = ligne.split()
			
			if lis[1]== 'limit':
				init_table1=lis
				init_table = []
				init_table.append(init_table1[0])
				init_table.append(init_table1[4])
				init_table.append(init_table1[5])
				list_limit.append(init_table)
				
			else:
				del lis[1]
				
				lis.insert(0,i)
				liste.append(lis)
				i = i+1
				
	for ligne in lect :
		
		if ligne[2:9]=="Matches":
			matrice=lect[j+2:]
			
			matrix=[]
			
			for row in matrice:
				
				if row[3:8] == "TOTAL" :
					break
				else:
					num = row.split()
					matrix.append(int(num[1]))
					
		j= j+1
		
	number_of_matrix = max(matrix)

	liste_object=[]
	
	if verbose >=3 :
		TimeWarn("	Information extracted from the input")
	
	
	# every site is converted in siteft class object
	for j in liste :
		if j[7] == ".":
			j[7]=0
		if eval(j[8]) <= pval_threshold and lth_score_threshold<= float(j[7]) <= uth_score_threshold :
			j[0] = Siteft(j[1],j[2],j[3],j[4],j[5],j[6],j[7],j[8],j[9],j[10])
			liste_object.append(j[0])
	
	if verbose >=3 :
		TimeWarn("	Site object instantiated")
	
	
	lecture.close()
	return liste_object,list_limit,number_of_matrix
	
	
def read_stdin(stdin,pval_threshold,lth_score_threshold, uth_score_threshold,verbose):
	
	# Print information on standard error
	if verbose >= 3 :
		TimeWarn("Starting to parse")
	
	l= stdin
	l=l.strip()
	l=l.split("\n")
	l_tab = list(set(l))
	liste=[]
	i = 1
	list_limit=[]
	j=0
	# treatement line by line and extraction of interesting information in a list
	for ligne in l_tab :
		if ligne[0] != ";" and ligne[0]!="#":
			ligne = ligne.split("\t")
			if ligne[1]== 'limit':
				init_table1=ligne
				init_table = []
				init_table.append(init_table1[0])
				init_table.append(init_table1[4])
				init_table.append(init_table1[5])
				list_limit.append(init_table)
			else :
				del ligne[1]
				ligne.insert(0,i)
				liste.append(ligne)
				i = i+1
	
	for ligne in l :
		if ligne[2:9]=="Matches":
			matrice=l[j+2:]
			#print(matrice)
			matrix=[]
			
			for row in matrice:
				
				if row[3:8] == "TOTAL" :
					break
				else:
					num = row.split()
					#print(num[1])
					matrix.append(int(num[1]))
					#print(matrix)
		j= j+1
	
	number_of_matrix = max(matrix)
	
	del liste[0]
	liste_object=[]
	
	if verbose >=3 :
		TimeWarn("	Information extracted from the input")
		
	# every site is converted in siteft class object
	for j in liste :
		if j[7] == ".":
			j[7]=0
		if eval(j[8]) <= pval_threshold and lth_score_threshold<= float(j[7]) <= uth_score_threshold :
			j[0] = Siteft(j[1],j[2],j[3],j[4],j[5],j[6],j[7],j[8],j[9],j[10])
			liste_object.append(j[0])
			
	if verbose >=3 :
		TimeWarn("	Site object instantiated")
	
	return liste_object, list_limit, number_of_matrix

def read_bed(file_name,lth_score_threshold,uth_score_threshold,verbose):
	"""Read bed file.
	Return a list of site 
	
	this function will take as arguments :
		file_name = filehandle of read file in format bed
		lth_score_threshold = minimal score for read sites
		uth_score_threshold = maximal score for read sites
	"""
	
	if verbose >= 3 :
		TimeWarn("Starting to parse")
	
	# opening and reading of file
	lecture=open(file_name,'r')
	l=lecture.readlines()
	liste=[]
	i = 1
	e=0
	init_table=[]
	
	# treatment line by line and extraction of interesting information in a list
	for ligne in l :
		e=e+1
		if ligne[0:13]=="browser dense":
			ligne1=l[e-1:]
	
	
	del ligne1[0:2]
	
	for u in ligne1:
		if u[0:7] == "browser":
			init_table1=u
			init_table2 = []
			init_table1= init_table1.split()
			init_table1=init_table1[-1].split(":")
			init_table2.append(init_table1[0])
			init_table1=init_table1[-1].split("-")
			init_table2.append(init_table1[0])
			init_table2.append(init_table1[1])
			init_table.append(init_table2)
			
		else :
			tab = u.split()
			tab.insert(0,i)
			liste.append(tab)
			i = i+1
			
	liste_object=[]
	
	if verbose >=3 :
		TimeWarn("	Information extracted from the input")
	
	
	# every site is converted in site class object
	matrix = []
	for j in liste :
		if j[5] == ".":
			j[5]=0
		if lth_score_threshold<= float(j[5]) <= uth_score_threshold:
			matrix.append(j[4])
			j[0] = Site(j[1],j[4],j[6],j[2],j[3],j[5])
			liste_object.append(j[0])
	
	matrix=list(set(matrix))
	number_of_matrix = len(matrix)
	
	if verbose >=3 :
		TimeWarn("	Site object instantiated")
		
	lecture.close()
	
	return liste_object, init_table, number_of_matrix
	
	
	
def write_crer(objet,filehandle):
	""" This function write  a file CRER by CRER in feature format
	
	it takes two arguments: 
		object = CRER object
		filehandle = the name of the outfile
	"""
	
	
	filehandle.write("%s	%s	%s	%s	%d	%d	%d	%.2g	%.2e	%.2g	%.2g	%.2g	%d	%d\n" % (objet.seq_id(),"CRER",objet.ft_name(),objet.strand(), objet.start(),objet.end(),objet.hit_sum(),objet.crer_sig(),objet.crer_eval(),objet.crer_pval(),objet.hit_pval_product(),objet.weight_sum(),objet.crer_size(),objet.crer_center_size()))
	if filehandle != sys.stdout :
		filehandle.close()

def TimeWarn(step):
	""" This function write in standard err the time and the step associated to follow the script execution 
	
	it take one argument
		step = name of the step. type = string
	""" 
	# Date extraction 
	T = str(datetime.datetime.now())
	# Write this information in standard err
	sys.stderr.write("%s	%s\n"% (T,step))
# For generate documentation
#pydoc.help(crer_scan)


#===============================================================================
# Main of crer_scan
#===============================================================================
if __name__=='__main__':
	
	sys.warnoptions = "ignore"
	
	time_start= time.strftime("%Y-%d-%m.%H:%M:%S")
	seconds_start= time.time()
	
	#===========================================================================
	# Arguments definition 
	# (thanks to argparse module)
	#===========================================================================
	
	parse = argparse.ArgumentParser()
	
	parse.add_argument('-i', action='store', dest='file',default=False,help = """If not specified, input is read from STDIN """)
	parse.add_argument('-s', action='store_true', dest='sort', default=False, help ="sort list of site recommanded if input file is not producted by RSAT matrix-scan")
	parse.add_argument('-in_format',action='store',dest='format',default='ft',help= "input_format. Default: ft (produced by RSAT matrix-scan and dna-pattern). Supported: ft, bed")
	
	parse.add_argument('-lth_crer_size',action='store',dest='lth_crer_size',default=30,type= float, help = "minimal size of the enriched region (in bp). Default: min size = 30bp")
	parse.add_argument('-uth_crer_size',action='store',dest='uth_crer_size',default=500, type= float, help = "maximal size of the enriched region (in bp). Default:  max site = 500bp")
	
	parse.add_argument('-lth_crer_sites',action='store',dest='lth_crer_sites', default=2, type = float, help="minimal number of sites covered by the enriched region. Default: min number of sites = 2")
	parse.add_argument('-uth_crer_sites',action='store',dest='uth_crer_sites',default= 1000, type = float, help ="maximal number of sites covered by the enriched region")
	
	parse.add_argument('-lth_crer_sites_distance',action='store',dest='lth_crer_sites_distance',default=1,type = float,  help = "distance between successive sites to be considered. A minimal inter-site distance can be used to prevent overlap between redundant matrices. Default = min distance = 1")
	parse.add_argument('-uth_crer_sites_distance',action='store',dest='uth_crer_sites_distance',default=1000,type = float,  help ="distance between successive sites to be considered. A maximal inter-site distance can be used to prevent merging distinct modules into a single one. Note: the maximal inter-site distance is one of the most influential parameters in cluster-buster (default 35bp)")
	
	parse.add_argument('-uth_crer_pval',action='store',dest='uth_crer_pval',default = 1e-4, type= float, help = "maximal binomial p-value")
	
	parse.add_argument('-uth_crer_eval',action='store',dest='uth_crer_eval',default = 1e-4, type= float, help = "maximal e-value")
	
	parse.add_argument('-lth_crer_sig',action='store',dest='lth_crer_significance', type = float, default= 2,help = "minimal binomial significance")
	
	parse.add_argument('-uth_site_pval',action='store',dest='uth_site_pval',default= 1e-4, type=float, help = "maximal p-value of sites to be considered" )
	
	parse.add_argument('-lth_score',action='store',dest='lth_site_score',default= 0,type= float, help = "minimal site score to be considered")
	parse.add_argument('-uth_score',action='store',dest='uth_site_score',default= 1000, type= float, help = "maximal site score to be considered")
	
	parse.add_argument('-uth_overlap',action='store',dest='uth_overlap',default= 1, type= int, help = "maximal overlap to define two distinct sites")
	
	parse.add_argument('-return_limits',action='store_true',dest='return_limits',default= False, help = "return the limit of the sequence")
	
	parse.add_argument('-return_limits_filtered',action='store_true',dest='return_limits_filtered',default= False, help = "return the limit filtered of the sequence")
	
	parse.add_argument('-nopval',action='store_true',dest='nopval',default= False, help = "compute crer without p value")
	
	parse.add_argument('-o',action = 'store', dest='outfile',help ="Output file in ft format")
	
	parse.add_argument('-v',action = 'store', dest = 'verbosity',default=1,type = int ,help = "Verbosity")
	
	parse.add_argument('-pre_table',action = 'store_true', dest = 'pre_table',default= False,help = "compute a table where is all possible p_value ")
	
	# Arguments extraction
	args = parse.parse_args()
	
	# Variables for arguments extraction
	file,sort,format,outfile= args.file,args.sort,args.format,args.outfile
	lth_size, lth_sites, lth_dist, lth_sig = args.lth_crer_size, args.lth_crer_sites, args.lth_crer_sites_distance, args.lth_crer_significance
	uth_size, uth_sites, uth_dist, uth_pval,= args.uth_crer_size, args.uth_crer_sites, args.uth_crer_sites_distance, args.uth_crer_pval,
	pval_threshold= args.uth_site_pval
	lth_score_threshold,uth_score_threshold = args.lth_site_score, args.uth_site_score
	verbose = args.verbosity
	uth_overlap = args.uth_overlap
	uth_eval = args.uth_crer_eval
	return_limits= args.return_limits
	return_limits_filtered = args.return_limits_filtered
	nopval = args.nopval
	pre_table=args.pre_table
	
	if outfile:
		expression_reg=re.compile("[\(\) ]")
		corresp=expression_reg.search(outfile)
		
		if corresp:
			raise IOError ("Forbidden character in the outfile's name")
	# Print the start of the program in standard err
	if verbose >= 2:
		TimeWarn("Starting")
	
	#===========================================================================
	# Reading file
	#===========================================================================
	
	# Test if a value is return in argument : file
	step1_start=time.time()
	if file == False :
		lecture = sys.stdin.read()
		list_Site, init_table,number_of_matrix = read_stdin(lecture,pval_threshold, lth_score_threshold, uth_score_threshold,verbose)
		
	else:
		# Test file format 
		if format=="ft": #test of the input file format 
			
			if verbose >= 3:
				TimeWarn("	FT format recognized")
			
			try :
				list_Site,init_table,number_of_matrix = read_features(file,pval_threshold, lth_score_threshold, uth_score_threshold,verbose)
				
			except (UnboundLocalError,IndexError,NameError,ValueError,):
				sys.stderr.write("Parsing error :\nFile input not in ft format")
				sys.exit()
			except ( IOError):
				sys.stderr.write("Parsing error :\nFile input not found ")
				sys.exit()
				
		elif format=="bed":
			
			if verbose >= 3:
				TimeWarn("	bed format recognized")
			
			try :
				list_Site,init_table,number_of_matrix = read_bed(file,lth_score_threshold, uth_score_threshold,verbose)
				
			except (UnboundLocalError,IndexError,NameError,ValueError):
				sys.stderr.write("Parsing error :\nFile input not in bed format")
				sys.exit()
			except ( IOError):
				sys.stderr.write("Parsing error :\nFile input not found ")
				sys.exit()
		
		else :
			raise IOError ("Unknowing format")
	
	
	# Time when the program finish to parse
	if verbose >= 2:
		TimeWarn("Parsing Finished")
	step1_end=time.time()
	
	#==========================================================================
	# Beginning of output file
	#==========================================================================
	
	if verbose >=3 :
		TimeWarn("Starting to write beginning of the file")
	
	# Definition of the output file name
	if outfile :
		name = outfile
		try:
			filehandle = open(name,'w') # writing the beginning of the outfile
		except IOError:
			sys.stderr.write("Invalid path for Outfile: No such file or directory\n\n")
			filehandle = sys.stdout
			outfile = False
	
	# By default write in the standard output
	else:
		filehandle = sys.stdout
	
	
	
	
	
	# writing the beginning of the output file	
	filehandle.write("; crer-scan  -i %s " % args.file)
	if sort:
		filehandle.write("-s ")
	if pre_table:
		filehandle.write("-pre_table ")
	filehandle.write("-in_format %s " % args.format)
	filehandle.write("-lth_crer_size %g " % args.lth_crer_size)
	filehandle.write("-uth_crer_size %g " % args.uth_crer_size)
	filehandle.write("-lth_crer_sites %g " % args.lth_crer_sites)
	filehandle.write("-uth_crer_sites %g " % args.uth_crer_sites)
	filehandle.write("-lth_crer_sites_distance %g " % args.lth_crer_sites_distance)
	filehandle.write("-uth_crer_sites_distance %g " % args.uth_crer_sites_distance)
	filehandle.write("-uth_crer_eval %g " % args.uth_crer_eval)
	filehandle.write("-uth_crer_pval %g " % args.uth_crer_pval)
	filehandle.write("-lth_crer_sig %s " % args.lth_crer_significance)
	filehandle.write("-uth_site_pval %g " % args.uth_site_pval)
	filehandle.write("-lth_score %g " % args.lth_site_score)
	filehandle.write("-uth_score %g " % args.uth_site_score)
	filehandle.write("-uth_overlap %g "% args.uth_overlap)
	if nopval :
		filehandle.write("-nopval " )
	if return_limits:
		filehandle.write("- return_limits ")
	if return_limits_filtered:
		filehandle.write("- return_limits_filtered ")
	filehandle.write("-o %s " % args.outfile)
	filehandle.write("-v %d \n" % verbose)


	filehandle.write("; Input file\n;	input	%s\n" % file)
	filehandle.write("; File format	%s\n" % format )
	filehandle.write(";Thresholds	lower	upper\n")
	filehandle.write(";	crer_size	%s	%s\n" % (lth_size, uth_size))
	filehandle.write(";	crer_sites	%s	%s\n" % (lth_sites, uth_sites))
	filehandle.write(";	crer_sites_distance	%s	%s\n" % (lth_dist, uth_dist))
	filehandle.write(";	crer_pval	%s\n" % uth_pval)
	filehandle.write(";	crer_sig	%s\n" % lth_sig)
	filehandle.write(";	crer_eval	%s\n" % uth_eval)
	filehandle.write(";	site_pval	Na	%s\n" % (pval_threshold))
	filehandle.write(";	site_score	%s	%s\n" % (lth_score_threshold, uth_score_threshold))
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
	
	
	
	filehandle.write("#seq_id	ft_type	ft_name	strand	start	end	hit_sum	crer_sig	crer_eval	crer_pval	hit_pval_product	weight_sum	crer_size	crer_center_size\n")
	
	if outfile:
		
		filehandle.close()
	
	if verbose >=3 :
		TimeWarn("End writing beginning of the file")
	
	#===========================================================================
	# Creation and writing of CRER
	#===========================================================================
	sum_lengths = 0
	# Print the start of the crer computation in standard err
	if verbose >= 2:
		TimeWarn("Computation of crer")
	# Create a site_set object with an empty list to put our site objects 
	T= Site_set([], False)
	step2_start=time.time()
	
	# Computation of pre p-value
	
	length_list_Site = len(list_Site)
	p = scipy.stats.binom_test(1,number_of_matrix,pval_threshold)
	
	if pre_table:
		if verbose >= 2:
			TimeWarn("Computation of table of all possible p_value")
			
		val_pval = {}
		for nb_expect_size in range(1,uth_size+1):
			N=nb_expect_size*2
			
			for nb_expect_sites in range(1,nb_expect_size+1):
				
				# Compute p-value for the number of positions with at least one site
				# pvalue positions = crer pvalue
				crer_pval = scipy.stats.binom_test(nb_expect_sites,N,p)
				val_pval[nb_expect_size,nb_expect_sites]= crer_pval

	# Add site objects in T which is a site_set object
	for a in list_Site:
		T.add_site(a)
	
	nbr_crer= 0
	nb_binom_computation = 0
	nb_sig_computation = 0
	T_count = 0
	tmp = tempfile.TemporaryFile()
	
	if verbose >=3 :
		TimeWarn("	Instanciated Site set on same seq and Site set classes")
	
	
	# Treat each list of siteset_on_same_seq object
	for j in T.liste_site_set():
		
		seq_id = j.seq_id()
		
		for seq in init_table :
			if seq[0]== seq_id:
				
				if return_limits :
					if not nopval :
						tmp.write("%s	limit	START_END	D	%s	%s	.	0	0	0	0	0	0\n" % (seq[0],seq[1], seq[2]))
					if nopval:
						if outfile:
							#opening of the file as alterable
							filehandle = open(name,'a')
						else :
							filehandle = sys.stdout
						
						filehandle.write("%s	limit	START_END	D	%s	%s	.	0	0	0	0	0	0\n" % (seq[0],seq[1], seq[2]))
				
				sum_lengths = sum_lengths + (int(seq[2]) - int(seq[1]) + 1)
				
				
		# Test if we want to sort the sites
		# and get extract the list of site objects
		if sort == True :
			
			list_Site_one_seq= j.get_sorted_sites()
			
			if verbose >=3 :
				TimeWarn("	Sort sites")
				
		else:
			list_Site_one_seq= j.get_sites()
		
		# Treat each index of site object in the site list as the first position of the crer
		for i in range(0,len(list_Site_one_seq)-1):
			# Define the beginning characteristics of the crer 
			# The number of sites is one (there is only the ith site)
			# The score sum of the crer is the score of the ith site
			# The start position is the start position of the ith site
			nb_sites= 1
			
			weight_crer = list_Site_one_seq[i].weight()
			
				
			crer_start_center = list_Site_one_seq[i].center()
			
			if format == 'ft':
				pval_prod = list_Site_one_seq[i].Pval()
				
			else:
				pval_prod = 0
			
			
			
			# Treat for each ith site, each jth site representing the following sites could be keep or not in the crer, 
			# the first is the ith+1 site.
			for j in range(i+1,len(list_Site_one_seq)):
				
				# Compute the start of the crer
				crer_start = min(list_Site_one_seq[i].start(),list_Site_one_seq[j].start())
				
				# Compute overlap between two sites as a percent of overlap based on minimal site size. 
				site_size_1 = list_Site_one_seq[j-1].end()-list_Site_one_seq[j-1].start()+1
				site_size_2 = list_Site_one_seq[j].end()-list_Site_one_seq[j].start()+1
				minimal_site_size = min (site_size_1,site_size_2)
				overlap =  (min(list_Site_one_seq[j].end(),list_Site_one_seq[j-1].end()) - max(list_Site_one_seq[j].start(),list_Site_one_seq[j-1].start())) / minimal_site_size
				
				if overlap < uth_overlap:
					# Compute the size of the crer center to center 
					# with the max of the end position between the ith site end position and the jth site end position to avoid the inclusion of jth site in ith site
					size = max(list_Site_one_seq[j].center(),list_Site_one_seq[j-1].center()) - crer_start_center +1
					
					# Verification of the size of the crer
					# Default between 30 and 500 bp
					# If the size is not correct the actual the jth site is not in the actual crer
					if lth_size<= size <= uth_size :
						
						if verbose >=3 :
							TimeWarn("	size good")
						# Compute the distance between following sites in the crer
						# Define by the start position of the jth site and the previous site, jth-1 site.
						# ( We could try to estimate the overlapping)
						dist=list_Site_one_seq[j].center() - list_Site_one_seq[j-1].center() 
						
						
						# Verification of the distance between following sites
						# Default minimal distance is 1
						# If the distance is not correct the jth site is not keeping in the crer 
						if lth_dist <= dist <= uth_dist : 
							
							if verbose >=3 :
								TimeWarn("	distance good")
							# If the distance is correct we keep this site and we increment this number
							nb_sites = nb_sites + 1 
							
							# Compute of the weight sum with addition of the jth site score
							weight_crer = weight_crer + list_Site_one_seq[j].weight()
							# Compute pval_product as multiplication of sites p-values in the crer
							# Only possible for features files because in bed files we don't have site p-values.
							if format == 'ft':
								pval_prod = pval_prod * list_Site_one_seq[j].Pval()
							# Verification of the number of sites 
							# Default minimal value is 2
							# if it is not correct we don't keep the jth site in the crer
							if lth_sites <= nb_sites <= uth_sites :
								
								if verbose >=3 :
										TimeWarn("	Nb_sites good")
								
								# Compute of the density of the crer
								
								density_crer = nb_sites/size
								
								# Compute the expected number of sites
								exp_sites = pval_prod * size
							
								if nb_sites >= exp_sites:
									
									N=size*2
									
									# Probability to observe a hit for at least one matrix at a given position*
									if not nopval :
										try:
										
											if pre_table:
												crer_pval = val_pval[size,nb_sites]
											
											else:
												# Compute p-value for the number of positions with at least one site
												# pvalue positions = crer pvalue
												crer_pval = scipy.stats.binom_test(nb_sites,N,p)
												
												nb_binom_computation += 1
											
												if verbose >=3 :
												
													TimeWarn("		%d binomial computation" % nb_binom_computation)
	
											
										except (RuntimeWarning,ValueError):
											if verbose >= 3 :
												sys.stderr.write("Error in binomial and significance computation\n")
											crer_sig = 0
											crer_pval = 0
										
										if crer_pval <= uth_pval and crer_pval<= uth_eval:
											
											# Writing of this crer in the output file
											crer_end = max(list_Site_one_seq[j].end(),list_Site_one_seq[i].end())
											crer_end_center = max(list_Site_one_seq[j].center(),list_Site_one_seq[i].center())
											
											crer_size = crer_end-crer_start+1
										
											T_count += 1
											nbr_crer += 1
											data = "%s	%s	%s	%d	%d	%d	%.2g	%.2g	%.2g	%.2g	%d	%d\n" % (list_Site_one_seq[i].seq_id(),"crer","DR", crer_start,crer_end,nb_sites,0,crer_pval,pval_prod,weight_crer,crer_size,size)
										
												
											tmp.write(data)
										
											if verbose >=3 :
												TimeWarn("	Printing %d crer in tmp file" % nbr_crer)
										
									if nopval :
										crer_sig = 0
										crer_pval = 0
										nb_binom_computation = 0
										e_val = 0
										crer_end = max(list_Site_one_seq[j].end(),list_Site_one_seq[i].end())
										crer_end_center = max(list_Site_one_seq[j].center(),list_Site_one_seq[i].center())
											
										crer_size = crer_end-crer_start+1
										
										crer = CRER(list_Site_one_seq[i].seq_id(),"crer","DR", crer_start,crer_end,nb_sites,crer_sig,e_val,crer_pval,pval_prod,weight_crer,crer_size,size)
										nbr_crer += 1
										
										if verbose >=3 :
										
											TimeWarn("	Instantiated crer object")
										
										if outfile:
											#opening of the file as alterable
											filehandle = open(name,'a')
										else :
											filehandle = sys.stdout
										
										if return_limits_filtered:
											for seq in init_table :
												if seq[0]== crer.seq_id():
													filehandle.write("%s	limit	START_END	D	%s	%s	.	0	0	0	0	0	0	0\n" % (seq[0],seq[1], seq[2]))
													init_table.remove(seq)
										
										write_crer(crer,filehandle)
										
										
	if not nopval:
		tmp.seek(0)
		read = tmp.read().strip().split('\n')
		nbr_crer = 0
	
		try:
			for ligne in read :
				crer_tmp = ligne.split('\t')
			
				if crer_tmp[1] == 'limit':
					if outfile:
						filehandle = open(name,'a')#opening of the file as alterable
					else :
						filehandle = sys.stdout
					
					filehandle.write("%s	limit	START_END	D	%s	%s	.	0	0	0	0	0	0	0\n" % (crer_tmp[0],crer_tmp[4], crer_tmp[5]))
			
				else:
					e_val = float(crer_tmp[7])* T_count
					
					# For this moment we define sig_base as one and after we can compute this.
					# Computation of crer significance
					crer_sig = -math.log10(e_val)
				
					if verbose >=3 :
						nb_sig_computation += 1
						TimeWarn("		%d significance computation" % nb_sig_computation)
					
					if lth_sig <= crer_sig :
						crer = CRER(crer_tmp[0],crer_tmp[1],crer_tmp[2],crer_tmp[3],crer_tmp[4],crer_tmp[5],crer_sig,e_val,crer_tmp[7],crer_tmp[8],crer_tmp[9],crer_tmp[10],crer_tmp[11])
						
						nbr_crer += 1
					
						if verbose >=3 :
					
							TimeWarn("	Instantiated crer object")
					
						if outfile:
						
							#opening of the file as alterable
							filehandle = open(name,'a')
						else :
						
							filehandle = sys.stdout
					
						if return_limits_filtered:
							for seq in init_table :
								if seq[0]== crer.seq_id():
									filehandle.write("%s	limit	START_END	D	%s	%s	.	0	0	0	0	0	0	0\n" % (seq[0],seq[1], seq[2]))
									init_table.remove(seq)
					
						write_crer(crer,filehandle)
					
						# Print the time every five crer on standard err
						if verbose == 2 and nbr_crer%5 == 0:
							TimeWarn("	Printing %d crer" % nbr_crer)
						
						if verbose >=3 :
							TimeWarn("	Printing %d crer" % nbr_crer)
				
		except (IndexError):
			sys.stderr.write("No Crer found\n")
					
				
				

	if outfile:
		filehandle = open(name,'a')#opening of the file as alterable
	else :
		filehandle = sys.stdout
		
		
	if verbose >=3 :
		TimeWarn("	Computation of crer finished")
	
	if verbose >=3 :
		TimeWarn("Writing the end of output file")
	
	# Time when crer computation is finished
	step2_end=time.time()

	
	nb_binom_max= uth_size*(uth_size+1)/2
	
	
	filehandle.write("; Number of sites scanned	%d\n" % length_list_Site)
	filehandle.write("; Sum of sequence lengths	%d\n" % sum_lengths)
	filehandle.write("; Number of Crer	%d\n" % nbr_crer)
	filehandle.write("; Number of binomial	%d\n"% nb_binom_computation)
	filehandle.write("; Maximal number of binomial	%d\n"% nb_binom_max)
	
	time_end= time.strftime("%Y-%d-%m.%H:%M:%S")
	filehandle.write("; Job started	%s\n" % time_start)
	filehandle.write("; Job done	%s\n" % time_end)
	seconds_end=time.time() 
	filehandle.write("; Seconds	%f\n" % (seconds_end - seconds_start))

	if outfile:
		filehandle.close()
	
	# Print time when the script is finished on standard err
	if verbose >= 2:
		TimeWarn("Finished")

