
# coding: utf-8

## Python client for infer-operons

# This iPhtyon notebook presents a Python client for the SOAP/WSDL Web services supported by the Regulatory Sequences Analysis Tools (RSAT, http://www.rsat.eu/).

### Define a function that collects operons from the remote RSAT server

# In[2]:

#!/usr/bin/python

################################################################
## This script provides an example of the invocation of the Web
## service interface to the Regulatory Sequence Analysis Tools (RSAT;
## http://rsat.eu/).
##
## The script calls the service "infer-operon"

## Let us define a function, that we will later use to address different types of queries to the RSAT server
def collect_operons(organism, min_gene_nb,distance,query,all_genes,return_fields):
    """Collect the operons from a remote RSAT server, by using the SOAP/WSDL Web services to the method infer-operons.
    
    @param organism: query organism (must be supported on the remote RSAT server)
    @param min_gene_nb: minimal number of genes to report an operon 
    @param distance: threshold on the length of intergenic regions, to decide whether they are considered as within-operon or between-operson. 
    @param query: a list query of gene names
    @param all_genes: instead oof using a query list, return operons for all the genes of the query organism
    @param return_fields: a string contaning the list of (columns) fields to export, separated by commas
    @return: the method returns a table with the operon information
    """
    ## Import the Python SOAPpy library, required to use SOAP interfaces.
    import os, sys, SOAPpy
    
    ## Define the URL of the Web service
    url = "http://rsat.ulb.ac.be/rsat/web_services/RSATWS.cgi"

    ## Open a connection to the server
    server = SOAPpy.SOAPProxy(url)
         
    ## Create a dictionary (associative array), which will be used to
    ## pass all the parameters to the Web service.
    request = {'organism' : organism,
           'distance' : distance,
           'min_gene_nb' : min_gene_nb,
           'query' : query,
           'all' : all_genes,
           'return': return_fields}

    server.soapaction = 'urn:RSATWS#infer_operon' ## define the method that will be called on the server
    server.namespace =  'urn:RSATWS' ## Define the name space
    
    ## Send the request to the server and collect the result
    result = server.infer_operon(request)

    ## Print out the details of the command that was executed on the remote server
    print "; Server command: " + result.command

    ## Print out the result (the operon description table)
    return result.client


# Maintnenant que nous avons défini la fonction "collect_operons()", nous pouvons l'utilser à multiple reprises, dans des situations diverses (organismes différents, gènes de requête,...).

### Infer operons for a selected gene

# A first way to use the function collect_operons() is to specify a list of query genes.
# 

# In[3]:

result_string = collect_operons("Bacillus_subtilis_168_uid57675", 1, 55, ["hisB", "epsJ"], 0, "query,name,operon,upstr_dist,gene_nb")
print result_string


### Collect a complete table of operons for all the genes of the query organism

# In[4]:

result_string = collect_operons("Bacillus_subtilis_168_uid57675", 1, 55, [], 1, "query,name,operon,upstr_dist,gene_nb")
print result_string


# We will now draw a histogram of the operon lenghts. For this, we need to extract the information from the table.
# - ignore the comment line (starging with ';')
# - ignore the header line (starting with #)
# - cut the 5th column, containing the gene_nb field.
# 
# Documentation about the hist() function:
# http://matplotlib.org/api/pyplot_api.html

# In[43]:

genes_per_operon = {}

## Split the string 
line_nb = 0
for line in result_string.splitlines() :
    line_nb = line_nb + 1
    ## Skip comment lines
    if ((line[0] != ";") & (line[0] != '#')) :
        fields = line.split("\t")
        operon = fields[2]
        gene_nb = int(fields[4])
        genes_per_operon[operon] = gene_nb

#print "Collected operon lengths\t" + str(len(genes_per_operon.values))
gene_numbers = genes_per_operon.values()

print gene_numbers[0:10]


# In[66]:

import matplotlib.pyplot as plt
import numpy as np



max_gene_number = max(gene_numbers)
mu, sigma = 100, 15
x = mu + sigma * np.random.randn(10000)
hist, bins = np.histogram(gene_numbers, bins=max_gene_number-1)
width = 0.7 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2-0.5
plt.bar(center, hist, align='center', width=width)
plt.grid()
plt.xlabel("Number of genes within an operon")
plt.ylabel("Number of operons")
plt.show()


# In[39]:

import matplotlib.pyplot as plt
import numpy as np

mu, sigma = 100, 15
x = mu + sigma * np.random.randn(10000)
hist, bins = np.histogram(x, bins=50)
width = 0.7 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
plt.bar(center, hist, align='center', width=width)
plt.show()


# In[ ]:



