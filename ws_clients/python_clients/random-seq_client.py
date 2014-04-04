#!/usr/bin/python

################################################################
## This script provides an example of the invocation of the Web
## service interface to the Regulatory Sequence Analysis Tools (RSAT;
## http://rsat.eu/).
##
## The script calls the service "random-sequence"

class RandomSeqRequest:

    def __init__(self):

        self.organism = None
        self.sequence_length = None
        self.repetition = None
        self.bg_model = None
        self.oligo_length = None

if __name__ == '__main__':

    ## Import the Python SOAPpy library, required to use SOAP interfaces.
    import os, sys, SOAPpy
    
    ## Optional: define a proxy for the http port (in most cases you don't need this)
    if os.environ.has_key("http_proxy"):
        my_http_proxy=os.environ["http_proxy"].replace("http://","")
    else:
        my_http_proxy=None

    ## Define the URL of the Web service
    url = "http://rsat.ulb.ac.be/rsat/web_services/RSATWS.cgi"

    ## Open a connection to the server
    server = SOAPpy.SOAPProxy(url, http_proxy = my_http_proxy)
 
    ## Specify the options for the SOAP server
    server.config.dumpSoapOutput = 0
    server.config.dumpSoapInput = 0
    server.config.debug = 0 ## Beware: the debug option is very verbosy
        
    ## Define the parameters for the random-seq command that will run
    ## on the RSAT server
    organism = "Mus_musculus_EnsEMBL"
    sequence_length = 200
    repetition = 3
    bg_model = "upstream-noorf"
    oligo_length = 3

    ## Create an associative array (dictionary), which will be used to
    ## pass all the parameters to the Web service.
    request = {'organism' : organism,
           'sequence_length' : sequence_length,
           'repetition' : repetition,
           'bg_model' : bg_model,
           'oligo_length' : oligo_length}

    server.soapaction = 'urn:RSATWS#random_seq' ## define the method that will be called on the server
    server.namespace =  'urn:RSATWS' ## Define the name space
    
    ## Send the request to the server and collect the result
    result = server.random_seq(request)

    ## Print out the details of the command that was executed on the remote server
    print "; Server command: " + result.command

    ## Print out the result (the random sequence)
    print result.client
