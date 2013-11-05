#! /usr/bin/python

class GeneInfoRequest:

    def __init__(self):

        self.organism = None
        self.query = None
        self.noquery = None
        self.desrc = None
        self.full = None
        self.feattype = None

if __name__ == '__main__':

    import os, sys, SOAPpy
    
    if os.environ.has_key("http_proxy"):
        my_http_proxy=os.environ["http_proxy"].replace("http://","")
    else:
        my_http_proxy=None

    url = "http://rsat.scmbb.ulb.ac.be/rsat/web_services/RSATWS.cgi"
    server = SOAPpy.SOAPProxy(url, http_proxy = my_http_proxy)

    server.config.dumpSoapOutput = 0
    server.config.dumpSoapInput = 0
    server.config.debug = 1
        
    organism = "Escherichia_coli_K12"
    query = ["metA"]

    req = {'organism' : organism,
           'query' : query,
           'full' : 0,
           'descr' : 0,
           'noquery' : 0}

    server.soapaction = 'urn:RSATWS#gene_info'
    server.namespace =  'urn:RSATWS'
    
    res = server.gene_info(req)

    print res.command
    print res.client
