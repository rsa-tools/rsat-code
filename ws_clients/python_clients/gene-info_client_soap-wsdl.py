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

    organism = "Escherichia_coli_K12"
    query = ["methionine", "purine"]
    full = 0
    noquery = 0
    descr = 0
    feattype = "CDS"

    url = "http://rsat.scmbb.ulb.ac.be/rsat/web_services/RSATWS.wsdl"
    #url = "./RSATWS.wsdl" # hacked wsdl
    #url = "./RSATWS.wsdl.javaOK" # hacked wsdl

    server = SOAPpy.WSDL.Proxy(url, http_proxy = my_http_proxy)

    # get some info on server possibilities
    #methods  = server.methods.keys()
    #inpar = [(a, x.name, x.type) for a in methods for x in server.methods[a].inparams]
    #for m, name, type in inpar:
    #    print m, name, type
    #
    #outpar = [(a, x.name, x.type) for a in methods for x in server.methods[a].outparams]
    #for m, name, type in outpar:
    #    print m, name, type

    server.soapproxy.config.dumpSoapOutput = 1
    server.soapproxy.config.dumpSoapInput = 1
    server.soapproxy.config.debug = 0   

    req = GeneInfoRequest()
    req.organism = organism
    req.query = query
    req.full = 0
    req.descr = 1
    #req.feattype = "CDS"
    #req.noquery = 0

    res = server.gene_info(req)

    print res.command
    print res.client

