#! /usr/bin/python

class RandomSeqRequest:

    def __init__(self):

        self.organism = None
        self.sequence_length = None
        self.repetition = None
        self.bg_model = None
        self.oligo_length = None

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
        
    organism = "Mus_musculus_EnsEMBL"
    sequence_length = 2000
    repetition = 4
    bg_model = "upstream-noorf"
    oligo_length = 3

    req = {'organism' : organism,
           'sequence_length' : sequence_length,
           'repetition' : repetition,
           'bg_model' : bg_model,
           'oligo_length' : oligo_length}

    server.soapaction = 'urn:RSATWS#random_seq'
    server.namespace =  'urn:RSATWS'
    
    res = server.random_seq(req)

    print res.command
    print res.client
