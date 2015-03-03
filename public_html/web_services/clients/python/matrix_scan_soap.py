###########################################################'
# 
# Client to download matrix-scan results from RSAT server.
#
#
# usage: matrix_scan_soap.py [-h] -s <SEQUENCE_FILE> -m <MATRIX_FILE> -u
#                             <P-VALUE>
# 
# 
# Required arguments:
#   -h, --help            show this help message and exit
#   -s <SEQUENCE_FILE>, --sequence_file <SEQUENCE_FILE>
#                         Sequence in FASTA format.
#   -m <MATRIX_FILE>, --matrix_file <MATRIX_FILE>
#                         Matrix in transfac format.
#   -u <P-VALUE>, --uth_pval <P-VALUE>
#                         Upper threshold for p-value.
# 
# Examples :
# 
# Print help
#  python matrix_scan_soap.py -h
#
# Run test files
#  python matrix_scan_soap.py -m param_matrix.txt -s param_sequence.txt -u 0.001
# 
# 
# Version 0.2 - 21/07/2014 - Celine HERNANDEZ, RSAT team
# 
###########################################################'




## Tested with python 2.7.3 and suds version 0.4.1
matrixscanClientVersion = '0.2 - 21/07/2014'




###########################################################'
## Define log options for suds

# Import log package
import logging
# Configure default log package to INFO
logging.basicConfig(level=logging.INFO)
# Configure log of suds clients to DEBUG for verbose output concerning Client request
logging.getLogger('suds.client').setLevel(logging.ERROR)




###########################################################'
## Create the SOAP client to request the RSAT service


# Load Client class from suds
from suds.client import Client
# Define URL for RSAT services 
#wsdlUrl =  'http://rsat.ulb.ac.be/rsat/web_services/RSATWS.wsdl'
wsdlUrl =  'http://pedagogix-tagc.univ-mrs.fr/rsat/web_services/RSATWS.wsdl'
# Create the client
client = Client(wsdlUrl)
# Need service interface to perform requests
rsat_service = client.service

# Define client header (optional)
import os, platform
userAgent = 'RSAT-Client/v%s (%s; Python %s; %s)' % (
    matrixscanClientVersion, 
    os.path.basename( __file__ ),
    platform.python_version(), 
    platform.system()
)
httpHeaders = {'User-agent': userAgent}
client.set_options(headers=httpHeaders)
client.set_options(timeout=300)

## Define a function to make a service perform the desired request using provided arguments

def call_matrix_scan(service, args):
  result = rsat_service.matrix_scan(args)
  return result



  
###########################################################'
## Parameters for RSAT tool matrix-scan


## Read and process command-line

# Field to be filled
sequence_file = ''
matrix_file = ''
uth_pval = ''

# See https://docs.python.org/2/library/argparse.html
import argparse

# Create the parser
parser = argparse.ArgumentParser(description='Client to download matrix-scan results from RSAT server.', epilog='Version '+matrixscanClientVersion)

# List desired arguments
parser.add_argument('-s', '--sequence_file', metavar='<SEQUENCE_FILE>', type=argparse.FileType('r'), nargs=1, help='Sequence in FASTA format.', required=True)
parser.add_argument('-m', '--matrix_file', metavar='<MATRIX_FILE>', type=argparse.FileType('r'), nargs=1, help='Matrix in transfac format.', required=True)
parser.add_argument('-u', '--uth_pval', metavar='<P-VALUE>', type=float, nargs=1, help='Upper threshold for p-value.', required=True)

# Parse command line
args = parser.parse_args()

# Get values to be used (read files and format uth_pval)
sequence = args.sequence_file[0].read()
matrix = args.matrix_file[0].read()
uth = ['pval '+str(args.uth_pval[0])]


# Wrap all arguments into a named list (dictionary), including some hard-coded characters
#http://rsat.ulb.ac.be/web_services/RSATWS_documentation.xml
arguments = {
	'sequence' : sequence, 
	'matrix' : matrix,
	'uth' : uth,
	'quick' : 1,
	'str' : 2,
	'origin' : 'start',
	
	'background_input' : 1, # this option requires 'markov'
	'markov' : 1,
	
	'pseudo' : 1,
	'n_treatment' : 'score',
	'matrix_format' : 'transfac',
	'background_pseudo' : 0.01,
	#'return_fields' : 'pval'
}




###########################################################'
## Perform SOAP request on RSAT server

result = call_matrix_scan(rsat_service, arguments)




###########################################################'
## Display request results

print '###############################################'
print 'Command performed on server'
print ''
print result.command
print ''
print '###############################################'
print 'Result'
print ''
print result.client

#import sys; sys.exit(0)