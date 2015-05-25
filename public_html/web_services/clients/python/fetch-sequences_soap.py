###########################################################'
# 
# Client to download matrix-scan results from RSAT server.
#
#
# usage: fetch-sequences_soap.py [-h] -b <BED_FILE> -g <GENOME_NAME> 
# 
# 
# Required arguments:
#   -h, --help            show this help message and exit
#
#   -b <BED_FILE>, --bed_file <BED_FILE>
#                         File in BED format.
#   -g <GENOME_NAME>, --genome <GENOME_NAME>
#                         Genome on which BED file is to be matched. Example: mm9, hg19.
# 
# Examples :
# 
# Print help
#  python fetch-sequences_soap.py -h
#
# 
# 
# Version 0.1 - 04/09/2014 - Adapted from Celine HERNANDEZ, RSAT team
# 
###########################################################'




## Tested with python 2.7.3 and suds version 0.4.1
fetchsequencesClientVersion = '0.1 - 04/09/2014'




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
    fetchsequencesClientVersion, 
    os.path.basename( __file__ ),
    platform.python_version(), 
    platform.system()
)
httpHeaders = {'User-agent': userAgent}
client.set_options(headers=httpHeaders)
client.set_options(timeout=300)

## Define a function to make a service perform the desired request using provided arguments

def call_fetch_sequences(service, args):
  result = rsat_service.fetch_sequences(args)
  return result

  
###########################################################'
## Parameters for RSAT tool fetch-sequences


## Read and process command-line

# Field to be filled
#sequence_file = ''
#matrix_file = ''
#uth_pval = ''


# See https://docs.python.org/2/library/argparse.html
import argparse

# Create the parser
parser = argparse.ArgumentParser(description='Client to download fetch-sequences results from RSAT server.', epilog='Version '+fetchsequencesClientVersion)

# List desired arguments
parser.add_argument('-b', '--bed_file', metavar='<BED_FILE>', type=argparse.FileType('r'), nargs=1, help='File in BED format.', required=True)
parser.add_argument('-g', '--genome', metavar='<GENOME_NAME>', type=str, nargs=1, help='Genome on which BED file is to be matched.', required=True)

# Parse command line
args = parser.parse_args()

# Get values to be used (read files and format uth_pval)
# Read content of the bed file
bed_filecontent = args.bed_file[0].read()
genome_4bed = args.genome[0]

# Wrap all arguments into a named list (dictionary), including some hard-coded characters
#http://rsat.ulb.ac.be/web_services/RSATWS_documentation.xml
arguments = {
		'input' : bed_filecontent,
		'genome' : genome_4bed,
		'header_format' : 'galaxy',
		'output' : 'client'
}


###########################################################'
## Perform SOAP request on RSAT server

result = call_fetch_sequences(rsat_service, arguments)




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