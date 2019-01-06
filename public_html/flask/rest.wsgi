import sys
sys.path.insert(0,'/usr/local/lib/python2.7/site-packages')
sys.path.insert(0,'/Users/thnguyen/rsat/public_html/flask')

from rest_server import app as application
import os

os.environ['RSAT'] = '/Users/thnguyen/rsat'
os.environ['RSA_OUTPUT_CONTEXT'] = 'RSATWS'
os.environ['rsat_www'] = 'http://rsatlocal/rsat/'
