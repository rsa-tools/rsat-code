import sys
sys.path.insert(0,'/usr/local/lib/python2.7/site-packages')
sys.path.insert(0,'/home/rsat/rsat/public_html/flask')
sys.path.insert(0, '/home/rsat/miniconda3/lib/python3.7/site-packages')

from rest_server import app as application
import os

os.environ['RSAT'] = '/home/rsat/rsat'
os.environ['RSA_OUTPUT_CONTEXT'] = 'RSATWS'
os.environ['rsat_www'] = 'http://rsat-tagc.univ-mrs.fr/rsat/'
