import sys,os,re
flask_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(flask_dir)
sys.path.append(flask_dir+'/flask')

from rest_server import app as application
