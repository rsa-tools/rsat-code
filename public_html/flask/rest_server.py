#!/usr/bin/env python
## #!/usr/bin/python
from flask import Flask, jsonify, abort, request, make_response, url_for, Blueprint
from flask_restplus import Api
import os,sys,re
import requests

#api_v1 = Blueprint('api', __name__, url_prefix='/v1')
#api = Api(api_v1, version='1.0', title='RSAT RESTful API', description='RSAT RESTful API')

app = Flask(__name__)
#app.register_blueprint(api_v1)
api = Api(app, title='RSAT RESTful API', description='RSAT RESTful API')

from services import *
import utils

os.environ['RSA_OUTPUT_CONTEXT'] = 'RSATWS'

api.representations['text/plain'] = utils.output_txt

if __name__ == '__main__':
    app.run(debug = True)
