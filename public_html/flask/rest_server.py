from flask import Flask, jsonify, abort, request, make_response, url_for
from flask_restplus import Api
import os,sys,re
import requests

app = Flask(__name__, static_url_path = "")
api = Api(app, title='RSAT RESTful API', description='RSAT RESTful API')

from services import *

os.environ['RSA_OUTPUT_CONTEXT'] = 'RSATWS'

@app.errorhandler(400)
def not_found(error):
    return make_response(jsonify( { 'error': 'Bad request' } ), 400)
 
@app.errorhandler(404)
def not_found(error):
    return make_response(jsonify( { 'error': 'Not found' } ), 404)

def output_txt(data,code,headers=None):
	"""
	Output format as text. Header's content-type:text/plain 
	"""
	resp = make_response(data['output'],code)
	resp.headers.extend(headers or {})
	return resp

api.representations['text/plain'] = output_txt

if __name__ == '__main__':
    app.run(debug = True)
