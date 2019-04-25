from flask import Flask, jsonify, abort, request, make_response, url_for
from flask_restplus import Resource, reqparse, fields
from werkzeug.datastructures import FileStorage
import yaml
from subprocess import check_output, Popen, PIPE
import os,sys,re
import requests

service_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(service_dir + '/../lib')
sys.path.append(service_dir + '/../')
import utils
from rest_server import app,api

### Read parameters from yaml file
(descr, get_parser, post_parser) = utils.read_parameters_from_yml(api, service_dir+'/retrieve_variation_seq.yml')

ns = api.namespace('retrieve-variation-seq', description=descr)

################################################################
### Get information about polymorphic variations
@ns.route('/')
class VariationInfo(Resource):
	@api.expect(get_parser)
	def get(self):
    		data = get_parser.parse_args()
		return self._run(data)
	
	@api.expect(post_parser)
	def post(self):
		data = post_parser.parse_args()
		return self._run(data)
	
	def _run(self, data):
		output_choice = 'display'
		fileupload_parameters = ['i','m']
		command = utils.rsat_bin + '/retrieve-variation-seq'
		for param in data:
			if data[param] is not None and data[param] != '' and '_input_string' not in param and param not in fileupload_parameters:
				command += ' -' + param + ' ' + str(data[param])
		command += utils.parse_fileupload_parameters(data, fileupload_parameters, 'retrieve-variation-seq', '', True, ',')
		return utils.run_command(command, output_choice, 'variation-scan', 'varBed','')
