################################################################
### Generate random sequences

## Import libraries
from flask import Flask, jsonify, abort, request, make_response, url_for
from flask_restplus import Resource, reqparse, fields
from werkzeug.datastructures import FileStorage
import yaml
from subprocess import check_output, Popen, PIPE
import os,sys,re
import requests

## Specify directories
service_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(service_dir + '/../lib')
sys.path.append(service_dir + '/../')
import utils
from rest_server import app,api

## Specify RSAT tool to be used
tool = 'convert-background-model'

## Read parameters from yaml file
(descr, get_parser, post_parser) = utils.read_parameters_from_yml(api, service_dir+'/' + tool.replace('-', '_') + '.yml')

ns = api.namespace(tool, description=descr)

## Define the route of the entry point (namesace)
@ns.route('/', methods=['POST','GET'])
## @api.doc(params={'species':'Species name, e.g. Homo_sapiens','assembly':'Assembly name, e.g. GRCh38'})
class RandomSeq(Resource):
	@api.expect(get_parser)

        ## Support for GET queries
	def get(self):
        	data = get_parser.parse_args()
		if data['content-type'] == 'text/plain':
			resp = self._run(data)
			return utils.output_txt(resp,200)
		return self._run(data)
    	
        ## Support for POST queries
	@api.expect(post_parser)
	def post(self):
		data = []
		if request.headers.get('Content-type') == 'application/json':
			data = request.get_json(force=True)
		else:
			data = post_parser.parse_args()
		return self._run(data)
	
        ## Run the query
	def _run(self,data):
		output_choice = 'display'
		fileupload_parameters = ['i']
		exclude = fileupload_parameters + ['content-type']
		for x in fileupload_parameters:
			exclude = exclude + [x + '_string', x + '_string_type']
		command = utils.perl_scripts + '/' + tool
		result_dir = utils.make_tmp_dir(tool)
		
		boolean_var = ['noov', 'ovlp', '1str', '2str', 'to2str']
		for param in data:
		    if param in boolean_var:
		        if data[param] == True:
		            command += ' -' + param
		        elif data[param] == False:
		            continue
		    elif data[param] is not None and data[param] != '' and param not in exclude:
				command += ' -' + param + ' ' + str(data[param])
		command += utils.parse_fileupload_parameters(data, fileupload_parameters, tool, result_dir, ',') 
		return utils.run_command(command, output_choice, tool, 'tab', result_dir) 
