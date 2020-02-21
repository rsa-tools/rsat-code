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

tool = 'matrix-clustering'
### Read parameters from yaml file
(descr, get_parser, post_parser) = utils.read_parameters_from_yml(api, service_dir+'/' + tool.replace('-','_') +'.yml')

ns = api.namespace(tool, description=descr)

################################################################
### Get information about polymorphic variations
@ns.route('/',methods=['POST','GET'])
class MatrixClustering(Resource):
	@api.expect(get_parser)
	def get(self):
    		data = get_parser.parse_args()
		if data['content-type'] == 'text/plain':
			resp = self._run(data)
			return utils.output_txt(resp,200)
		return self._run(data)
	
	@api.expect(post_parser)
	def post(self):
		data = []
		if request.headers.get('Content-type') == 'application/json':
			data = request.get_json(force=True)
		else:
			data = post_parser.parse_args()
		return self._run(data)
	
	def _run(self, data):
		output_choice = 'display'
		(boolean_var, fileupload_parameters) = utils.get_boolean_file_params(service_dir+'/' + tool.replace('-','_') +'.yml')
		exclude = fileupload_parameters + ['content-type']
		for x in fileupload_parameters:
			exclude = exclude + [x + '_string', x + '_string_type']
		command = utils.perl_scripts + '/' + tool
		result_dir = utils.make_tmp_dir(tool)
		
		matrix_title = dict()
		matrix_format = dict()
		
		for param in data:
		    if 'matrix_title_' in param and data[param] != '':
		        p = re.compile('matrix_title_(\d)')
		        num = p.search(param)
		        matrix_title[int(num.group(1)) - 1] = data[param]
		    elif 'matrix_format_' in param and data[param] != '':
		        p = re.compile('matrix_format_(\d)')
		        num = p.search(param)
		        matrix_format[int(num.group(1)) - 1] = data[param]
		    elif param == 'o':
		        command += ' -o ' + result_dir + '/' + data[param]
		    elif param in boolean_var:
		        if data[param] == True:
		            command += ' -' + param
		        elif data[param] == False:
		            continue
		    elif data[param] is not None and data[param] != '' and param not in exclude:
		        if 'uth_' in param:
		            uth_type = param.split('_', 1)
		            command += ' -uth ' + uth_type[1] + ' ' + str(data[param])
		        elif 'lth_' in param:
		            lth_type = param.split('_', 1)
		            command += ' -lth ' + lth_type[1] + ' ' + str(data[param])
		        else: 
			        command += ' -' + param + ' ' + str(data[param])
		upload_str = utils.parse_fileupload_parameters(data, fileupload_parameters, tool, result_dir, ',')
		upload_str = re.sub(r"matrix_\d+", "matrix", upload_str)
		matrix_file = re.split('\s*-matrix\s*', upload_str)
		upload_str = ''
		for key in matrix_title:
		    if matrix_title[key] != '':
		        upload_str += ' -matrix ' + matrix_title[key] + ' ' + matrix_file[key+1] + ' ' + matrix_format[key]
		command += upload_str
		command += ' -outdir ' + result_dir + ' -prefix ' + tool 
		return utils.run_command_background(command, tool, result_dir, data['o']+'_SUMMARY.html')
