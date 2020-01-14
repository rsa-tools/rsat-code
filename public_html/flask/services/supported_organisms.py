from flask import Flask, jsonify, abort, request, make_response, url_for
from flask_restplus import Resource, fields, reqparse
from werkzeug.datastructures import FileStorage
import yaml
import os,sys,re
import requests

service_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(service_dir + '/../lib')
sys.path.append(service_dir + '/../')
import utils
from rest_server import app, api

tool = 'supported-organisms'
### Read parameters from yaml file
(descr, get_parser, post_parser) = utils.read_parameters_from_yml(api, service_dir+'/'+tool.replace('-','_') + '.yml')
ns = api.namespace(tool, description=descr)

@ns.route('/', methods=['POST','GET'])
class SupportedOrganisms(Resource):
    @api.doc(description='Get supported organisms')
    @api.expect(get_parser)
    def get(self):
        "Get list of organisms supported by the server"
        data = get_parser.parse_args()
        if data['content-type'] == 'text/plain':
            resp = self._run(data)
            return utils.output_txt(resp,200)
        return self._run(data)
    
    @api.expect(post_parser)
    def post(self):
        "Get list of organisms supported by the server"
        data = []
        if request.headers.get('Content-type') == 'application/json':
            data = request.get_json(force=True)
        else:
            data = post_parser.parse_args()
        return self._run(data)
    
    def _run(self,data):
        output_choice = 'display'
        command = utils.perl_scripts + '/' + tool
	for param in data:
		if data[param] == True:
			command += ' -' + param
		elif data[param] == False:
			continue
		elif data[param] is not None and data[param] != '' and param != 'content-type':
			command += ' -' + param + ' ' + str(data[param])
        return utils.run_command(command, output_choice, tool, 'txt','')
