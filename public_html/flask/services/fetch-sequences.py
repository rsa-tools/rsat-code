from flask import Flask, jsonify, abort, request, make_response, url_for
from flask_restplus import Resource, reqparse, fields
from werkzeug.datastructures import FileStorage
import yaml
import os,sys,re
import requests

service_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(service_dir + '/../lib')
sys.path.append(service_dir + '/../')
import utils
from rest_server import app,api

## Read yaml file
#yaml_data = {}
#with open(service_dir+'/fetch-sequences.yml') as f:
#    docs = yaml.load_all(f)
#    for doc in docs:
#        for k, v in doc.items():
#                yaml_data[k] = v

#get_parser = api.parser()
#for param in yaml_data['parameters']:
#    if param['type'] != 'file':
#    	get_parser.add_argument(param['name'], type=utils.param_types[param['type']], required=param.get('required','false'), help=param.get('description',''), default=param.get('default',''))

#post_parser = api.parser()
#for param in yaml_data['parameters']:
#    if param['type'] == 'file':
#        post_parser.add_argument(param['name'], type=FileStorage, help=param.get('description',''), default=param.get('default',''), location='files')
#    else:
#        post_parser.add_argument(param['name'], type=utils.param_types[param['type']], required=param.get('required','false'), help=param.get('description',''), default=param.get('default',''), location='form')

(descr,get_parser,post_parser) = utils.read_parameters_from_yml(api,service_dir + '/fetch-sequences.yml')

ns = api.namespace('fetch-sequences', description=descr)

@ns.route('/')
class FetchSequences(Resource):
    @api.expect(get_parser)
    def get(self):
        output_choice = 'display'
        data = get_parser.parse_args()

        command = utils.perl_scripts + '/fetch-sequences'

        for param in data:
            if data[param] is not None and data[param] != '':
                command += ' -' + param + ' ' + data[param]
        
        return utils.run_command(command, output_choice, 'fetch-sequences', 'fasta', ''), 200
        
    @api.expect(post_parser)
    def post(self):
        output_choice = 'display'
        data = post_parser.parse_args()
        
        command = utils.perl_scripts + '/fetch-sequences'
        
        fileupload_parameters = ['i']
        
        for param in data:
            param_arg = param.replace('_url','').replace('_text','').replace('_piping','')
            if data[param] is not None and data[param] != '' and not param_arg in fileupload_parameters:
                command += ' -' + param + ' ' + data[param]
        command += utils.parse_fileupload_parameters(data,fileupload_parameters,'fetch-sequences','',True,',')

        return utils.run_command(command, output_choice, 'fetch-sequences', 'fasta', ''), 200
