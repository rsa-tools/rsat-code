from flask import Flask, jsonify, abort, request, make_response, url_for
from flask_restplus import Resource, fields, reqparse
import json
from subprocess import check_output, Popen, PIPE
import os,sys,re
import requests

service_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(service_dir + '/../lib')
sys.path.append(service_dir + '/../')
import utils
from rest_server import app, api

ns = api.namespace('Supported Organisms', description='supported-organisms operations')

parser = reqparse.RequestParser()
parser.add_argument('group', type=str, help='selected_group')
parser.add_argument('format', type=str, help='output format')

ns = api.namespace('supported-organisms', description='supported-organisms operations')

@ns.route('/')
class SupportedOrganisms(Resource):
    @api.doc(description='Get supported organisms')
    @api.expect(parser)
    def get(self):
        "Get list of organisms supported by the server"
            #if request.method == 'PUT':
            #data = request.get_json(force=True) or request.form
            #elif request.method == 'GET':
        data = parser.parse_args()
        return self._call(data,'')
    
    @api.expect(parser)
    def post(self):
        "Get list of organisms supported by the server"
        #data = parser.parse_args()
        data = parser.parse_args()
        files = request.files
        return self._call(data,files)
    
    def _call(self,data,files):
        output_choice = 'display'
        command = utils.perl_scripts + '/supported-organisms'
        if 'group' in data:
            command += ' -group ' + data['group']
        if data['format'] is not None:
            command += ' -format ' + data['format']
        if 'depth' in data:
            command += ' -depth ' + data['depth']
        if 'taxon' in data:
            command += ' -taxon ' + data['taxon']
        if 'unique_species' in data:
            command += ' -unique_species'
        if 'unique_genus' in data:
            command += ' -unique_genus'
        return utils.run_command(command, output_choice, 'supported-organisms', 'txt','')
