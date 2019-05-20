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

post_parser = reqparse.RequestParser()
post_parser.add_argument('group', type=str, help='selected_group', location='form')
post_parser.add_argument('format', type=str, help='output format', location='form')

get_parser = reqparse.RequestParser()
get_parser.add_argument('group', type=str, help='selected_group')
get_parser.add_argument('format', type=str, help='output format')
get_parser.add_argument('content-type', type=str, help='content-type of response', default='text/plain')

ns = api.namespace('supported-organisms', description='supported-organisms operations')

@ns.route('/')
class SupportedOrganisms(Resource):
    @api.doc(description='Get supported organisms')
    @api.expect(get_parser)
    def get(self):
        "Get list of organisms supported by the server"
        data = get_parser.parse_args()
        if data['content-type'] == 'text/plain':
            resp = self._call(data)
            return utils.output_txt(resp,200)
        return self._call(data)
    
    @api.expect(post_parser)
    def post(self):
        "Get list of organisms supported by the server"
        data = dict()
        if request.headers.get('Content-type') == 'application/json':
            data = request.get_json(force=True)
        else:
            data = post_parser.parse_args()
        return self._call(data)
    
    def _call(self,data):
        output_choice = 'display'
        command = utils.perl_scripts + '/supported-organisms'
        if 'group' in data and data['group'] is not None:
            command += ' -group ' + data['group']
        if data['format'] is not None:
            command += ' -format ' + data['format']
        if 'depth' in data and data['depth'] is not None:
            command += ' -depth ' + data['depth']
        if 'taxon' in data and data['taxon'] is not None:
            command += ' -taxon ' + data['taxon']
        if 'unique_species' in data and data['unique_species'] is not None:
            command += ' -unique_species'
        if 'unique_genus' in data and data['unique_genus'] is not None:
            command += ' -unique_genus'
        return utils.run_command(command, output_choice, 'supported-organisms', 'txt','')
