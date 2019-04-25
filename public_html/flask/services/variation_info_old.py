from flask import Flask, jsonify, abort, request, make_response, url_for
from flask_restplus import Resource, reqparse, fields
from subprocess import check_output, Popen, PIPE
import os,sys,re
import requests

service_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(service_dir + '/../lib')
sys.path.append(service_dir + '/../')
import utils
from rest_server import app

################################################################
### Get information about polymorphic variations
@app.route('/variation-info', methods=['POST','GET'])
def get_variation_info():
    files = ''
    output_choice = 'display'
    if request.method == 'POST':
        data = request.form or request.get_json(force=True)
        files = request.files
    elif request.method == 'GET':
        data = request.args
    
    command = utils.perl_scripts + '/variation-info'
    if 'h' in data: # help message and list options
        command += ' -h'
    mandatory_parameters = ['species','assembly']
    optional_parameters = ['species_suffix','release','format','type','col']
    default_param_values = {'format':'id'}
    fileupload_parameters = ['i']

    ## Read regular parameters
    parameters = utils.read_parameters(data, mandatory_parameters, optional_parameters, default_param_values)
    if parameters['error'] != 0:
        return parameters['error_message']
    command += parameters['arguments']
    
    ## Upload input file if specified
    input_files = utils.read_fileupload_parameters(data, files, fileupload_parameters, 'variation-info', '', True, ',')
    if input_files['error'] != 0:
        return input_files['error_message']
    command += input_files['arguments']

    return jsonify(utils.run_command(command, output_choice, 'variation-info', 'varBed', ''))
