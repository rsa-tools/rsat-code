from flask import Flask, jsonify, abort, request, make_response, url_for
import json
from subprocess import check_output, Popen, PIPE
import os,sys,re
import requests

service_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(service_dir + '/../lib')
sys.path.append(service_dir + '/../')
import utils
from rest_server import app

### retrieve variation sequence
@app.route('/retrieve-variation-seq', methods = ['POST', 'GET'])
def retrieve_variation_seq():
    files = ''
    output_choice = 'display'
    if request.method == 'POST':
        data = request.form or request.get_json(force=True)
        files = request.files
    elif request.method == 'GET':
        data = request.args
    
    command = utils.rsat_bin + '/retrieve-variation-seq'
    if 'h' in data: # help message and list options
        command += ' -h'
    mandatory_parameters = ['species','assembly']
    optional_parameters = ['species_suffix','release','format','mml','col']
    default_param_values = {'format':'varBed', 'mml':'30', 'col':'1'}
    fileupload_parameters = ['i']

    ## Read regular parameters
    parameters = utils.read_parameters(data, mandatory_parameters, optional_parameters, default_param_values)
    if parameters['error'] != 0:
        return parameters['error_message']
    command += parameters['arguments']
    
    ## Upload input file if specified
    input_files = utils.read_fileupload_parameters(data, files, fileupload_parameters, 'retrieve-variation-seq', '')
    if input_files['error'] != 0:
        return input_files['error_message']
    command += input_files['arguments']

    return jsonify(utils.run_command(command, output_choice, 'retrieve-variation-seq', 'varSeq', ''))
