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

################################################################
### variation-scan
@app.route('/variation-scan', methods=['POST','GET'])
def variation_scan():
    output_choice = 'display'
    files = ''
    if request.method == 'POST':
        data = request.form or request.get_json(force=True)
        files = request.files
    elif request.method == 'GET':
        data = request.args
    command = utils.rsat_bin + '/variation-scan'

    if 'h' in data: # help message and list options
        command += ' -h'
    mandatory_parameters = ['m_format','bg']
    optional_parameters = ['top_matrices','mml','top_variation','lth','uth','calc_distrib','distrib_dir']
    default_param_values = {}
    fileupload_parameters = ['i','m']
    ## Read regular parameters
    parameters = utils.read_parameters(data, mandatory_parameters, optional_parameters, default_param_values)
    if parameters['error'] != 0:
        return parameters['error_message']
    command += parameters['arguments']
    ## Upload input file if specified
    tmp_dir = utils.make_tmp_dir('variation_scan')
    input_files = utils.read_fileupload_parameters(data, files, fileupload_parameters, 'variation_scan', tmp_dir)
    if input_files['error'] != 0:
        return input_files['error_message']
    command += input_files['arguments']
    
    return jsonify(utils.run_command(command, output_choice, 'variation_scan', 'fasta', tmp_dir))
