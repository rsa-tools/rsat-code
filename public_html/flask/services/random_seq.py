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

#### random_seq
@app.route('/random-seq', methods = ['POST', 'GET'])
def get_random_seq():
    output_choice = 'email'
    if request.method == 'POST':
        data = request.get_json(force=True) or request.form
    elif request.method == 'GET':
        data = request.args
        output_choice = 'display'
    command = utils.perl_scripts + '/random-seq'

    if 'h' in data: # help message
        command += ' -h ' + data['h']
    if 'help' in data: # list options
        command += ' -help '
    if 'l' in data: # sequence length
        command += ' -l '  + data['l']
    if 'n' in data: # number of sequences
        command += ' -n '  + data['n']

    return utils.run_command(command, output_choice, 'random-seq', 'fasta','')
