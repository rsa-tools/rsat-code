#!flask/bin/python
from flask import Flask, jsonify, abort, request, make_response, url_for
import json
from subprocess import check_output, Popen, PIPE
from tempfile import mkstemp, mkdtemp
import os
from time import localtime
import datetime
import pwd

app = Flask(__name__, static_url_path = "")

perl_scripts = '/Users/thnguyen/rsat/perl-scripts'

@app.errorhandler(400)
def not_found(error):
    return make_response(jsonify( { 'error': 'Bad request' } ), 400)

@app.errorhandler(404)
def not_found(error):
    return make_response(jsonify( { 'error': 'Not found' } ), 404)

@app.route('/', methods=['GET'])
def index():
    return "Hello World"

#### supported_organisms
@app.route('/supported-organisms', methods = ['POST'])
def get_supported_organisms():
    data = request.get_json(force=True)
    command = perl_scripts + '/supported-organisms'
    if 'group' in data:
        command += ' -group ' + data['group']
    if 'format' in data:
        command += ' -format ' + data['format']
    if 'depth' in data:
        command += ' -depth ' + data['depth']
    return run_command(command, data['output'], 'supported-organisms', 'tab')


def run_command(command, output_choice, method_name, out_format):
    ### execute command
    p = Popen(command, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
    (child_stdin, child_stdout, child_stderr) = (p.stdin, p.stdout, p.stderr)
    result = ''
    for line in iter(child_stdout.readline, ''):
        result += line
    error = ''
    for line in iter(child_stderr.readline, ''):
        if 'WARNING' in line:
            result += line
        error += line
    
    #### write to file
    (fd, temp_path) = make_tmp_file(method_name, out_format, dir='')
    with open(temp_path, 'w') as f:
        f.write(result)
    os.close(fd)

    if error != '':
        result = 'Server Error: ' + error
    return jsonify( {'output' : result, 'command' : command, 'server' : temp_path} )

def make_tmp_file(method_name, out_format, dir=''):
    if dir == '':
        tmp_dir = make_tmp_dir(method_name)
    else:
        tmp_dir = dir
    
    dt = datetime.datetime.now()
    pref = method_name + '_' + dt.strftime('%Y') + '-' + dt.strftime('%m') + '-' + dt.strftime('%d') + '.' + dt.strftime('%H')+dt.strftime('%M')+dt.strftime('%S') + '_'

    fd, temp_path = mkstemp('.' + out_format, pref, tmp_dir)
    return (fd, temp_path)

def make_tmp_dir(method_name):
    user = pwd.getpwuid(os.getuid())[0]
    tmp_dir = os.environ['RSAT'] + '/public_html/tmp/' + user
    
    dt = datetime.datetime.now()
    tmp_dir += '/' + dt.strftime('%Y') + '/' + dt.strftime('%m') + '/' + dt.strftime('%d') + '/'
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
        os.chmod(tmp_dir,755)
    pref = method_name + '_' + dt.strftime('%Y') + '-' + dt.strftime('%m') + '-' + dt.strftime('%d') + '.' + dt.strftime('%H')+dt.strftime('%M')+dt.strftime('%S') + '_'

    dir_path = mkdtemp('',pref,tmp_dir)
    return dir_path



if __name__ == '__main__':
    app.run(debug = True)
