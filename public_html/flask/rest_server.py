#!flask/bin/python
from flask import Flask, jsonify, abort, request, make_response, url_for
import json
from subprocess import check_output, Popen, PIPE
from tempfile import mkstemp, mkdtemp, NamedTemporaryFile
import os,sys,re
import stat
from time import localtime
import datetime
import pwd
import subprocess
import requests

app = Flask(__name__, static_url_path = "")

flask_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(flask_dir + '/lib')
import util

## Read environment variables
props_file = flask_dir + '/../../RSAT_config.props'
util.get_environ_vars(props_file)
os.environ['RSA_OUTPUT_CONTEXT'] = 'RSATWS'
rsat_home = os.environ['RSAT']
rsat_bin = os.environ['RSAT_BIN']
perl_scripts = rsat_home + '/perl-scripts'
public_html = rsat_home + '/public_html'

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
@app.route('/supported-organisms', methods = ['POST', 'GET'])
def get_supported_organisms():
    output_choice = 'email'
    if request.method == 'POST':
        data = request.get_json(force=True) or request.form
    elif request.method == 'GET':
        data = request.args
        output_choice = 'display'
    command = perl_scripts + '/supported-organisms'
    if 'group' in data:
        command += ' -group ' + data['group']
    if 'format' in data:
        command += ' -format ' + data['format']
    if 'depth' in data:
        command += ' -depth ' + data['depth']
    if 'taxon' in data:
        command += ' -taxon ' + data['taxon']
    if 'unique_species' in data:
        command += ' -unique_species'
    if 'unique_genus' in data:
        command += ' -unique_genus'

    return run_command(command, output_choice, 'supported-organisms', 'txt','')

#### random_seq
@app.route('/random-seq', methods = ['POST', 'GET'])
def get_random_seq():
    output_choice = 'email'
    if request.method == 'POST':
        data = request.get_json(force=True) or request.form
    elif request.method == 'GET':
        data = request.args
        output_choice = 'display'
    command = perl_scripts + '/random-seq'

    if 'h' in data: # help message
        command += ' -h '
    if 'help' in data: # list options
        command += ' -help '
    if 'l' in data: # sequence length
        command += ' -l '  + data['l']
    if 'n' in data: # number of sequences
        command += ' -n '  + data['n']

    return run_command(command, output_choice, 'random-seq', 'fasta','')

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

    command = perl_scripts + '/variation-info'
    if 'h' in data: # help message and list options
	command += ' -h'
    mandatory_parameters = ['species','assembly']
    optional_parameters = ['species_suffix','release','format','type','col']
    default_param_values = {'format':'id'}
    fileupload_parameters = ['i']

    ## Read regular parameters
    parameters = read_parameters(data, mandatory_parameters, optional_parameters, default_param_values)
    if parameters['error'] != 0:
	return parameters['error_message']
    command += parameters['arguments']

    ## Upload input file if specified
    input_files = read_fileupload_parameters(data, files, fileupload_parameters, 'variation-info', '', True, ',')
    if input_files['error'] != 0:
  	return input_files['error_message']
    command += input_files['arguments']
	
    return run_command(command, output_choice, 'variation-info', 'varBed', '' ,parameters['content_type'])

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

    command = rsat_bin + '/retrieve-variation-seq'
    if 'h' in data: # help message and list options
	command += ' -h'
    mandatory_parameters = ['species','assembly']
    optional_parameters = ['species_suffix','release','format','mml','col']
    default_param_values = {'format':'varBed', 'mml':'30', 'col':'1'}
    fileupload_parameters = ['i']

    ## Read regular parameters
    parameters = read_parameters(data, mandatory_parameters, optional_parameters, default_param_values)
    if parameters['error'] != 0:
	return parameters['error_message']
    command += parameters['arguments']

    ## Upload input file if specified
    input_files = read_fileupload_parameters(data, files, fileupload_parameters, 'retrieve-variation-seq', '')
    if input_files['error'] != 0:
        return input_files['error_message']
    command += input_files['arguments']

    return run_command(command, output_choice, 'retrieve-variation-seq', 'varSeq', '' ,parameters['content_type'])


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
    command = rsat_bin + '/variation-scan'

    if 'h' in data: # help message and list options
	command += ' -h'
    mandatory_parameters = ['m_format','bg']
    optional_parameters = ['top_matrices','mml','top_variation','lth','uth','calc_distrib','distrib_dir']
    default_param_values = {}
    fileupload_parameters = ['i','m']
    ## Read regular parameters
    parameters = read_parameters(data, mandatory_parameters, optional_parameters, default_param_values)
    if parameters['error'] != 0:
	return parameters['error_message']
    command += parameters['arguments']
    ## Upload input file if specified
    tmp_dir = make_tmp_dir('variation_scan')
    input_files = read_fileupload_parameters(data, files, fileupload_parameters, 'variation_scan', tmp_dir)
    if input_files['error'] != 0:
	if not 'u' in data:	
  	    return input_files['error_message']
    else:
	if not 'u' in data:
	    command += input_files['arguments']
	
    return run_command(command, output_choice, 'variation_scan', 'fasta', tmp_dir,parameters['content_type'])
    
#### fetch_sequence
#### return: {'output' : error, 'command' : command, 'server' : result file URL}
@app.route('/fetch-sequences', methods = ['POST', 'GET'])
def get_sequences():
    output_choice = 'email'
    files = ''
    if request.method == 'POST':
        data = request.form or request.get_json(force=True)
        files = request.files
    elif request.method == 'GET':
        data = request.args
        output_choice = 'display'
    command = perl_scripts + '/fetch-sequences'

    if 'h' in data: # help message and list options
	command += ' -h'
    mandatory_parameters = ['genome']
    optional_parameters = ['u','header_format','upstr_ext','downstr_ext','extend','reference','top','chunk']
    default_param_values = {'header_format':'UCSC'}
    fileupload_parameters = ['i']
    ## Read regular parameters
    parameters = read_parameters(data, mandatory_parameters, optional_parameters, default_param_values)
    if parameters['error'] != 0:
	return parameters['error_message']
    command += parameters['arguments']
    ## Upload input file if specified
    tmp_dir = make_tmp_dir('fetch-sequences')
    input_files = read_fileupload_parameters(data, files, fileupload_parameters, 'fetch-sequences', tmp_dir)
    if input_files['error'] != 0:
	if not 'u' in data:	
  	    return input_files['error_message']
    else:
	if not 'u' in data:
	    command += input_files['arguments']
	
    return run_command(command, output_choice, 'fetch-sequences', 'fasta', tmp_dir,parameters['content_type'])

    
### peak-motifs
### return: {'output' : errors/warnings, 'command' : command, 'server' : synthesis file URL}
@app.route('/peak-motifs', methods=['POST', 'GET'])
def peak_motifs():
    if request.method == 'POST':
        data = request.form or request.get_json(force=True)
        files = request.files
    elif request.method == 'GET':
        data = request.args
        
    command = perl_scripts + '/peak-motifs'
	### get parameters
    mandatory_parameters = []
    optional_parameters = ['max_seq_len','noov','str','title','image_format','soure','task','disco'
 'nmotifs','top_peaks','minol','maxol','markov','min_markov','max_markov','ci','r_plot']
    default_param_values = {'disco':'oligos,positions','task':'purge,seqlen,composition,disco,merge_motifs,split_motifs,motifs_vs_motifs,timelog,archive,synthesis,small_summary,motifs_vs_db,scan'}
    fileupload_parameters = ['i']
    optional_fileupload_parameters = ['ctrl','ref_motifs']
    ### Regular parameters
    parameters = read_parameters(data, mandatory_parameters, optional_parameters, default_param_values)
    if parameters['error'] != 0:
	return parameters['error_message']
    command += parameters['arguments']
    
    if 'motif_db' in data:
        motif_db = data['motif_db'].replace("'","")
        motif_db = motif_db.replace('"','')
        dbs = motif_db.split(",")
        matrix_db = supported_motif_database()
        for db in dbs:
            if db in matrix_db:
                command += " -motif_db " + matrix_db[db]['name'] + " " + matrix_db[db]['format'] + " " + public_html + "/motif_databases/" + matrix_db[db]['file']

    ## Upload input file if specified
    dir = make_tmp_dir("peak-motifs")
    input_files = read_fileupload_parameters(data, files, fileupload_parameters, 'peak-motifs', dir)
    if input_files['error'] != 0:
	return input_files['error_message']
    command += input_files['arguments']
    ## Upload optional file if specified
    #control file
    opt_files = read_fileupload_parameters(data, files, optional_fileupload_parameters, 'peak-motifs', dir, mandatory=0)
    command += opt_files['arguments']
		
	# outdir
    command += " -outdir '" + dir + "' -prefix peak-motifs &"
    result = ''
    os.system(command)
	#path to synthesis page	
    server_dir = dir + "/peak-motifs_synthesis.html"
	## change path to url
    server_dir = server_dir.replace(public_html, os.environ['rsat_www'])
    return jsonify( {'output' : result, 'command' : command, 'server' : server_dir} )

def supported_motif_database():
	matrix_db = {}
	mat_db = public_html + "/motif_databases/db_matrix_files.tab"
	with open(mat_db,"r") as fh:
		for line in fh:
			if ";" in line and line.index(";") == 0: #skip comment line
				continue
			if "#" in line and line.index("#") == 0: #skip header line
				continue
			if line in ['\n', '\r\n']: #skip empty line
				continue
			t = line.split("\t")
			db_name = t[0]
			matrix_db[db_name] = {}
			matrix_db[db_name]['name'] = t[0]
			matrix_db[db_name]['format'] = t[1]
			matrix_db[db_name]['file'] = t[2]
			matrix_db[db_name]['descr'] = t[3]
			matrix_db[db_name]['version'] = t[4]
			matrix_db[db_name]['url'] = t[5]
	return matrix_db

## Regular parameters
def read_parameters(data, mandatory, optional, default_values):
    parameters = {'error':0,
		   'error_message':'Valid parameters',
   		   'arguments':''
    }
    missing_parameters = list(set(mandatory) - set(data))
    if missing_parameters:
	parameters['error'] = 500
	parameters['error_message'] = 'Missing parameters: ' + ', '.join(missing_parameters)
    for param in mandatory + optional:
	if param in data:
	    val = data[param].replace('"','')
	    val = val.replace("'","")
	    # check whether param is boolean
	    if re.match('^true$',val.lower()):
	        parameters['arguments'] += ' -' + param
	    elif not re.match('^true$',val.lower()) and not re.match('^false$',val.lower()):
		parameters['arguments'] += ' -' + param + ' ' + data[param]
    	elif param in default_values:
	    parameters['arguments'] += ' -' + param + ' ' + default_values[param]
    ## content_type
    parameters['content_type'] = data.get('content_type','json')
    return parameters


def read_fileupload_parameters(data, files, fileupload_parameters, prefix, dir='', mandatory=True, split_char=''):
    """Upload input files from the parameters of a POST query or load
    equivalent data from the URL (GET). 

    Upload input file if specified (POST only) else check if a
    parameter with the same name is specified in the data (GET), and if
    so write a local file with its content.
    
    :param data: list. Parameters submitted with the REST query
    :param files: dictionary. Local file names for each input file
    :param fileupload_parameters: list. Names of the paramters corresponding to input files
    :param mandatory: boolean. Variable indicating whether the files are mandatory
    :param split_char: string. Character used to split the parameter from a GET arguemnt into separate lines to print in the file.

    :return: dictionary with the keys 'arguments' (arguments to add to the command-line), 'error' (error code) and 'error_message'.

    """
    input_files = {'error':0,
		   'error_message':'Input file OK',
   		   'arguments':''
    }
   
    missing_parameters = []
    for param in fileupload_parameters:
 	(fd,tmp_input_file_name) = make_tmp_file(prefix, '',dir);

        # Upload files if specified (POST)
	if param in files:
	    input_file = files[param]
	    input_file.save(tmp_input_file_name)
	    input_files['arguments'] += ' -' + param + ' ' + tmp_input_file_name

        # Get values provided in the URL and store them in file (GET)
        elif param in data:
	    input = data[param]
            #            input = "HELLO,BOUM"
#            input.strip()
#            input.replace(",", "\n")
            if (split_char != ''):
                input.replace(split_char, '\n')
            with open(tmp_input_file_name, 'w') as f:
                f.write(input)
                f.close()
	    input_files['arguments'] += ' -' + param + ' ' + tmp_input_file_name

	# Download a file specified by an URL
    	elif param+'_url' in data:
    	    readURL = requests.get(data[param+'_url'])
    	    with open(tmp_input_file_name, 'w') as f:
                f.write(readURL.text)
                f.close()
	    input_files['arguments'] += ' -' + param + ' ' + tmp_input_file_name
            
        # Get the path of a file piped from a previous query 
	elif param+'_piping' in data:
	    input_files['arguments'] += ' -' + param + ' ' + data[param+'_piping']

        # If none of the above solutions worked, indicate the missing parameter
        else:
	    missing_parameters.append(param)

    # Report the list of missing files with an error of type 500
    if missing_parameters and mandatory:
	input_files['error'] = 500
	input_files['error_message'] = 'Missing input file(s): -' + ', -'.join(missing_parameters) 

    return input_files


def run_command(command, output_choice, method_name, out_format, out_dir='', content_type='json'):
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
    (fd, temp_path) = make_tmp_file(method_name, out_format, dir=out_dir)
    os.chmod(temp_path, 0777)
    with open(temp_path, 'w') as f:
        f.write(result)
    os.close(fd)
    if error != '':
        result = 'Server Error/Warning: ' + error
    output = error
    if(output_choice == 'display'):
        output = result
	# change path to url
    output = ''
    if output_choice == 'email':
        output = error
    elif output_choice == 'display':
        output = result
    response = {}
    response['command'] = command.replace(rsat_home, '$RSAT')
    response['output'] = output
    response['result_path'] = temp_path.replace(rsat_home, '$RSAT')
    response['result_url'] =  temp_path.replace(public_html, os.environ['rsat_www'])
    if content_type == 'json':
    	return jsonify(response)
    elif content_type == 'text':
	response = make_response(output)
	response.headers['Content-type'] = 'text/plain'
	return response
    return 'This server does not support ' + content_type + ' content-type'

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
        os.chmod(tmp_dir,0777)
    pref = method_name + '_' + dt.strftime('%Y') + '-' + dt.strftime('%m') + '-' + dt.strftime('%d') + '.' + dt.strftime('%H')+dt.strftime('%M')+dt.strftime('%S') + '_'

    dir_path = mkdtemp('',pref,tmp_dir)
    os.chmod(dir_path,0777)
    return dir_path

if __name__ == '__main__':
    app.run(debug = True)
