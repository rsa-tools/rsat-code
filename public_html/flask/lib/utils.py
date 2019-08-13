import os,sys,re
import json
#from subprocess import check_output, Popen, PIPE
import subprocess
from tempfile import mkstemp, mkdtemp, NamedTemporaryFile
import datetime
import pwd
import requests
from flask import jsonify, make_response
from werkzeug.datastructures import FileStorage
import yaml
from flask_restplus import Resource, reqparse, fields, inputs

def get_environ_vars(props_file):
    """Add the RSAT environment variables defined in the RSAT_config.props file into 
    the list of the flask server's environment variables.
    
    :param props_file: path to the RSAT_config.props file
    
    """
    with open(props_file) as f:
        for line in f:
            if not re.match('^#',line) and not re.match('\n',line):
                fields = line.strip().split('=')
                if len(fields) == 2:
                    os.environ[fields[0]] = fields[1]
        
    f.close()

"""Global variables used by all the services
"""
util_dir = os.path.dirname(os.path.abspath(__file__))
props_file = util_dir + '/../../../RSAT_config.props'
get_environ_vars(props_file)

rsat_home = os.environ['RSAT']
rsat_bin = os.environ['RSAT_BIN']
perl_scripts = rsat_home + '/perl-scripts'
public_html = rsat_home + '/public_html'

param_types = {'string':str, 'int':int, 'file':FileStorage, 'boolean':inputs.boolean}

def hide_RSAT_path(path):
    """Hide the absolute RSAT path base, replace it by $RSAT alias.
        
    :param path: the path to hide the RSAT base dir
    
    :return: the path with $RSAT alias as base.
    
    """
    path = path.replace(os.environ['RSAT'],'$RSAT')
    return path

def make_url(path):
    """Make the URL link for a server path to make it accessible for the client.
    
    :param path: the path that needs to be changed into URL link
    
    :return: the URL link for the path.
    
    """
    path = path.replace(os.environ['RSAT'] + '/public_html', os.environ['rsat_www'])
    return path

def supported_motif_database():
    """Get the list of the motif databases supported by the server.
        
    Read the file db_matrix_files.tab to get all the supported motif databases.
    
    :return: dictionary with database name with its descriptions: name, format, local file, description, version, url
    
    """
    matrix_db = {}
    mat_db = os.environ['RSAT'] + '/public_html/motif_databases/db_matrix_files.tab'
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

def get_backgroundfile(org, oligo_len, background="upstream-noorf"):
	"""Call the function ExpectedFreqFile in perl-scripts to get the name of the background file.
	:param org: organism's or taxon's name
	:param oligo_len: oligo length
	:param background: type of background, for ex. upstream-noorf
	
	:return: name of background file
	"""
	command = "perl -e 'use lib \"" + perl_scripts + "/lib/\"; use RSAT::OrganismManager; $x=&RSAT::OrganismManager::ExpectedFreqFile(" + org + ","+ str(oligo_len) +",\""+ background+"\", str=>\"-1str\", noov=>\"ovlp\"); print $x;'"
	p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	file = p.stdout.readline()
	return file

def read_parameters_from_yml(api,yml_file):
	"""Read the parameters from yaml file
	"""
	get_parser = api.parser()
	post_parser = api.parser()
	descr = ''
	yaml_data = {}
	with open(yml_file) as f:
		docs = yaml.load_all(f)
		for doc in docs:
			for k,v in doc.items():
				yaml_data[k] = v
	descr = yaml_data['descr']
	for param in yaml_data['parameters']:
		if param['type'] != 'file':
			choices = ''
			if param['type'] == 'boolean':
				choices = (True,False)

			if 'choices' in param and param['choices'] != '':
				if param['type'] == 'string':
					choices = tuple(x.strip() for x in param['choices'].split(','))
				elif param['type'] == 'int':
					choices = tuple(int(x.strip()) for x in param['choices'].split(','))	
			get_parser.add_argument(param['name'], type=param_types[param.get('type','string')], choices=choices, required=param.get('required',False), help=param.get('description',''), default=param.get('default',''))
			post_parser.add_argument(param['name'], type=param_types[param.get('type','string')], choices=choices, required=param.get('required',False), help=param.get('description',''), default=param.get('default',''), location='form')
		else:
			help_str = 'Alternative to the file upload option, to enable specifying an input from another location than the client computer. The value can be the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows) or the content of the input. '
			help_str_type = 'url: URL (Web address) to the input file; piping: result file from other tool; text: input content'
			post_parser.add_argument(param['name'], type=FileStorage, required=param.get('required',False), help=param.get('description',''), location='files')
			get_parser.add_argument(param['name']+'_string', type=str, help=help_str)
			get_parser.add_argument(param['name']+'_string_type', type=str, choices=('url','piping','text'), default='url',help=help_str_type)
			post_parser.add_argument(param['name']+'_string', type=str, help=help_str, location='form')
			post_parser.add_argument(param['name']+'_string_type', type=str, choices=('url','piping','text'), default='url',help=help_str_type,location='form')
	get_parser.add_argument('content-type', type=str, help='Response content-type', default='text/plain')
	return (descr, get_parser, post_parser)

def read_parameters(data, mandatory, optional, default_values):
    """Read the parameters (except upload files) from a POST query (can be in JSON format or form-data format) or read
    equivalent data from the URL of a GET query.
    
    Check if the parameter is a boolean, if so only add it's name (and not its value) to the arguments string,
    else add both name and value to the arguments string.
    
    :param data: list. Parameters submitted with the REST query
    :param mandatory: list. Names of the mandatory parameters using by the service
    :param optional: list. Names of the optional paramters using by the service
    :param default_values: dictionary. Names of the optional parameters and its default values if specified
    
    :return: dictionary with the keys 'arguments' (arguments to add to the command-line), 'error' (error code) and 'error_message'.
    
    """
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
            parameters['arguments'] += ' -' + param + ' ' + data[param]
        elif param in default_values:
            parameters['arguments'] += ' -' + param + ' ' + default_values[param]
    ## content_type
    parameters['content_type'] = data.get('content_type','json')
    return parameters

def parse_fileupload_parameters(data, fileupload_parameters, prefix, dir='', split_char=''):
    parameters = ''
    for param in fileupload_parameters:
        (fd,tmp_input_file_name) = make_tmp_file(prefix, '',dir);
        
        # Upload files if specified (POST)
        if param in data and data[param] is not None and data[param] != '':
            input_file = data[param]
            input_file.save(tmp_input_file_name)
            parameters += ' -' + param + ' ' + tmp_input_file_name
        
        # Get values provided in the URL and store them in file (GET)
        elif param+'_string' in data and data[param+'_string'] is not None and data[param+'_string'] != '':
	    if data[param+'_string_type'] == 'text':
                input = data[param+'_string']
                if (split_char != ''):
                    input = '\n'.join(x.strip() for x in input.split(split_char))
                with open(tmp_input_file_name, 'w') as f:
                    f.write(input)
                    f.close()
                parameters += ' -' + param + ' ' + tmp_input_file_name
        
            # Download a file specified by an URL
            elif data[param+'_string_type'] == 'url':
                try:
		    readURL = requests.get(data[param+'_string'])
                    with open(tmp_input_file_name, 'w') as f:
                        f.write(readURL.text)
                    	f.close()
                	parameters += ' -' + param + ' ' + tmp_input_file_name
        	except requests.exceptions.RequestException as e:
                    parameters += ' -' + param + ' ' + data[param + '_string']
	   # Get the path of a file piped from a previous query
            elif data[param+'_string_type'] == 'piping':
                parameters += ' -' + param + ' ' + data[param+'_string']
    return parameters

def read_fileupload_parameters(data, files, fileupload_parameters, prefix, dir='', mandatory=True, split_char=''):
    """Upload input files from the parameters of a POST query or load
    equivalent data from the URL (GET) or getting file piped from another tool.
        
    Upload input file if specified (POST only) else check if a
    parameter with the same name is specified in the data (GET), and if
    so write a local file with its content, else check if the input is provided from piping from another tool, if so use that file.
        
    :param data: list. Parameters submitted with the REST query
    :param files: dictionary. Local file name for each input file submitted with the REST query
    :param fileupload_parameters: list. Names of the paramters corresponding to input files using by the service
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
                input = input.replace(split_char, '\n')
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


def make_tmp_file(method_name, out_format, dir=''):
    """Create a temporary file on the server.
    
    :param method_name: name of the tool (service) that uses the temporary file, used as prefix for the file's name
    :param out_format: format of the temporary file
    :dir: directory in which the temporary file will be created. If not specified, a temporary directory will be created. See @make_tmp_dir.
    
    :return: file handler and file name.
    
    """
    if dir == '':
        tmp_dir = make_tmp_dir(method_name)
    else:
        tmp_dir = dir
    dt = datetime.datetime.now()
    pref = method_name + '_' + dt.strftime('%Y') + '-' + dt.strftime('%m') + '-' + dt.strftime('%d') + '.' + dt.strftime('%H')+dt.strftime('%M')+dt.strftime('%S') + '_'

    fd, temp_path = mkstemp('.' + out_format, pref, tmp_dir)
    return (fd, temp_path)

def make_tmp_dir(method_name):
    """Make a temporary directory on the server and make it readable and writable.
        
    :param method_name: name of the tool (service) that uses the temporary directory, used as prefix for the directory's name
    
    :return: name of the temporary directory.
    
    """
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

def run_command(command, output_choice, method_name, out_format, out_dir=''):
    """Execute the command of the service requested by a REST query.
        
    :param command: full command (name and arguments) to be executed
    :output_choice: specify what to return for the 'output' key of the response, can be the result or the warning/error produced by the command.
    :method_name: name of the tool (service)
    :out_format: format of the output file
    :out_dir: directory for the output files and folders. If not specified, a temporary directory will be created
    :content_type: content-type of the response return to the client. Default: json
    
    :return: a JSON object containing the full command, the result or the warning/error of the command and the link/path to it if the content-type is json. Or the result of the command if the content-type is text.
    
    """
    p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
    (child_stdin, child_stdout, child_stderr) = (p.stdin, p.stdout, p.stderr)
    result = ''
    for line in iter(child_stdout.readline, ''):
        result += line
    error = ''
    for line in iter(child_stderr.readline, ''):
        if 'WARNING' in line:
            result += line
        else:
	    error += line
    
    #### write to file
    (fd, temp_path) = make_tmp_file(method_name, out_format, dir=out_dir)
    os.chmod(temp_path, 0777)
    with open(temp_path, 'w') as f:
        f.write(result)
    os.close(fd)
    if error != '':
        result += error + '\n'
    if(output_choice == 'display'):
        output = result
    
    # change path to url
    output = ''
    if output_choice == 'email':
        output = error
    elif output_choice == 'display':
        output = result
    response = {}
    response['command'] = hide_RSAT_path(command)
    response['output'] = output
    response['result_path'] = hide_RSAT_path(temp_path)
    response['result_url'] = make_url(temp_path)
    return response

def output_txt(data,code,headers=None):
	"""
	Output format as text. Header's content-type:text/plain 
	"""
	resp = make_response(data['output'],code)
	resp.headers.extend(headers or {'Content-type':'text/plain'})
	return resp
