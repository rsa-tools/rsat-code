#!flask/bin/python
from flask import Flask, jsonify, abort, request, make_response, url_for
import json
from subprocess import check_output, Popen, PIPE
from tempfile import mkstemp, mkdtemp, NamedTemporaryFile
import os
import stat
from time import localtime
import datetime
import pwd
import subprocess
import requests

app = Flask(__name__, static_url_path = "")

perl_scripts = '/home/rsat/rsat/perl-scripts'
public_html = '/home/rsat/rsat/public_html'

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

    return run_command(command, output_choice, 'supported-organisms', 'txt','')

#### fetch_sequence
#### return: {'output' : error, 'command' : command, 'server' : result file URL}
@app.route('/fetch-sequences', methods = ['POST', 'GET'])
def get_sequences():
    output_choice = 'email'
    if request.method == 'POST':
        args = request.form or request.get_json(force=True)
        files = request.files
    elif request.method == 'GET':
        args = request.args
        output_choice = 'display'
    command = perl_scripts + '/fetch-sequences'
    tmp_dir = make_tmp_dir('fetch-sequences')
    (fd,tmp_input_file_name) = make_tmp_file('fetch-sequences', '', dir=tmp_dir);
    if 'input' in args:
        input = args['input']
        input.strip()
        with open(tmp_input_file_name, 'w') as f:
            f.write(input)
            f.close()
    elif 'tmp_input_file' in files:
        tmp_input_file = files['tmp_input_file']
        tmp_input_file.save(tmp_input_file_name)
    if 'input' in args or 'tmp_input_file' in files:
        command += " -i '" + tmp_input_file_name + "'"

    if 'url' in args:
        url = args['url']
        command += " -u '" + url + "'"
    if 'genome' in args:
        genome = args['genome']
        command += " -genome '" + genome + "'"
    if 'header_format' in args:
        header = args['header_format']
        command += " -header_format '" + header + "'"
    if 'upstr_ext' in args:
        upstr_ext = args['upstr_ext']
        command += " -upstr_ext '" + upstr_ext + "'"
    if 'downstr_ext' in args:
        downstr_ext = args['downstr_ext']
        command += " -downstr_ext '" + downstr_ext + "'"
    if 'extend' in args:
        extend = args['extend']
        command += " -extend '" + extend + "'"
    if 'reference' in args:
        reference = args['reference']
        command += " -reference '" + reference + "'"
    if 'top' in args:
        top = args['top']
        command += " -top '" + top + "'"
    if 'chunk' in args:
        chunk = args['chunk']
        command += " -chunk '" + chunk + "'"

    return run_command(command, output_choice, 'fetch-sequences', 'fasta', tmp_dir)

### peak-motifs
### return: {'output' : errors/warnings, 'command' : command, 'server' : synthesis file URL}
@app.route('/peak-motifs', methods=['POST', 'GET'])
def peak_motifs():
    if request.method == 'POST':
        args = request.form or request.get_json(force=True)
        files = request.files
    elif request.method == 'GET':
        args = request.args
        
    command = perl_scripts + '/peak-motifs'
	### get parameters
	
    if 'verbosity' in args:
        verbosity = args['verbosity'].replace("'","");
        verbosity = verbosity.replace('"', '');
        command += " -v '" + verbosity +"'";
    if 'max_seq_length' in args:
        max_seq_length = args['max_seq_length'].replace("'", "")
        max_seq_length = max_seq_length.replace('"', '')
        command += " -max_seq_len '" + max_seq_length + "'"
    if 'noov' in args:
        if args['noov'] == 1:
            command += " -noov"
    if 'str' in args:
        strand = args['str']
        if isinstance(strand, int):
            if strand == 1 or strand == 2:
                command += " -" + strand + "str"
            else:
                raise Exception("str value must be 1 or 2")
    if 'motif_db' in args:
        motif_db = args['motif_db'].replace("'","")
        motif_db = motif_db.replace('"','')
        dbs = motif_db.split(",")
        matrix_db = supported_motif_database()
        for db in dbs:
            if db in matrix_db:
                command += " -motif_db " + matrix_db[db]['name'] + " " + matrix_db[db]['format'] + " " + public_html + "/motif_databases/" + matrix_db[db]['file']
	
    if 'graph_title' in args:
         graph_title = args['graph_title'].replace("'","")
         graph_title = graph_title.replace('"',"")
         command += " -title '" + graph_title + "'"
    if 'image_format' in args:
         image_format = args['image_format'].replace("'", "")
         image_format = image_format.replace('"', '')
         command += " -img_format '" + image_format + "'"
    if 'source' in args:
         source = args['source'].replace("'", "")
         source = source.replace('"', '')
         command += " -source '" + source + "'"
    if 'task' in args:
         task = args['task'].replace("'", "")
         task = task.replace('"', '')
         command += " -task '" + task + "'"
    if 'disco' in args:
         disco = args['disco'].replace("'", "")
         disco = disco.replace('"', '')
         command += " -disco '" + disco + "'"
    if 'max_motif_number' in args:
         max_motif_number = args['max_motif_number'].replace("'", "")
         max_motif_number = max_motif_number.replace('"', '')
         command += " -nmotifs '" + max_motif_number + "'"
    if 'top_peaks' in args:
         top_peaks = args['top_peaks'].replace("'", "")
         top_peaks = top_peaks.replace('"', '')
         command += " -top_peaks '" + top_peaks + "'"
    if 'min_length' in args:
         min_length = args['min_length'].replace("'", "")
         min_length = min_length.replace('"', '')
         command += " -minol '" + min_length + "'"
    if 'max_length' in args:
         max_length = args['max_length'].replace("'", "")
         max_length = max_length.replace('"', '')
         command += " -maxol '" + max_length + "'"
    if 'markov' in args:
         markov = args['markov'].replace("'", "")
         markov = markov.replace('"', '')
         command += " -markov '" + markov + "'"
    if 'min_markov' in args:
         min_markov = args['min_markov'].replace("'", "")
         min_markov = min_markov.replace('"', '')
         command += " -min_markov '" + min_markov + "'"
    if 'max_markov' in args:
         max_markov = args['max_markov'].replace("'", "")
         max_markov = max_markov.replace('"', '')
         command += " -max_markov '" + max_markov + "'"
    if 'class_int' in args:
         class_int = args['class_int'].replace("'", "")
         class_int = class_int.replace('"', '')
         command += " -ci '" + class_int + "'"
	
	#sequence file
    dir = make_tmp_dir("peak-motifs")
    if oct(os.stat(dir).st_mode) != '040777':
        os.chmod(dir, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)
    (fd,tmp_test_infile_name) = make_tmp_file('peak-motifs', '', dir=dir)
    if 'test' in args:
        test = args['test'].rstrip()
        with open(tmp_test_infile_name, 'w') as f:
            f.write(test)
            f.close()

    if 'tmp_test_infile' in files:
        tmp_test_infile = files['tmp_test_infile']
        tmp_test_infile.save(tmp_test_infile_name)

    if 'tmp_test_infile_URL' in args:
        if args['tmp_test_infile_URL'] != '':
            readURL = requests.get(args['tmp_test_infile_URL'])
            with open(tmp_test_infile_name, 'w') as f:
                f.write(readURL.text)
                f.close()
    if 'test' in args or 'tmp_test_infile' in files or 'tmp_test_infile_URL' in args:
        command += " -i '" + tmp_test_infile_name + "'"
	
	#control file
    (fd,tmp_control_infile_name) = make_tmp_file("peak-motifs-ctr", "tab", dir=dir)
    if 'control' in args:
        control = args['control'].rstrip()
        with open(tmp_control_infile_name, 'w') as f:
            f.write(control)
            f.close()
    if 'tmp_control_infile' in files:
        tmp_control_infile = files['tmp_control_infile']
        tmp_control_infile.save(tmp_control_infile_name)
        
    if 'control' in args or 'tmp_control_infile' in files:
        command += " -ctrl '" + tmp_control_infile_name + "'"
	
	# ref_motif
    tmp_ref_motif = ''
    if 'ref_motif' in args:
        ref_motif = args['ref_motif'].rstrip()
        (fd, tmp_ref_motif) = make_tmp_file("peak-motifs_ref-motifs", "tab",dir=dir)
        with open(tmp_ref_motif, 'w') as f:
            f.write(ref_motif)
            f.close()
        tmp_ref_motif = tmp_ref_motif.rstrip()
        command += " -ref_motifs '" + tmp_ref_motif + "'"	

	# outdir
    command += " -outdir '" + dir + "' -prefix peak-motifs &"
	
	## execute command
	#p = Popen(command, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
	#(child_stdin, child_stdout, child_stderr) = (p.stdin, p.stdout, p.stderr)
    result = ''
	#for line in iter(child_stdout.readline, ''):
	#	result += line
	#error = ''
    #for line in iter(child_stderr.readline, ''):
    #	if 'WARNING' in line:
    #		result += line
    #	error += line

    #if error != '':
    #		result += '\nServer Error/Warnings: ' + error
    os.system(command)
	#path to synthesis page	
    server_dir = dir + "/peak-motifs_synthesis.html"
	## change path to url
    server_dir = server_dir.replace(os.environ['RSAT'] + "/public_html", os.environ['rsat_www'])
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
			
@app.route('/oligo_analysis', methods = ['POST'])
def oligo_analysis_pipeline():
    args = request.get_json(force=True)
    ## create config file
    tmp_path = make_tmp_dir('oligo_analysis_pipeline')
    config_json = dict()
    config_json['path'] = tmp_path
    config_json['scripts'] = perl_scripts
    
    if 'genome' in args:
        genome = args['genome']
        config_json['genome'] = genome

    #### read cne_coords
    if 'cne_coords' in args:
        cne_coords = args['cne_coords']
        cne_coords.strip()
        cne_coords_file = tmp_path + '/cne_coords.txt'
        with open(cne_coords_file, 'w') as f:
            f.write(cne_coords)
            f.close()
        config_json['cne_coords'] = cne_coords_file
    elif 'sequence' in args:
        sequence = args['sequence']
        sequence.strip()
        tmp_infile = tmp_path + '/sequence.fasta'
        with open(tmp_infile, 'w') as f:
            f.write(sequence)
            f.close()
        config_json['sequence'] = tmp_infile
        config_json['cne_coords'] = ''
    elif 'tmp_infile' in args:
        tmp_infile = tmp_path + '/sequence.fasta'
        with open(tmp_infile, 'w') as f:
            with open(args['tmp_infile'], 'r') as fi:
                lines = fi.readlines()
                for x in lines:
                    f.write(x)
                fi.close()
            f.close()
        config_json['sequence'] = tmp_infile
        config_json['cne_coords'] = ''
    #### other parameters
    if 'purge_seq' in args and args['purge_seq'] == 1:
        config_json['purge_seq'] = 1
    else:
        config_json['pruge_seq'] = 0
    if 'format' in args:
        config_json['format'] = args['format']
    if 'organism' in args:
        config_json['organism'] = args['organism']
    if 'background' in args:
        config_json['background'] = args['background']
    if 'stats' in args:
        config_json['return'] = args['stats']
    if 'noov' in args and args['noov'] == 1:
        config_json['noov'] = ' -noov'
    else:
        config_json['noov'] = ''
    if 'str' in args:
        if args['str'] == 1 or args['str'] == 2:
            config_json['strand'] = ' -' + str(args['str']) + 'str'
        else:
            raise Exception('str value must be 1 or 2')
    ### list of length
    if 'length' in args:
        length = args['length']
        config_json['length'] = []
        if type(length) is list:
            for x in length:
                config_json['length'].append(x)
        else:
            config_json['length'].append(length)
    ### list of lower threshold
    if 'lth' in args:
        lth_ref = args['lth']
        lth = ''
        if type(lth_ref) is list:
            for x in lth_ref:
                x = x.replace("'", '')
                x = x.replace('"','')
                x.strip()
                (lt1, lt2) = x.split(' ')
                lth += " -lth '" + lt1 + "' '" + lt2 + "'"
        else:
            lth_ref.strip()
            (lt1, lt2) = lth_ref.split(' ')
            lth += " -lth '" + lt1 + "' '" + lt2 + "'"
        config_json['lth'] = lth
    ### list of upper threshold
    if 'uth' in args:
        uth_ref = args['uth']
        uth = ''
        if type(uth_ref) is list:
            for x in uth_ref:
                x = x.replace("'", '')
                x = x.replace('"','')
                x.strip()
                (ut1, ut2) = x.split(' ')
                uth += " -uth '" + ut1 + "' '" + ut2 + "'"
        else:
            uth_ref.strip()
            (ut1, ut2) = uth_ref.split(' ')
            uth += " -uth '" + ut1 + "' '" + ut2 + "'"
        config_json['uth'] = uth

    if 'pseudo' in args:
        config_json['pseudo'] = ' -pseudo ' + str(args['pseudo'])
    else:
        config_json['pseudo'] = ''
    if 'max_asmb_nb' in args:
        config_json['max_asmb_nb'] = args['max_asmb_nb']
    if 'maxfl' in args:
        config_json['flanks'] = args['maxfl']
    if 'min_weight' in args:
        config_json['min_weight'] = args['min_weight']
    if 'output_prefix' in args:
        config_json['output_prefix'] = args['output_prefix']

    with open(tmp_path + '/config.json', 'w') as f:
        f.write(json.dumps(config_json))
##### end config_json

    logfile = tmp_path + '/snakemake.log'

    command = perl_scripts + '/snakemake --snakefile ' + perl_scripts + '/oligo_analysis_pipeline --configfile ' + tmp_path + '/config.json --core 5 --directory ' + tmp_path + ' oligo_all 2> ' + logfile
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
    (fd, temp_path) = make_tmp_file('oligo_analysis', 'txt', tmp_path)
    with open(temp_path, 'w') as f:
        f.write(result)
    os.close(fd)

    url_path = os.environ['rsat_www'] + tmp_path.split('public_html')[1]
    return_files = ''
    return_files += 'Sequences> ' + url_path + '/sequence.fasta' + '\n'
    if 'purge_seq' in args:
        return_files += 'Purged Sequence> ' + url_path + '/sequence.fasta.purged' + '\n'
    return_files += 'Merged oligos> ' + url_path + '/oligo_analysis.tab' + '\n'
    return_files += 'Assembly> ' + url_path + '/oligo_analysis.asmb' + '\n'
    return_files += "Significance matrices> " + url_path + "/oligo_analysis_pssm_sig_matrices.tf" + "\n"
    return_files += "Rescaled significance matrices> " + url_path + "/oligo_analysis_pssm_sig_matrices_rescaled.tf" + "\n"
    return_files += "Final matrices(transfac format)> " + url_path + "/oligo_analysis_pssm_count_matrices.tf" + "\n"
    return_files += "Final matrices (tab format)> " + url_path + "/oligo_analysis_pssm_count_matrices.txt" + "\n"
    return_files += "Converted matrices> " + url_path + "/oligo_analysis_pssm_count_matrices.txt_converted.tab"
    
    return jsonify({'command':command, 'server': return_files })


def run_command(command, output_choice, method_name, out_format, out_dir):
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
    #os.chmod(temp_path, stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH)
    os.chmod(temp_path, 0o0777)
    #os.system('chmod 777 ' + temp_path)
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
    return jsonify( {'output' : output, 'command' : command, 'server' : temp_path.replace(os.environ['RSAT'] + "/public_html", os.environ['rsat_www'])} )

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
		os.chmod(tmp_dir, 0o0777)
		#os.system('chmod -R 777 ' + tmp_dir)
	pref = method_name + '_' + dt.strftime('%Y') + '-' + dt.strftime('%m') + '-' + dt.strftime('%d') + '.' + dt.strftime('%H')+dt.strftime('%M')+dt.strftime('%S') + '_'
	dir_path = mkdtemp('',pref,tmp_dir)
	#os.chmod(dir_path, stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH)
	os.chmod(dir_path, 0o0777)
	#subprocess.call(["chmod", "a+w", dir_path])
	#os.system('chmod -R 777 ' + dir_path)
	return dir_path

if __name__ == '__main__':
    app.run(debug = True)
