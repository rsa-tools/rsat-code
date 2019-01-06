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

    return run_command(command, 'email', 'supported-organisms', 'txt')

@app.route('/fetch_sequences', methods = ['POST'])
def get_sequences():
    args = request.get_json(force=True)
    command = perl_scripts + '/fetch-sequences'
    if 'input' in args:
        input = args['input']
        input.strip()
        tmp_input_file = make_tmp_file('fetch-sequences', '', '')
        with open(tmp_input_file, 'w') as f:
            f.write(input)
            f.close()
    elif 'tmp_input_file' in args:
        tmp_input_file = args['tmp_input_file']
        tmp_input_file = tmp_input_file.replace("'", '')
        tmp_input_file = tmp_input_file.replace('"', '')

    tmp_input_file.strip()
    command += " -i '" + tmp_input_file + "'"

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

    return run_command(command, 'display', 'fetch-sequences', 'fasta')


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
