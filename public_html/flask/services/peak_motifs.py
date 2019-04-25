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


### peak-motifs
### return: {'output' : errors/warnings, 'command' : command, 'server' : synthesis file URL}
@app.route('/peak-motifs', methods=['POST', 'GET'])
def peak_motifs():
    if request.method == 'POST':
        data = request.form or request.get_json(force=True)
        files = request.files
    elif request.method == 'GET':
        data = request.args
    
    command = utils.perl_scripts + '/peak-motifs'
    ### get parameters
    mandatory_parameters = []
    optional_parameters = ['max_seq_len','noov','str','title','image_format','soure','task','disco'
                           'nmotifs','top_peaks','minol','maxol','markov','min_markov','max_markov','ci','r_plot']
    default_param_values = {'disco':'oligos,positions','task':'purge,seqlen,composition,disco,merge_motifs,split_motifs,motifs_vs_motifs,timelog,archive,synthesis,small_summary,motifs_vs_db,scan'}
    fileupload_parameters = ['i']
    optional_fileupload_parameters = ['ctrl','ref_motifs']
    ### Regular parameters
    parameters = utils.read_parameters(data, mandatory_parameters, optional_parameters, default_param_values)
    if parameters['error'] != 0:
       return parameters['error_message']
    command += parameters['arguments']

    if 'motif_db' in data:
        motif_db = data['motif_db'].replace("'","")
        motif_db = motif_db.replace('"','')
        dbs = motif_db.split(",")
        matrix_db = utils.supported_motif_database()
        for db in dbs:
            if db in matrix_db:
                command += " -motif_db " + matrix_db[db]['name'] + " " + matrix_db[db]['format'] + " " + public_html + "/motif_databases/" + matrix_db[db]['file']

    ## Upload input file if specified
    dir = utils.make_tmp_dir("peak-motifs")
    input_files = utils.read_fileupload_parameters(data, files, fileupload_parameters, 'peak-motifs', dir)
    if input_files['error'] != 0:
        return input_files['error_message']
    command += input_files['arguments']
    ## Upload optional file if specified
    #control file
    opt_files = utils.read_fileupload_parameters(data, files, optional_fileupload_parameters, 'peak-motifs', dir, mandatory=False)
    command += opt_files['arguments']
    
    # outdir
    command += " -outdir '" + dir + "' -prefix peak-motifs &"
    result = ''
    os.system(command)
    #path to synthesis page
    server_dir = dir + "/peak-motifs_synthesis.html"
    ## change path to url
    server_dir = utils.make_url(server_dir)
    return jsonify( {'output' : result, 'command' : command, 'server' : server_dir} )
