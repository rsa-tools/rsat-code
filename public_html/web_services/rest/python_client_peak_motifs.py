import requests
import json

rsat_server = "http://rsat-tagc.univ-mrs.fr/rest/"

###### EXAMPLE FOR SENDING JSON DATA #####
#data = json.dumps({'title':'Human_from_Jaspar','markov':'auto','disco':'oligos,positions','motif_db':'jaspar_core_nonredundant_vertebrates','tmp_test_infile_URL':'http://rsat-tagc.univ-mrs.fr/rsat/demo_files/Oct4_peaks_top1000.fa','task':'purge,seqlen,composition,disco,merge_motifs,split_motifs,motifs_vs_motifs,timelog,archive,synthesis,small_summary,motifs_vs_db,scan'})
#r = requests.post(rsat_server+'peak-motifs', headers={"Content-type":"application/json"}, data=data)
#result = json.loads(r.content)
#print("Result file URL : " + result['server'])
####### END EXAMPLE

######### EXAMPLE FOR SENDING FORM-DATA INCLUDING FILES AND OTHER TYPES OF DATA #######
data = {'title':'Human_from_Jaspar', 'markov':'auto', 'disco':'oligos,positions', 'motif_db':'jaspar_core_nonredundant_vertebrates','tmp_test_infile_URL':'http://rsat-tagc.univ-mrs.fr/rsat/demo_files/Oct4_peaks_top1000.fa','task':'purge,seqlen,composition,disco,merge_motifs,split_motifs,motifs_vs_motifs,timelog,archive,synthesis,small_summary,motifs_vs_db,scan'}

#files = {}
r = requests.Request('POST', rsat_server+'peak-motifs', data=data).prepare()
s = requests.Session()
resp = s.send(r)
print(resp.text)
######## END EXAMPLE



