import requests
import json

rsat_server = "http://rsat-tagc.univ-mrs.fr/rest/"

#variation-info
data = {'species':'Homo_sapiens', 'assembly':'GRCh38'}
files = {'i':open('/workspace/rsat/public_html/flask/test.txt','rb')}
r = requests.Request('POST', rsat_server+'variation-info', files=files, data=data).prepare()
s = requests.Session()
resp = s.send(r)
json_res = json.loads(resp.text)
print json_res['result_path']

#retrieve-variation-seq
data['i_input_string'] = json_res['result_path']
data['i_input_string_type'] = 'piping'
data['format'] = 'varBed'
r = requests.Request('POST', rsat_server+'retrieve-variation-seq', data=data).prepare()
s = requests.Session()
resp = s.send(r)
#print resp.text
json_res = json.loads(resp.text)
print json_res['result_path']

# variation-scan
data = {'species':'Homo sapiens', 'assembly':'GRCh38','m_input_string':'http://rsat-tagc.univ-mrs.fr/flask/variation_scan_matrix.tf', 'm_format':'tf','m_input_string_type':'url', 'i_input_string':json_res['result_path'],'i_input_string_type':'piping'}
headers = {'Accept':'text/plain'}
r = requests.Request('POST', rsat_server+'variation-scan', data=data, headers=headers).prepare()
s = requests.Session()
resp = s.send(r)
#json_res = json.loads(resp.text)
#print json_res['result_url']
print resp.text
