import requests
import json

rsat_server = "http://rsat-tagc.univ-mrs.fr/rest/fetch-sequences"

###### EXAMPLE FOR SENDING JSON DATA #####
#data = json.dumps({"genome":"mm9", "header_format":"galaxy", "url":"http://rsat-tagc.univ-mrs.fr/rsat/demo_files/fetch-sequences_Schmidt_2011_mm9_CEBPA_SWEMBL_R0.12_702peaks.bed"})
#r = requests.post(rsat_server, headers={"Content-type":"application/json"}, data=args)
#result = json.loads(r.content)
#print "Result file URL : " + result['server']
####### END EXAMPLE

######### EXAMPLE FOR SENDING FORM-DATA INCLUDING FILES AND OTHER TYPES OF DATA #######
data = {'genome':'mm9', 'header_format':'galaxy'}
files = {'i': open('/home/rsat/rsat/public_html/demo_files/fetch-sequences_Schmidt_2011_mm9_CEBPA_SWEMBL_R0.12_702peaks.bed', 'rb')}
headers = {'Accept':'text/plain'}
r = requests.Request('POST', rsat_server, files=files, data=data, headers=headers).prepare()
s = requests.Session()
resp = s.send(r)
print(resp.text)
######## END EXAMPLE



