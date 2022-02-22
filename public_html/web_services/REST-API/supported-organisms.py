#Begin by importing the Requests module
import requests, sys

#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"

#The endpoint indicates which RSAT resource you are interested in
ext = "/supported-organisms/" ##Returns the list of organisms supported on the local RSAT or on a remove server (opion -server).

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "format" : "tab", ##Output format. supported tab, tree, html_tree, newick
        "return" : "ID", ##String. Output fields. Supported ID,name,taxid,source,last_update,nb,seq_format,up_from,up_to,taxonomy,data,genome,genome_assembly,genome_version,download_date,variant_available,variant_source,variant_date,path_to_variant_files
        "taxon" : None, ##String. Selected_taxon. Only returns organisms belonging to a selected taxon.
        "group" : "Fungi", ##String. Selected_group. Only returns organisms belonging to a selected taxonomy group.
        "source" : "NCBI-Refseq", ##String. Selected_source. Only return organisms from a user-selected source(s).
        "depth" : None, ##Integer. Depth for exploring the taxonomic tree.
        "unique_species" : None, ##Boolean. Select at most one organism per species. This option aims at avoiding to be submerged by hundreds of strains sequenced for the same species.
        "unique_genus" : None, ##Boolean. Select at most one organism per genus.
        "source" : None ##String. Return the list of organisms supported on a remote RSAT server, via the Web services interfaces.
    } 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
# print (repr(r.json))