#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"

#The endpoint indicates which RSAT resource you are interested in
ext = "/feature-map/" ##Draws a graphical map of features (e.g. results of pattern matching) in a set of sequences

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "format" : None, ##String. Output image format. Supported:  jpg, gif, png, ps
        "from" : None, ##Integer. Lower limit of the positions represented on the graph
        "to" : None, ##Integer. Upper limit of the positions represented on the graph
        "title" : None, ##String. Generic Title for the feature map.
        "label" : None, ##String. Keylist define the info to display for each feature. Valid keys are id, strand, descr, pos
        "boxlabels" : None, ##String. Writes the label within the feature box (by default, the label is written outside of this box).
        "symbol" : None, ##Boolean. Associates a graphical symbol (i.e. rectangle, circle, buterfly, …) to each feature.
        "dot" : None, ##Boolean. A color dot is associated to each feature.
        "mlen" : None, ##Integer. Map length (in pixels). Default is 600.
        "mapthick" : None, ##Integer. Thickness refers to either width (for vertical maps) or height (horizintal maps).
        "mspacing" : None, ##Integer. The size of the border between maps (in pixel).
        "border" : None, ##Integer. Image border (default=10 pixels)
        "origin" : None, ##Integer. All coordinates are recalculated relative to origin.
        "legend" : None, ##Boolean. Draws a legend on the graph, showing the symbol associated to each distinct feature.
        "scalebar" : None, ##Boolean. Draw a scale bar on the left of the graph.
        "scalestep" : None, ##Integer. Step between annotations of the scale bar.
        "no_name" : None, ##Boolean. Do not print sequence names besides each sequence.
        "scorethick" : None, ##Boolean. Each feature is displayed with a thickness proportional to its score. Only positive scores are represented.
        "maxscore" : None, ##Integer. (Only valid when -scorethick is active). Maximal allowed score value. Higher score values are clipped for the drawing.
        "minscore" : None, ##Integer. (Only valid when -scorethick is active). minimal allowed score value. Features with smaller score are not displayed.
        "maxfthick" : None, ##Integer. Max feature thickness.
        "minfthick" : None, ##Integer. Min feature thickness.
        "htmap" : None, ##Boolean. HTML map.
        "mono" : None, ##Boolean. Monochrome palette (for printing on black/white printer).
        "aacolors" : None, ##Boolean. Amino acid specific colors. 
        "colors_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "colors_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "export_colors" : None, ##String. Export the feature-specific colors in a separate file as a color file.
        "bgcolor" : None, ##String. Background color in R,G,B format.
        "horizontal" : None, ##Boolean. Horizontal map (default).
        "vertical" : None, ##Boolean. Vertical map (default is horizontal).
        "select" : None, ##String. Id_list. Only display the features whose ID is in id_list. -select ‘gataag’,’gattag’
        "seq_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "seq_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "seqformat" : None, ##String. Format of the reference sequence file.
    } 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
# print (repr(r.json))
 
