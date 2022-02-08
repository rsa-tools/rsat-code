#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"

#The endpoint indicates which RSAT resource you are interested in
ext = "/draw-heatmap/" ##Draw a heatmap from a table

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "o" : None, ##String. Outputfile. name of the output file
        "html" : None, ##String. html_map_file. If a HTML map file is defined, draw-heatmap then produces a HTML file that loads the HEATMAP.
        "rownames" : None, ##Boolean. Use this option if the first column contain the row names.
        "no_text" : None, ##Boolean. Using this option, the values are not written in the cells of the heatmap
        "out_format" : None, ##String. Output format. Supported - png,jpeg
        "title" : None, ##String. Title for the graph (only works with option -r_plot so far). 
        "gradient" : None, ##String. Color of the intensity gradient of the heatmap. Default is grey. Supported - green, blue, red, fire, grey.
        "col_width" : None, ##Integer. Width of the columns (in pixel).
        "row_height" : None, ##Integer. Height of the rows (in pixel).
        "min" : None, ##Integer. Minimal value of the heatmap.
        "max" : None, ##Integer. Maximal value of the heatmap.
        "digits" : None, ##Boolean. Round the values to the specified number of digit
        "lines" : None, ##Boolean. Add black vertical and horizontal separations lines between the cells of the heatmap
        "chaos" : None, ##Boolean. The heatmap is a CHAOS Game Representation.
        "r_plot" : None ##Boolean. Use R to generate the heatmap, rather than using the Perl GD module.
    } 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
# print (repr(r.json))
 
