#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"
#The endpoint indicates which RSAT resource you are interested in
ext = "/XYgraph/" ##Draws a XY plot using the numeric values in selected columns of a tab-delimited file.

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server
        "i_string_type" : "text", ##Type of information provided by the input string (URL, piping, text)
        "format" : None, ##String. Output image format. Supported: gif, eps, pdf, png, jpg
        "title1" : None, ##String. first graph title. The title string should be embedded in double quotes if is contains spaces or special chars.
        "title2" : None, ##String. Second graph title
        "xleg1" : None, ##String. First X legend.
        "xleg2" : None, ##String. Second X legend.        
        "yleg1" : None, ##String. First Y legend.        
        "yleg2" : None, ##String. Second Y legend.       
        "xmax" : None, ##Integer. Maximal value represented on X axis.
        "ymax" : None, ##Integer. Maximal value represented on Y axis.
        "xmin" : None, ##Integer. Minimal value represented on X axis.
        "ymin" : None, ##Integer. Minimal value represented on Y axis.
        "same_limits" : None, ##Boolean. Use same limits (max, min) for the X and Y axes.
        "min" : None, ##Integer. Minimal value represented on both X and Y axis (combinates the effects of -xmin and -ymin).
        "max" : None, ##Integer. Maximal value represented on both X and Y axis (combinates the effects of -xmin and -ymin).
        "xgstep1" : None, ##Integer. 1st step value for the grid across X axis.
        "ygstep1" : None, ##Integer. 1st step value for the grid across Y axis.
        "gstep1" : None, ##Integer. 1st step value for the grid across both X and Y axis.
        "xgstep2" : None, ##Integer. 2nd step value for the grid across X axis.
        "ygstep2" : None, ##Integer. 2nd step value for the grid across Y axis
        "gstep2" : None, ##Integer. 2nd step value for the grid across both X and Y axis.
        "xsize" : None, ##Integer. Size of the X axis (in pixels). Default is 400.
        "ysize" : None, ##Integer. Size of the Y axis (in pixels). Default is 400.
        "size" : None, ##Integer. Size of the X and Y axis (in pixels).
        "pointsize" : None, ##Integer. Point size (in pixels).
        "columns" : None, ##Boolean. Data fields are in columns (default)
        "rows" : None, ##Boolean. Data fields are in row
        "xcol" : None, ##String. Column containing data for the X axis.
        "ycol" : None, ##String. Column containing data for the Y axis.
        "xlog" : None, ##Integer. X data are displayed on a logarithmic scale
        "ylog" : None, ##Integer. Y data are displayed on a logarithmic scale
        "log" : None, ##Integer. Same as combining -xlog N -ylog M
        "lines" : None, ##Boolean. Points are jointed by lines.
        "force_lines" : None, ##Boolean. Points are jointed by lines.
        "line" : None, ##Integer. Same as lines, but for the Nth column only
        "header" : None, ##Boolean. First line of the data file contains a column header
        "legend" : None, ##Boolean. Use the content of the first line from input file as legend for Y data.
        "histo" : None, ##Boolean. Histogram. The X data should in this case contain the middle position of ach class, and Y data the frequencies
        "fhisto" : None, ##Boolean. filled histogram.
        "hbox" : None, ##String. Highlight box. the coordinates of the highlighted box in pixels (left, top, right, bottom respectively).
        "tbox" : None, ##String. Threshold box. units of X and Y data (low_x, high_x, low_y, high_y respectively).
        "bg" : None, ##String. Background color.
        "mono" : None, ##Boolean. Monochrome. All dots are drawn in black, and a specific symbol is associated to each.
        "htmap" : None ##Boolean. An HTML document is automatically generated
        "lc" : None, ##Boolean. Label column
        "colors_string" : None, ##String. Input string specifying the query.
        "colors_string_type" : None, ##Type of information provided by the input string. Available values : text, url, piping
        "export_colors" : None, ##String. Color_file. Only working with eps and pdf formats
        "gp" : None, ##String. gnuplot additional commands. Ex. -gp ‘set size ratio 0.5’
        "hline" : None, ##String. Draw an horizontal line on the indicated position. Ex. -hline red 0 or -hline red 0,10,20,30
        "null" : None, ##Boolean. Indicate the NULL or NA character in the table to omit them.
        "vline" : None ##String. Draw a vertical line on the indicated position. Ex. -vline green 0 or -vline green 1,2,4,8,16
} 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
# print (repr(r.json))
 
