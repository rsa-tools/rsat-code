<html>
<head>
   <title>Network Analysis Tools - convert-graph</title>
   <link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
</head>
<body class="form">
<?php
   require ('functions.php');
//   require ('demo_dataset.php');

# variable definition
$default_scol = 1;
$default_tcol = 2;
$default_wcol = "";
$default_pathcol = "";
$default_layout = "none";
$default_ecolors = "";
$default_ewidth = "";
  
# PIPE VALUES
$pipe = $_REQUEST['pipe'];
$graph_file = $_REQUEST['graph_file'];
$graph_format = $_REQUEST['graph_format'];
  
$scol = $_REQUEST['scol'];
$tcol = $_REQUEST['tcol'];
$wcol = $_REQUEST['wcol'];  
$eccol = $_REQUEST['eccol'];
$sccol = $_REQUEST['sccol'];  
$tccol = $_REQUEST['tccol'];
$pathcol = $_REQUEST['pathcol'];
$distinct_path = $_REQUEST['distinct_path'];
  
# demo graph
$demo = $_REQUEST['demo'];
if ($demo == 1) {
  $demo_graph = storeFile("demo_files/string_yeast_gene_coexpression_names_gt_550.tab");
  $demo_ecolors = "selected";
  $demo_ewidth = "checked";
  $default_ecolors = "";
  $default_wcol = 3;
 }

title('convert-graph');

## Description
echo ("<center>Convert graphs between different formats. <br>Optionally,
	  apply a layout algorithm to position the nodes of the graph.\n");

## Authors
echo ("<br>This program was developed by
	  <a target=_blank
          href='http://www.bigre.ulb.ac.be/people/Members/sylvain'>Sylvain
          Broh&eacute;e</a> and <a target=_blank
          href=http://www.bigre.ulb.ac.be/Users/jvanheld/>Jacques van
          Helden</a>.
          </center>\n");

## Open the form
echo ("<form method='post' action='convert_graph.php'
	 enctype='multipart/form-data'>");

################################################################
## Input options
echo("<hr>");
echo("<p><b>Input options</b></p>");

## Input format   
echo("<a href='help.convert_graph.html#formats'><B>Input format</B></a>");

if (!$pipe) {
  echo("&nbsp;<select name='in_format'>
         <option selected value = 'tab'>tab-delimited
	 <option value = 'adj_matrix'>adjacency matrix
	 <option value = 'gml'>GML
	 <option value = 'path'>path(s)
   </select><br>");
 } else {
  echo ": $graph_format<br>";
  echo "<input type='hidden' NAME='in_format' VALUE='$graph_format'>";
 } 

## JvH I move output format belw, to better separate all the input parameters from all the output parameters

if ($pipe) {
  info_link("Graph uploaded from the previous treatment", rsat_path_to_url($graph_file));
  echo "<input type='hidden' NAME='pipe_graph_file' VALUE='$graph_file'>"; 
  echo "<br><table><tr><td><B><a href = 'help.convert_graph.html#wcol'>Column containing the weight or the label</a></B></td><td><input type = 'text' name='w_col' size = 1 value = ''></input></td></tr></table><br>";

} else {
  if ($demo) {
    demo("This demonstration graph consist in the top scoring edges of the yeast
         co-expression network downloaded from the integrative
         database <a href='http://string.embl.de/'
         target='top'>STRING</a> (<a
         href='http://www.ncbi.nlm.nih.gov/sites/entrez?cmd=Retrieve&db=PubMed&list_uids=17098935&dopt=AbstractPlus'
         target='top'>Von Mering et al, 2007)</a>.  The co-expression
         network contains 537 nodes and 4801 edges, weighted according
         to a STRING-specific co-expression score.");
  }

  ## Graph
  echo ("<b>Graph</b><br>");

  echo ("<textarea name='graph' rows='6' cols='65'>$demo_graph</textarea>
    <br>Upload graph from file : <br>
    <input type='file' name='graph_file' size='45' /><br>
    &nbsp;&nbsp;&nbsp;
  
    <br><a href = 'help.convert_graph.html#columns'><B>Column specification (only relevant for tab-delimited input)</b></a><br>
    <table>
    <tr><td><B><a href = 'help.convert_graph.html#scol'>Source</a></B></td><td><input type = 'text' name='s_col' value = '$default_scol' size = 1></input></td></tr>
    <tr><td><B><a href = 'help.convert_graph.html#scol'>Target</a></B></td><td><input type = 'text' name='t_col' value = '$default_tcol' size = 1></input></td></tr>
    <tr><td><B><a href = 'help.convert_graph.html#wcol'>Weight or label</a></B></td><td><input type = 'text' name='w_col' size = 1 value = '$default_wcol'></input></td></tr>
    <tr><td><B><a href = 'help.convert_graph.html#eccol'>Edge color</a></B></td><td><input type = 'text' name='ec_col' size = 1></input></td></tr>
    <tr><td><B><a href = 'help.convert_graph.html#sccol'>Source node color</a></B></td><td><input type = 'text' name='sc_col' size = 1></input></td></tr>
    <tr><td><B><a href = 'help.convert_graph.html#sccol'>Target node color</a></B></td><td><input type = 'text' name='tc_col' size = 1></input></td></tr>
    <tr><td><B><a href = 'help.convert_graph.html#sccol'>Path</a></B></td><td><input type = 'text' name='path_col' value = '$pathcol' size = 1></input> (only for path input format) </td></tr>    
    </table>");
 }


if ($graph_format == 'tab') {
  echo "<input type='hidden' NAME='s_col' VALUE='$scol'/>\n";
  echo "<input type='hidden' NAME='t_col' VALUE='$tcol'/>\n";
  echo "<input type='hidden' NAME='ec_col' VALUE='$eccol'/>\n";
  echo "<input type='hidden' NAME='sc_col' VALUE='$sccol'/>\n";
  echo "<input type='hidden' NAME='tc_col' VALUE='$tccol'/>\n";
 }
if ($graph_format == 'path') {
  echo "<input type='hidden' NAME='path_col' VALUE='$pathcol'/>\n";
  echo "<input type='hidden' NAME='distinct_path' VALUE='$distinct_path'/>\n";
 }

## Directed/undirected
echo (" <input type='checkbox' name='undirected' value='on'>&nbsp;<B><a href =
       'help.convert_graph.html#undirected'>Undirected graph</a></B>  (only
       relevant for adjacency matrix input and output)<br> ");



################################################################
## Output options

echo("<hr>");
echo("<p><b>Output options</b></p>");

## Output format   
echo("<B><a href = 'help.convert_graph.html#formats'>Output format</a></B>&nbsp;<select name='out_format'>");
echo ("<option value = 'tab'>tab-delimited
<option selected value = 'gml'>GML
<option value = 'dot'>DOT
<option value = 'adj_matrix'>adjacency matrix</select><br>");


## Layout
echo ("<br><B><a href = 'help.convert_graph.html#layout'>Layout</a></B>
    <select name='layout'>
      <option selected value='none'> No layout
      <option value='spring'>Spring embedding
      <option value='random'> Random
    </select> (only valid for GML output) <br> "
      );


if (!$pipe) {
  echo "<input type='checkbox' name='distinct_path' value='on' />&nbsp;<B><a href = 'help.convert_graph.html#distinctpath'>Distinct paths</a></B> (only for path input format)<br>";
 }

## Edge color scale
echo ("<B><a href = 'help.convert_graph.html#ecolors'>Edge color scale (according to weights)</a></B>
  <select name='ecolors'>
  <option $default_ecolors value = 'none'> None
  <option value = 'red'> Red gradient
  <option value = 'green'> Green gradient
  <option value = 'blue'> Blue gradient
  <option value = 'grey'> Grey gradient
  <option $demo_ecolors value = 'fire'> Yellow to red
  </select>");

## Edge width
echo("<br>
  <input type='checkbox' $demo_ewidth name='ewidth' value='on' />&nbsp;<B><a href = 'help.convert_graph.html#ewidth'>Edge width proportional to weight</a></B>");

## Action buttons
echo ("<hr>");
echo ("<br>
  <table class='formbutton'>
  <TD><input type='submit' name='.submit' value='GO' /></TD>
  <TD><B><A HREF='convert_graph_form.php?demo=0'>RESET</A></B></TD>
  <TD><B><A HREF='convert_graph_form.php?demo=1'>DEMO</A></B></TD>
  </form>
  <TD><B><A HREF='help.convert_graph.html'>MANUAL</A></B></TD>
  <TD><B><A target = '_blank' HREF='".checkNeatTutorial("tutorials/neat_tutorial/Network_visualization_forma.html")."'>TUTORIAL</A></B></TD>
  <TD><B><A HREF='mailto:sylvain@bigre.ulb.ac.be'>MAIL</A></B></TD>
  </TR></TABLE>
 ");

?>
</body>
</html>
