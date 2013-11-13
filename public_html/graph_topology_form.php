<html>
<head>
   <title>Network Analysis Tools - graph-topology</title>
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
  # variable definition
  $default_degree = "checked";
  $default_betweenness = "checked";
  $default_closeness = "checked";
  # PIPE VALUES
  $pipe = $_REQUEST['pipe'];
  $graph_file = $_REQUEST['graph_file'];
  $graph_format = $_REQUEST['graph_format'];
  $scol = $_REQUEST['scol'];
  $tcol = $_REQUEST['tcol'];
  $wcol = $_REQUEST['wcol'];
  $all_nodes_selected = 'checked';
  # demo graph
  $demo = $_REQUEST['demo'];
  if ($demo == 1) {
    $demo_graph = storeFile("demo_files/protein_interactions_uetz.tab");
    $all_nodes_selected = '';
    $all_nodes_selected = 'checked';
  }

  title('graph-topology');
  echo ("<center>Calculate various node topology statistics : the node degree, the betweenness and the closeness. In directed graph, it specifies whether the nodes are only source or target nodes.\n");
   echo ("<br>This program was developed by 
	  <a target=_blank href='http://www.bigre.ulb.ac.be/people/Members/sylvain'>Sylvain Broh&eacute;e</a> and
          <a target=_blank href=http://www.bigre.ulb.ac.be/Users/jvanheld/>Jacques van Helden</a>.
          </center>\n");
  echo ("<form method='post' action='graph_topology.php' enctype='multipart/form-data'>
  &nbsp;&nbsp;&nbsp;<a href = 'help.graph_topology.html#formats'><B>Input format</B></a>");
  if (!$pipe) {
    echo ("
    &nbsp;<select name='in_format'>
    <option selected value = 'tab'> tab-delimited format
    <option value = 'adj_matrix'> Adjacency matrix
    <option value = 'gml'> GML format
    </select><br><br>");
  } else {
    echo ": $graph_format<br>";
    echo "<input type='hidden' NAME='in_format' VALUE='$graph_format'>";
  }
  if (!$pipe) {
    if ($demo) {
      demo("This demonstration graph is the yeast two-hybrid dataset produced by <a target = 'top' href = 'http://www.ncbi.nlm.nih.gov/sites/entrez?db=pubmed&uid=10688190&cmd=showdetailview&indexed=google'>Uetz et al (2001)</a>. It consists in 865 interactions between 926 proteins. " );
    }
    echo ("<b>Graph</b><br>");
    echo ("<textarea name='graph' rows='6' cols='65'>$demo_graph</textarea>
    <br>Upload graph from file : <br>
    <input type='file' name='graph_file' size='45' /><br>
    &nbsp;&nbsp;&nbsp;
    <br><a href = 'help.graph_topology.html#columns'>Column specification (only relevant for tab-delimited input)</a><br>
    <table>
    <tr><td><B><a href = 'help.graph_topology.html#scol'>Source node</a></B></td><td><input type = 'text' name='s_col' value = '$default_scol' size = 1></input></td></tr>
    <tr><td><B><a href = 'help.graph_topology.html#scol'>Target node</a></B></td><td><input type = 'text' name='t_col' value = '$default_tcol' size = 1></input></td></tr>
    <tr><td><B><a href = 'help.graph_topology.html#wcol'>Weight or label column</a></B></td><td><input type = 'text' name='w_col' size = 1 value = '$default_wcol'></input></td></tr>
    </table>");
  } else {
    info_link("Graph uploaded from the previous treatment", rsat_path_to_url($graph_file));    echo "<input type='hidden' NAME='pipe_graph_file' VALUE='$graph_file'>";
  }
  if ($graph_format == 'tab') {
    echo "<input type='hidden' NAME='s_col' VALUE='$scol'/>\n";
    echo "<input type='hidden' NAME='t_col' VALUE='$tcol'/>\n";
    echo "<input type='hidden' NAME='w_col' VALUE='$wcol'/>\n";
  }
    echo("<hr><table>
  <tr><td  colspan = 2><b><a href = 'help.graph_topology.html#return_fields'>Return fields</a></b></td><td></td></tr>
  <tr><td><input type='checkbox' name='degree' value='on' $default_degree/></td><td><B>Degree</B></td></tr>
  <tr><td><input type='checkbox' name='betweenness' value='on' $default_betweenness/></td><td><B>Betweenness</B></td><td rowspan = 2><b>Warning</b> May be time consuming depending on the number of nodes of the graph</td></tr>
  <tr><td><input type='checkbox' name='closeness' value='on' $default_closeness/></td><td><B>Closeness</B></td><tr>
  </table><br>");
  
 echo("
  <input type='checkbox' name='directed' value='on' />&nbsp;<B><a href = 'help.graph_topology.html#undirected'>Directed graph</a></B><br>");
   
  echo("<hr>
  <br><br><b><a href = 'help.graph_topology.html#nodes'>Node selection</a></b><br>
  <br><table>
  <tr><td><input type = 'radio' $all_nodes_selected name='allnodes' value = 'all'/></td><td>All nodes</td>
  <tr><td><input type='radio' $selection_nodes_selected name='allnodes' value='selection' /></td>
  <td></select>
  List of nodes<br></td></tr>
  <tr><td></td><td><textarea name='nodes' rows='6' cols='65'></textarea>
  <br>Upload nodes from file : <br>
  <input type='file' name='nodes_file' size='45' /></td></tr></table><br>");
  

  
  echo ("
  <ul><ul><table class='formbutton'>
  <TD><input type='submit' name='.submit' value='GO' /></TD>
  <TD><B><A HREF='graph_topology_form.php?demo=0'>RESET</A></B></TD>
  <TD><B><A HREF='graph_topology_form.php?demo=1'>DEMO</A></B></TD>
  </form>
  <TD><B><A HREF='help.graph_topology.html'>MANUAL</A></B></TD>
  <TD><B><A target = '_blank' HREF='".checkNeatTutorial("tutorials/neat_tutorial/Node_degree_statistics.html")."'>TUTORIAL</A></B></TD>
  <TD><B><A HREF='mailto:sylvain@bigre.ulb.ac.be'>MAIL</A></B></TD>
  </TR></TABLE></ul></ul>");


?>
</body>
</html>
