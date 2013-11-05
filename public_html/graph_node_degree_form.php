<html>
<head>
   <title>GrA-tools - graph-node-degree</title>
   <link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
</head>
<body class="form">
<?php
  require ('functions.php');
  require ('demo_dataset.php');
  # variable definition
  $default_scol = 1;
  $default_tcol = 2;
  $default_wcol = "";
  # PIPE VALUES
  $pipe = $_REQUEST['pipe'];
  $graph_file = $_REQUEST['graph_file'];
  $graph_format = $_REQUEST['graph_format'];
  $scol = $_REQUEST['scol'];
  $tcol = $_REQUEST['tcol'];
  $wcol = $_REQUEST['wcol'];
  $eccol = $_REQUEST['eccol'];
  $all_nodes_selected = 'checked';
  # demo graph
  $demo = $_REQUEST['demo'];
  if ($demo == 1) {
    $demo_graph = $uetz;
    $all_nodes_selected = '';
    $all_nodes_selected = 'checked';
  }

  title('graph-node-degree');
  echo ("<center>Calculate the node degree of each node and specifies if this node is a
    seed or a target node.\n");
   echo ("<br>This program was developed by 
	  <a target=_blank href='http://www.bigre.ulb.ac.be/people/Members/sylvain'>Sylvain Broh&eacute;e</a> and
          <a target=_blank href=http://www.bigre.ulb.ac.be/Users/jvanheld/>Jacques van Helden</a>.
          </center>\n");
  echo ("<form method='post' action='graph_node_degree.php' enctype='multipart/form-data'>
  &nbsp;&nbsp;&nbsp;<a href = 'help.graph_node_degree.html#formats'><B>Input format</B></a>");
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
    <br><a href = 'help.graph_node_degree.html#columns'>Column specification (only relevant for tab-delimited input)</a><br>
    <table>
    <tr><td><B><a href = 'help.graph_node_degree.html#scol'>Source node</a></B></td><td><input type = 'text' name='s_col' value = '$default_scol' size = 1></input></td></tr>
    <tr><td><B><a href = 'help.graph_node_degree.html#scol'>Target node</a></B></td><td><input type = 'text' name='t_col' value = '$default_tcol' size = 1></input></td></tr>
    </table>");
  } else {
    info_link("Graph uploaded from the previous treatment", rsat_path_to_url($graph_file));    echo "<input type='hidden' NAME='pipe_graph_file' VALUE='$graph_file'>";
  }
  if ($graph_format == 'tab') {
    echo "<input type='hidden' NAME='s_col' VALUE='$scol'/>\n";
    echo "<input type='hidden' NAME='t_col' VALUE='$tcol'/>\n";
  }
  echo("<hr>
  <br><br><b><a href = 'help.graph_node_degree.html#nodes'>Node selection</a></b><br>
  <br><table>
  <tr><td><input type = 'radio' $all_nodes_selected name='allnodes' value = 'all'/></td><td>All nodes</td>
  <tr><td><input type='radio' $selection_nodes_selected name='allnodes' value='selection' /></td>
  <td></select>
  List of nodes<br></td></tr>
  <tr><td></td><td><textarea name='nodes' rows='6' cols='65'></textarea>
  <br>Upload nodes from file : <br>
  <input type='file' name='nodes_file' size='45' /></td></tr></table>");
  echo ("
  <ul><ul><table class='formbutton'>
  <TD><input type='submit' name='.submit' value='GO' /></TD>
  <TD><B><A HREF='graph_node_degree_form.php?demo=0'>RESET</A></B></TD>
  <TD><B><A HREF='graph_node_degree_form.php?demo=1'>DEMO</A></B></TD>
  </form>
  <TD><B><A HREF='help.graph_node_degree.html'>MANUAL</A></B></TD>
  <TD><B><A target = '_blank' HREF='".checkNeatTutorial("tutorials/neat_tutorial/Node_degree_statistics.html")."'>TUTORIAL</A></B></TD>
  <TD><B><A HREF='mailto:sylvain@bigre.ulb.ac.be'>MAIL</A></B></TD>
  </TR></TABLE></ul></ul>");


?>
</body>
</html>
