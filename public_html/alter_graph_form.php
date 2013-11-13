<html>
<head>
   <title>Network Analysis Tools - alter-graph</title>
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
  $default_layout = "checked";
  $default_ecolors = "";
  $default_ewidth = "";

  # PIPE VALUES
  $pipe = $_REQUEST['pipe'];
  $graph_file = $_REQUEST['graph_file'];
  $graph_format = $_REQUEST['graph_format'];
  
  $scol = $_REQUEST['scol'];
  $tcol = $_REQUEST['tcol'];
  $wcol = $_REQUEST['wcol'];  
  
  
  # demo graph
  $demo = $_REQUEST['demo'];
  if ($demo == 1) {
    $demo_graph = storeFile("demo_files/protein_interactions_uetz.tab");
    $demo_ecolors = "selected";
    $demo_ewidth = "checked";
    $default_ecolors = "";
    $default_wcol = "";
    $demo_rm_nodes = "";
    $demo_add_nodes = "5%";
    $demo_rm_edges = "";
    $demo_add_edges = "1000";
  }

  title('alter-graph');
   echo ("<center>Alter a graph either by adding or removing edges or nodes.\n");
   echo ("<br>This program was developed by 
	  <a target=_blank href='http://www.bigre.ulb.ac.be/people/Members/sylvain'>Sylvain Broh&eacute;e</a> and
          <a target=_blank href=http://www.bigre.ulb.ac.be/Users/jvanheld/>Jacques van Helden</a>.
          </center>\n");
   
   echo ("<form method='post' action='alter_graph.php' enctype='multipart/form-data'>
   &nbsp;&nbsp;&nbsp;<a href = 'help.alter_graph.html#formats'><B>Input format</B></a>");

  if (!$pipe) {
    echo("&nbsp;<select name='in_format'>
    <option selected value = 'tab'> tab-delimited format
    <option value = 'adj_matrix'> Adjacency matrix
    <option value = 'gml'> GML format
    </select><br>");
  } else {
    echo ": $graph_format<br>";
    echo "<input type='hidden' NAME='in_format' VALUE='$graph_format'>";
  } 
   
  echo("<br>
  &nbsp;&nbsp;&nbsp;<B><a href = 'help.alter_graph.html#formats'>Output format</a></B>&nbsp;<select name='out_format'>
  <option selected value = 'tab'> tab-delimited format
  <option value = 'gml'> GML format
  <option value = 'adj_matrix'> Adjacency matrix
  </select><br>");
  echo("<br>");
  if (!$pipe) {
    if ($demo) {
      demo("This demonstration graph is the yeast two-hybrid dataset produced by <a target = 'top' href = 'http://www.ncbi.nlm.nih.gov/sites/entrez?db=pubmed&uid=10688190&cmd=showdetailview&indexed=google'>Uetz et al (2001)</a>. It consists in 865 interactions between 926 proteins. " );
    }
    echo ("<b>Graph</b><br>");
    echo ("<textarea name='graph' rows='6' cols='65'>$demo_graph</textarea>
    <br>Upload graph from file : <br>
    <input type='file' name='graph_file' size='45' /><br>
    &nbsp;&nbsp;&nbsp;
    <br><a href = 'help.alter_graph.html#columns'>Column specification (only relevant for tab-delimited input)</a><br>
    <table>
    <tr><td><B><a href = 'help.alter_graph.html#scol'>Source node</a></B></td><td><input type = 'text' name='s_col' value = '$default_scol' size = 1></input></td></tr>
    <tr><td><B><a href = 'help.alter_graph.html#scol'>Target node</a></B></td><td><input type = 'text' name='t_col' value = '$default_tcol' size = 1></input></td></tr>
    <tr><td><B><a href = 'help.alter_graph.html#wcol'>Weight or label column</a></B></td><td><input type = 'text' name='w_col' size = 1 value = '$default_wcol'></input></td></tr>
    </table>");
  } else {
    info_link("Graph uploaded from the previous treatment", rsat_path_to_url($graph_file));
    echo "<input type='hidden' NAME='pipe_graph_file' VALUE='$graph_file'>";
  }
  if ($graph_format == 'tab') {
    echo "<input type='hidden' NAME='s_col' VALUE='$scol'/>\n";
    echo "<input type='hidden' NAME='t_col' VALUE='$tcol'/>\n";
    echo "<input type='hidden' NAME='w_col' VALUE='$wcol'/>\n";
  }
  echo("
    <br><a href = 'help.alter_graph.html#alteration_format'>Graph alterations (values can be entered as discrete numbers or as percentages)</a><br>
    <table>
    <tr><td><B><a href = 'help.alter_graph.html#alteration_format'>Nodes to remove</a></B></td><td><input type = 'text' name='rm_nodes' value = '$demo_rm_nodes' size = 4></input></td></tr>
    <tr><td><B><a href = 'help.alter_graph.html#alteration_format'>Edges to remove</a></B></td><td><input type = 'text' name='rm_edges' value = '$demo_rm_edges' size = 4></input></td></tr>
    <tr><td><B><a href = 'help.alter_graph.html#alteration_format'>Nodes to add</a></B></td><td><input type = 'text' name='add_nodes' value = '$demo_add_nodes' size = 4></input></td></tr>
    <tr><td><B><a href = 'help.alter_graph.html#alteration_format'>Edges to add</a></B></td><td><input type = 'text' name='add_edges' value = '$demo_add_edges' size = 4></input></td></tr>
    </table>");
    echo ("<b>Targets (nodes to be removed)</b><br>");
    echo ("<textarea name='target' rows='6' cols='65'>$targets</textarea><br>
  
  
  <input type='checkbox' name='duplicate' value='on' />&nbsp;<B><a href = 'help.random_graph.html#duplicate'>Allow duplicated edges</a></B><br>
  <input type='checkbox' name='self' value='on' />&nbsp;<B><a href = 'help.random_graph.html#self'>Allow self loops (edges having the same source and target node)</a></B><br>
  <input type='checkbox' name='directed' value='on' />&nbsp;<B><a href = 'help.random_graph.html#duplicate'>Network considered as directed</a></B><br>
 
  
  
  
  
  <ul><ul><table class='formbutton'>
  <TD><input type='submit' name='.submit' value='GO' /></TD>
  <TD><B><A HREF='alter_graph_form.php?demo=0'>RESET</A></B></TD>
  <TD><B><A HREF='alter_graph_form.php?demo=1'>DEMO</A></B></TD>
  </form>

  <TD><B><A HREF='help.alter_graph.html'>MANUAL</A></B></TD>
  <TD><B><A target = '_blank' HREF='".checkNeatTutorial("tutorials/neat_tutorial/Influence_graph_alteration.html")."'>TUTORIAL</A></B></TD>

  <TD><B><A HREF='mailto:sylvain@bigre.ulb.ac.be'>MAIL</A></B></TD>
  </TR></TABLE></ul></ul>
 ");

?>
</body>
</html>
