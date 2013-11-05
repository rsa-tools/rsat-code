<html>
<head>
   <title>Network Analysis Tools - graph-neighbours</title>
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
  $selection_nodes_selected = '';
  $all_nodes_selected = 'checked';
  $default_step = 1;
  # demo graph
  $demo = $_REQUEST['demo'];
  if ($demo == 1) {
    $demo_graph = storeFile("demo_files/protein_interactions_gavin_2006_names.tab");
//     $demo_graph = storeFile("demo_files/rdm_test.tab");
    $demo_nodes = $gavin_names_nodes;
    $selection_nodes_selected = 'checked';
    $all_nodes_selected = '';
  }

  title('graph-neighbours');
  echo ("<center>Extract from a graph the neighbourhood of a set of seed nodes.\n");
   echo ("<br>This program was developed by 
	  <a target=_blank href=http://www.bigre.ulb.ac.be/Users/sylvain/>Sylvain Broh&eacute;e</a> and
          <a target=_blank href=http://www.bigre.ulb.ac.be/Users/jvanheld/>Jacques van Helden</a>.
          </center>\n");
  echo ("<form method='post' action='graph_neighbours.php' enctype='multipart/form-data'>
  &nbsp;&nbsp;&nbsp;<a href = 'help.graph_neighbours.html#formats'><B>Input format</B></a>");
  if (!$pipe) {
    echo ("
    &nbsp;<select name='in_format'>
    <option selected value = 'tab'> tab-delimited format
    <option value = 'adj_matrix'> Adjacency matrix
    <option value = 'gml'> GML format
    </select><br><br>");
  } else {
    
    echo "<input type='hidden' NAME='in_format' VALUE='$graph_format'>";
  }
  echo("&nbsp;&nbsp;&nbsp;<a href = 'help.graph_neighbours.html#out_format'><B>Output format</B></a>");
  echo("&nbsp;<select name='stats'>
    <option selected value = '0'> Classification output
    <option value = '1'> One line for each seed node (weighted graph only)</select><br><br>
    ");
  echo("&nbsp;&nbsp;&nbsp;<B><a href = 'help.display_graph.html#layout'>Distance from seed node</a></B><input type='text' name='steps' value = '$default_step' size = 1/><br>");

  if (!$pipe) {
    if ($demo) {
    demo("This demonstration graph consists in the yeast co-immunopreciptation interaction dataset described in <a href = 'http://www.ncbi.nlm.nih.gov/sites/entrez?Db=pubmed&Cmd=ShowDetailView&TermToSearch=16429126&ordinalpos=1&itool=EntrezSystem2.PEntrez.Pubmed.Pubmed_ResultsPanel.Pubmed_RVDocSum'>Gavin et al (2006)</a>. It contains 1430 nodes and 6531 edges. We will look for the neighbours of some polypeptides to see if they are linked to other proteins of similar functions.");    }
    echo ("<b>Graph</b><br>");
    echo ("<textarea name='graph' rows='6' cols='65'>$demo_graph</textarea>
    <br>Upload graph from file : <br>
    <input type='file' name='graph_file' size='45' /><br>
    &nbsp;&nbsp;&nbsp;
    <br><a href = 'help.graph_neighbours.html#columns'>Column specification (only relevant for tab-delimited input)</a><br>
    <table>
    <tr><td><B><a href = 'help.graph_neighbours.html#scol'>Source node</a></B></td><td><input type = 'text' name='s_col' value = '$default_scol' size = 1></input></td></tr>
    <tr><td><B><a href = 'help.graph_neighbours.html#scol'>Target node</a></B></td><td><input type = 'text' name='t_col' value = '$default_tcol' size = 1></input></td></tr>
    <tr><td><B><a href = 'help.graph_neighbours.html#wcol'>Weight or label column</a></B></td><td><input type = 'text' name='w_col' size = 1 value = '$default_wcol'></input></td></tr>
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
    echo("<a href = 'help.graph_neighbours.html#direction'><B>Direction of the neighbourhood</B></a>");
    echo ("
    &nbsp;<select name='direction'>
    <option selected value = 'all'> All directions
    <option value = 'out'> Out neighbours
    <option value = 'in'> In neighbours
    </select><br><br>");

  echo("<br><br><input type='checkbox' name='self' value='on' checked/>&nbsp;<B><a href = 'help.display_graph.html#layout'>Include each node in its neighborhood (with a distance of zero)</a></B><br>");

  echo("<hr>
  <br><br><b><a href = 'help.graph_neighbours.html#seeds'>Seed node selection</a></b><br>
  <br><table>
  <tr><td><input type = 'radio' $all_nodes_selected name='allseeds' value = 'all'/></td><td>All nodes</td>
  <tr><td><input type='radio' $selection_nodes_selected name='allseeds' value='selection' /></select></td><td>Node selection</td><tr><td></td><td>
  <textarea name='seeds' rows='6' cols='65'>$demo_nodes</textarea>
  <br>Upload seed nodes from file : <br>
  <input type='file' name='seeds_file' size='45' /></td></tr></table>");
  echo ("
  <ul><ul><table class='formbutton'>
  <TD><input type='submit' name='.submit' value='GO' /></TD>
  <TD><B><A HREF='graph_neighbours_form.php?demo=0'>RESET</A></B></TD>
  <TD><B><A HREF='graph_neighbours_form.php?demo=1'>DEMO</A></B></TD>
  </form>

  <TD><B><A HREF='help.graph_neighbours.html'>MANUAL</A></B></TD>
  <TD><B><A target = '_blank' HREF='".checkNeatTutorial("tutorials/neat_tutorial/Study_neighborhood_nodes.html")."'>TUTORIAL</A></B></TD>

  <TD><B><A HREF='mailto:sylvain@bigre.ulb.ac.be'>MAIL</A></B></TD>
  </TR></TABLE></ul></ul>");


?>
</body>
</html>
