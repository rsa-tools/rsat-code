<html>
<head>
   <title>Network Analysis Tools - random-graph</title>
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
  $default_random_type = "selected";
  # demo graph
  $demo = $_REQUEST['demo'];
  if ($demo == 1) {
    $demo_graph = storeFile("demo_files/protein_interactions_uetz.tab");
    $demo_random_type = "selected";
    $default_random_type = "";
  }
  # PIPE VALUES
  $pipe = $_REQUEST['pipe'];
  $graph_file = $_REQUEST['graph_file'];
  $graph_format = $_REQUEST['graph_format'];
  $scol = $_REQUEST['scol'];
  $tcol = $_REQUEST['tcol'];
  $wcol = $_REQUEST['wcol'];
  $eccol = $_REQUEST['eccol'];
  
  title('random-graph');
  echo ("<center>Generate random graphs either from an existing graph or from scratch 
    according to different randomization procedures.\n");
   echo ("<br>This program was developed by 
	  <a target=_blank href=http://www.bigre.ulb.ac.be/Users/sylvain/>Sylvain Broh&eacute;e</a> and
          <a target=_blank href=http://www.bigre.ulb.ac.be/Users/jvanheld/>Jacques van Helden</a>.
          </center>\n");   
  echo ("<form method='post' action='random_graph.php' enctype='multipart/form-data'>");
  echo ("&nbsp;&nbsp;&nbsp;<B><a href = 'help.random_graph.html#formats'>Input format</a></B>&nbsp;");

  if (!$pipe) {
  echo("
    &nbsp;<select name='in_format'>
    <option selected value = 'tab'> tab-delimited format
    <option value = 'adj_matrix'> Adjacency matrix
    <option value = 'gml'> GML format
    </select><br><br>");
  } else {
    echo ": $graph_format<br>";
    echo "<input type='hidden' NAME='in_format' VALUE='$graph_format'>";
  }
  echo ("&nbsp;&nbsp;&nbsp;<B><a href = 'help.random_graph.html#formats'>Output format</a></B>&nbsp;<select name='out_format'>
  <option selected value = 'tab'> tab-delimited format
  <option value = 'gml'> GML format
  <option value = 'adj_matrix'> Adjacency matrix
  </select><br><br>");
  if (!$pipe) {
    if ($demo) {
      demo("This demonstration graph is the yeast two-hybrid dataset produced by <a target = 'top' href = 'http://www.ncbi.nlm.nih.gov/sites/entrez?db=pubmed&uid=10688190&cmd=showdetailview&indexed=google'>Uetz et al (2001)</a>. It consists in 865 interactions between 926 proteins.");
    }
    echo ("<b>Graph</b><br>
    <textarea name='graph' rows='6' cols='65'>$demo_graph</textarea>
    <br>Upload graph from file : <br>
    <input type='file' name='graph_file' size='45' /><br>
    &nbsp;&nbsp;&nbsp;
    <br><a href = 'help.random_graph.html#columns'>Column specification (only relevant for tab-delimited input)</a><br>
    <table>
    <tr><td><B><a href = 'help.random_graph.html#scol'>Source node</a></B></td><td><input type = 'text' name='s_col' value = '$default_scol' size = 1></input></td></tr>
    <tr><td><B><a href = 'help.random_graph.html#scol'>Target node</a></B></td><td><input type = 'text' name='t_col' value = '$default_tcol' size = 1></input></td></tr>
    <tr><td><B><a href = 'help.random_graph.html#wcol'>Weight or label column</a></B></td><td><input type = 'text' name='w_col' size = 1></input></td></tr>
    </table>");
  } else {
    info_link("Graph uploaded from the previous treatment", rsat_path_to_url($graph_file));
    echo "<input type='hidden' NAME='pipe_graph_file' VALUE='$graph_file'>";
  }
  echo ("<br>
  <input type='checkbox' name='directed' value='on' />&nbsp;<B><a href = 'help.random_graph.html#directed'>Directed graph</a></B><br>
  <table><tr>
  <B><td><a href = 'help.random_graph.html#random_type'>Randomization type</a></B>&nbsp;<select name='random_type'>
  <option $default_random_type value = 'scratch'> Erdos-Renyi from scratch
  <option $demo_random_type value = 'ER'> Erdos-Renyi randomization (ER)
  <option value = 'node_degree'> Node degree conservation
  <option value = 'node_degree_distrib'> Node degree distribution conservation 
  </select></td>
  <td>number of nodes</td><td><input type = 'text' name='nodes' size = 4></td> 
  <td>number of edges</td><td><input type = 'text' name='edges' size = 4></td>
  <td><i>(node and edge numbers must be provided only for randomization type 'from scratch')</i></td>
  </tr><tr>
  <td><input type='checkbox' name='normal' value='on' />&nbsp;<B><a href = 'help.random_graph.html#normal'>Normal distribution of weights</a></B></td>
  <td>Mean</td><td><input type = 'text' name='mean' size = 4></input></td>
  <td>Standard deviation</td><td><input type = 'text' name='stddev' size = 4></input></td>
  <td><i>(mean and standard deviation may be provided only for
  randomization types 'ER' and for 'from scratch')</i></td>
  </tr></table>
  <br>
  <input type='checkbox' name='duplicate' value='on' />&nbsp;<B><a href = 'help.random_graph.html#duplicate'>Allow duplicated edges</a></B><br>
  <input type='checkbox' name='self' value='on' />&nbsp;<B><a href = 'help.random_graph.html#self'>Allow self loops</a></B><br>
  <input type='checkbox' name='col_conservation' value='on' />&nbsp;<B><a href = 'help.random_graph.html#col_conservation'>Source and target nodes stay source and target nodes in the randomized graph (only for ER randomization type)</a></B><br>
  <input type='checkbox' name='no_single' value='on' />&nbsp;<B><a href = 'help.random_graph.html#no_single'>Prevent the graph from containing nodes with no neighbour (only for ER randomization type)</a></B><br>
  <ul><ul><table class='formbutton'>
  <TD><input type='submit' name='.submit' value='GO' /></TD>
  <TD><B><A HREF='random_graph_form.php?demo=0'>RESET</A></B></TD>
  <TD><B><A HREF='random_graph_form.php?demo=1'>DEMO</A></B></TD>
  </form>
  <TD><B><A HREF='help.random_graph.html'>MANUAL</A></B></TD>
  <TD><B><A target = '_top' HREF='".checkNeatTutorial("tutorials/neat_tutorial/Influence_graph_alteration.html")."'>TUTORIAL</A></B></TD>

  <TD><B><A HREF='mailto:sylvain@bigre.ulb.ac.be'>MAIL</A></B></TD>
  </TR></TABLE></ul></ul>
 ");

?>
</body>
</html>
