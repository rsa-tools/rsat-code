<html>
<head>
   <title>Network Analysis Tools - graph-get-clusters</title>
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
  # demo graph
  $demo = $_REQUEST['demo'];
  if ($demo == 1) {
    $demo_graph = storeFile("demo_files/protein_interactions_gavin_2006.tab");
    $demo_clusters = storeFile("demo_files/gavin_mcl_clusters_inf2.1.tab");
  }
  # PIPE VALUES
  $pipe = $_REQUEST['pipe'];
  $graph_file = $_REQUEST['graph_file'];
  $graph_format = $_REQUEST['graph_format'];
  $scol = $_REQUEST['scol'];
  $tcol = $_REQUEST['tcol'];
  $wcol = $_REQUEST['wcol'];
  $cluster_file = $_REQUEST['cluster_file'];
  title('graph-get-clusters');
  echo ("<center>Compares a graph with a classification/clustering file.\n");
   echo ("<br>This program was developed by 
	  <a target=_blank href=http://www.bigre.ulb.ac.be/Users/sylvain/>Sylvain Broh&eacute;e</a> and
          <a target=_blank href=http://www.bigre.ulb.ac.be/Users/jvanheld/>Jacques van Helden</a>.
          </center>\n");
  echo ("<form method='post' action='graph_get_clusters.php' enctype='multipart/form-data'>
  Graph<br>
  &nbsp;&nbsp;&nbsp;<a href = 'help.graph_get_clusters.html#formats'><B>Input format</B></a>&nbsp;");
  if (!$pipe) {
    echo ("  
    <select name='in_format'>
    <option selected value = 'tab'> tab-delimited format
    <option value = 'adj_matrix'> Adjacency matrix
    <option value = 'gml'> GML format
    </select> <br><br>");
  } else {
    echo ": $graph_format<br>";
    echo "<input type='hidden' NAME='in_format' VALUE='$graph_format'>";
  }
  if (!$pipe) {
    if ($demo) {
      demo("This demonstration consists in the comparaison of the yeast co-immunoprecipitation interaction dataset from <a href = 'http://www.ncbi.nlm.nih.gov/sites/entrez?Db=pubmed&Cmd=ShowDetailView&TermToSearch=16429126&ordinalpos=1&itool=EntrezSystem2.PEntrez.Pubmed.Pubmed_ResultsPanel.Pubmed_RVDocSum' target = 'top'>Gavin (2006)</a> with the clusters produced by the <a href = 'http://micans.org/mcl/'>MCL</a> algorithm (inflation 2.1)");
    }
    echo("
      <b>Graph</b><br>
      <textarea name='graph' rows='6' cols='65'>$demo_graph</textarea>
      <br>Upload graph from file : <br>
      <input type='file' name='graph_file' size='45' /><br>
      &nbsp;&nbsp;&nbsp;
      <br><a href = 'help.graph_get_clusters.html#columns'>Column specification (only relevant for tab-delimited input)</a><br>
      <table>
      <tr><td><B><a href = 'help.graph_get_clusters.html#scol'>Source node</a></B></td><td><input type = 'text' name='s_col' value = '$default_scol' size = 1></input></td></tr>
      <tr><td><B><a href = 'help.graph_get_clusters.html#scol'>Target node</a></B></td><td><input type = 'text' name='t_col' value = '$default_tcol' size = 1></input></td></tr>
      <tr><td><B><a href = 'help.graph_get_clusters.html#wcol'>Weight or label column</a></B></td><td><input type = 'text' name='w_col' size = 1></input></td></tr>
      </table>"
    );
  } else {
    info_link("Graph uploaded from the previous treatment", rsat_path_to_url($graph_file));    
    echo "<input type='hidden' NAME='pipe_graph_file' VALUE='$graph_file'>";
    echo "<input type='hidden' NAME='s_col' VALUE='$scol'>";
    echo "<input type='hidden' NAME='t_col' VALUE='$tcol'>";
    echo "<input type='hidden' NAME='w_col' VALUE='$wcol'>";
  }
  echo("<br>
  </select><br><br><b>Clusters (or list of nodes in case of induced graph)</b><br>");
  if ($cluster_file == "") {
    echo("
    <textarea name='clusters' rows='6' cols='65'>$demo_clusters</textarea>
    <br>Upload clusters (or nodes for induction) from file : <br>
    <input type='file' name='clusters_file' size='45' /><br>");
  } else {
    echo("
    <input type='hidden' name='pipe_clusters_file' value = '$cluster_file' /><br>");
    info_link("Clusters (or nodes) uploaded from the previous treatment", rsat_path_to_url($cluster_file));    
  }
  
  echo("
  <B><a href = 'help.compare_graphs.html#return'>Output type</B></a>&nbsp;<select name='return'>
  <option value = 'table'> node-cluster connections
  <option selected value = 'clusters'> intra-cluster edges
  <option value = 'graph'> annotated graph (all edges)
  </select><br><br>");
  echo("
    <B><a href = 'help.graph_get_clusters.html#formats'>Output format (not relevant for node-cluster connection output format)</a></B>&nbsp;<select name='out_format'>
    <option selected value = 'tab'> tab-delimited format
    <option value = 'gml'> GML format
    <option value = 'adj_matrix'> Adjacency matrix
    </select><br><br>");
 
    

  ## This option is useful to duplicate nodes belonging to more than one cluster but for the moment
  ## this option does not work anymore (-distinct)
  #echo("<input type='checkbox' name='distinct' value='on' />&nbsp;<B><a href = 'help.convert_graph.html#distinct'>Duplicate the nodes belonging to more than one cluster</a></B><br>");
  echo ("
  <ul><ul><table class='formbutton'>
  <TD><input type='submit' name='.submit' value='GO' /></TD>
  <TD><B><A HREF='graph_get_clusters_form.php?demo=0'>RESET</A></B></TD>
  <TD><B><A HREF='graph_get_clusters_form.php?demo=1'>DEMO</A></B></TD>
  </form>
  <TD><B><A HREF='help.graph_get_clusters.html'>MANUAL</A></B></TD>
  <TD><B><A HREF='mailto:sylvain@bigre.ulb.ac.be'>MAIL</A></B></TD>
  </TR></TABLE></ul></ul>");


?>
</body>
</html>
