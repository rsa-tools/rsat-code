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
$default_decimals = 2;
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
  $stat = $_REQUEST['stat'];
  $decimals = $_REQUEST['decimals'];

  title('graph-cluster-membership');
  echo ("<center>Map a clustering result onto a graph, and compute the membership degree
    between each node and each cluster, on the basis of egdes linking this
    node to the cluster.\n");
   echo ("<br>This program was developed by 
	  <a target=_blank href=http://www.bigre.ulb.ac.be/people/Members/gipsi/>Gipsi Lima-Mendez</a>
	  With the help of:
	  <a target=_blank href=http://www.bigre.ulb.ac.be/Users/sylvain/>Sylvain Broh&eacute;e</a> and
          <a target=_blank href=http://www.bigre.ulb.ac.be/Users/jvanheld/>Jacques van Helden</a>.
          </center>\n");
  echo ("<form method='post' action='graph_cluster_membership.php' enctype='multipart/form-data'>
  Graph<br>
  &nbsp;&nbsp;&nbsp;<a href = 'help.graph_cluster_membership.html#formats'><B>Input format</B></a>&nbsp;");
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
      demo("This demonstration consists on mapping of the clusters produced by the <a href = 'http://micans.org/mcl/'>MCL</a> algorithm (inflation 2.1) of the yeast co-immunoprecipitation interaction dataset from <a href = 'http://www.ncbi.nlm.nih.gov/sites/entrez?Db=pubmed&Cmd=ShowDetailView&TermToSearch=16429126&ordinalpos=1&itool=EntrezSystem2.PEntrez.Pubmed.Pubmed_ResultsPanel.Pubmed_RVDocSum' target = 'top'>Gavin (2006)</a> onto the corresponding graph");
    }
    echo("
      <b>Graph</b><br>
      <textarea name='graph' rows='6' cols='65'>$demo_graph</textarea>
      <br>Upload graph from file : <br>
      <input type='file' name='graph_file' size='45' /><br>
      &nbsp;&nbsp;&nbsp;
      <br><a href = 'help.graph_cluster_membership.html#columns'>Column specification (only relevant for tab-delimited input)</a><br>
      <table>
      <tr><td><B><a href = 'help.graph_cluster_membership.html#scol'>Source node</a></B></td><td><input type = 'text' name='s_col' value = '$default_scol' size = 1></input></td></tr>
      <tr><td><B><a href = 'help.graph_cluster_membership.html#scol'>Target node</a></B></td><td><input type = 'text' name='t_col' value = '$default_tcol' size = 1></input></td></tr>
      <tr><td><B><a href = 'help.graph_cluster_membership.html#wcol'>Weight column</a></B></td><td><input type = 'text' name='w_col' size = 1></input></td></tr>
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
  </select><br><br><b>Clusters</b><br>");
  if ($cluster_file == "") {
    echo("
    <textarea name='clusters' rows='6' cols='65'>$demo_clusters</textarea>
    <br>Upload clusters from file : <br>
    <input type='file' name='clusters_file' size='45' /><br>");
  } else {
    echo("
    <input type='hidden' name='pipe_clusters_file' value = '$cluster_file' /><br>");
    info_link("Clusters uploaded from the previous treatment", rsat_path_to_url($cluster_file));    
  }
  echo ("<br><br><form method='post' action='graph_cluster_membership.php' enctype='multipart/form-data'>
  <b>Membership Matrix<br>
  &nbsp;&nbsp;&nbsp;<a href = 'help.graph_cluster_membership.html#stats'><B>Stat</B></a>&nbsp;");
  if ($demo) {
    echo ("  
    <select name='stat'>
    <option selected value = 'edge'> edge
    <option value = 'reledge'> relative edge
    <option value = 'weight'> weight
    <option value = 'relw'> relative weight
    </select> <br><br>");
  }
  else {
    echo ("  
    <select name='stat'>
    <option selected value = 'weight'> weight
    <option value = 'relw'> relative weight
    <option value = 'edge'> edge
<option value = 'reledge'> relative edge
    </select> <br><br>");
 }
  echo ("<table>
      <tr><td><B><a href = 'help.graph_cluster_membership.html#decimals'>Number of Decimals</a></B></td><td><input type = 'text' name='decimals' value = '$default_decimals' size = 1></input></td></tr>
");
 
  echo ("
  <ul><ul><table class='formbutton'>
  <TD><input type='submit' name='.submit' value='GO' /></TD>
  <TD><B><A HREF='graph_cluster_membership_form.php?demo=0'>RESET</A></B></TD>
  <TD><B><A HREF='graph_cluster_membership_form.php?demo=1'>DEMO</A></B></TD>
  </form>
  <TD><B><A HREF='help.graph_cluster_membership.html'>MANUAL</A></B></TD>
  <TD><B><A HREF='mailto:gipsi@bigre.ulb.ac.be'>MAIL</A></B></TD>
  </TR></TABLE></ul></ul>");


?>
</body>
</html>
