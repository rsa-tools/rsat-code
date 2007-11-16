<html>
<head>
   <title>GrA-tools - graph-get-clusters</title>
   <link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
</head>
<body class="form">
<?
  require ('functions.php');
  require ('demo_dataset.php');
  # variable definition
  $default_scol = 1;
  $default_tcol = 2;
  $default_wcol = "";
  $default_layout = "checked";
  # demo graph
  $demo = $_REQUEST['demo'];
  if ($demo == 1) {
    $demo_graph = $gavin2006;
    $demo_clusters = $gavin_clusters_mcl_2_1;
  }

  title('graph-get-clusters');
  echo ("<center>Compares a graph with a classification/clustering file.</center>\n");
  echo ("<form method='post' action='graph_get_clusters.php' enctype='multipart/form-data'>
  &nbsp;&nbsp;&nbsp;<a href = 'help.graph_get_clusters.html#formats'><B>Input format</B></a>&nbsp;<select name='in_format'>
  <option selected value = 'tab'> tab-delimited format
  <option value = 'adj_matrix'> Adjacency matrix
  <option value = 'gml'> GML format
  </select><br><br><b>Graph</b><br>
  <textarea name='graph' rows='6' cols='65'>$demo_graph</textarea>
  <br>Upload graph from file : <br>
  <input type='file' name='graph_file' size='45' /><br>
  &nbsp;&nbsp;&nbsp;
  
  <br><a href = 'help.graph_get_clusters.html#columns'>Column specification (only relevant for tab-delimited input)</a><br>
  <table>
  <tr><td><B><a href = 'help.graph_get_clusters.html#scol'>Source node</a></B></td><td><input type = 'text' name='s_col' value = '$default_scol' size = 1></input></td></tr>
  <tr><td><B><a href = 'help.graph_get_clusters.html#scol'>Target node</a></B></td><td><input type = 'text' name='t_col' value = '$default_tcol' size = 1></input></td></tr>
  <tr><td><B><a href = 'help.graph_get_clusters.html#wcol'>Weight or label column</a></B></td><td><input type = 'text' name='w_col' size = 1></input></td></tr>
  </table>
  <br>
  </select><br><br><b>Clusters</b><br>
  <textarea name='clusters' rows='6' cols='65'>$demo_clusters</textarea>
  <br>Upload clusters from file : <br>
  <input type='file' name='clusters_file' size='45' /><br>
  <B><a href = 'help.graph_get_clusters.html#formats'>Output format (only useful for clusters output)</a></B>&nbsp;<select name='out_format'>
  <option value = 'tab'> tab-delimited format
  <option selected value = 'gml'> GML format
  <option value = 'adj_matrix'> Adjacency matrix
  </select><br><br>");
  echo("
  <B><a href = 'help.compare_graphs.html#return'>Output</B></a>&nbsp;<select name='return'>
  <option selected value = 'table'> contingency table
  <option selected value = 'clusters'> intra-cluster edges
  <option value = 'graph'> annotated graph (all edges)
  </select><br><br>");
  ## This option is useful to duplicate nodes belonging to more than one cluster but for the moment
  ## this option does not work anymore (-distinct)
  #echo("<input type='checkbox' name='distinct' value='on' />&nbsp;<B><a href = 'help.convert_graph.html#distinct'>Duplicate the nodes belonging to more than one cluster</a></B><br>");
  echo("<input type='checkbox' name='induced' value='on' />&nbsp;<B><a href = 'help.convert_graph.html#induced'>Induce the graph with the nodes of the cluster file</a></B><br>");
  echo ("
  <ul><ul><table class='formbutton'>
  <TD><input type='submit' name='.submit' value='GO' /></TD>
  <TD><B><A HREF='graph_get_clusters_form.php?demo=0'>RESET</A></B></TD>
  <TD><B><A HREF='graph_get_clusters_form.php?demo=1'>DEMO</A></B></TD>
  </form>
  <TD><B><A HREF='help.graph_get_clusters.html'>MANUAL</A></B></TD>
  <TD><B><A HREF='mailto:sylvain@scmbb.ulb.ac.be'>MAIL</A></B></TD>
  </TR></TABLE></ul></ul>");


?>
</body>
</html>