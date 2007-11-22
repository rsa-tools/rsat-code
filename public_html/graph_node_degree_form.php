<html>
<head>
   <title>GrA-tools - graph-node-degree</title>
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
  # demo graph
  $demo = $_REQUEST['demo'];
  if ($demo == 1) {
    $demo_graph = $uetz;
    $demo_nodes = $uetz_nodes;
  }

  title('graph-node-degree');
  echo ("<center>Calculate the node degree of each node and specifies if this node is a
    seed or a target node.</center>\n");
  echo ("<form method='post' action='graph_node_degree.php' enctype='multipart/form-data'>
  &nbsp;&nbsp;&nbsp;<a href = 'help.graph_node_degree.html#formats'><B>Input format</B></a>&nbsp;<select name='in_format'>
  <option selected value = 'tab'> tab-delimited format
  <option value = 'adj_matrix'> Adjacency matrix
  <option value = 'gml'> GML format
  </select><br><br><b>Graph</b><br>");
  if ($demo) {
    demo("This demonstration graph is the yeast two-hybrid dataset produced by <a target = 'top' href = 'http://www.ncbi.nlm.nih.gov/sites/entrez?db=pubmed&uid=10688190&cmd=showdetailview&indexed=google'>Uetz et al (2001)</a>. It consists in 865 interactions between 926 proteins. The nodes for which the degree is seeked were chosen in this network amongst the more or the less connected nodes." );
  }
  echo ("<textarea name='graph' rows='6' cols='65'>$demo_graph</textarea>
  <br>Upload graph from file : <br>
  <input type='file' name='graph_file' size='45' /><br>
  &nbsp;&nbsp;&nbsp;
  
  <br><a href = 'help.graph_node_degree.html#columns'>Column specification (only relevant for tab-delimited input)</a><br>
  <table>
  <tr><td><B><a href = 'help.graph_node_degree.html#scol'>Source node</a></B></td><td><input type = 'text' name='s_col' value = '$default_scol' size = 1></input></td></tr>
  <tr><td><B><a href = 'help.graph_node_degree.html#scol'>Target node</a></B></td><td><input type = 'text' name='t_col' value = '$default_tcol' size = 1></input></td></tr>
  </table>
  <br>
  </select><br><br><b><a href = 'help.graph_node_degree.html#nodes'>Nodes</a></b><br>
  <textarea name='nodes' rows='6' cols='65'>$demo_nodes</textarea>
  <br>Upload nodes from file : <br>
  <input type='file' name='nodes_file' size='45' /><br>");
  echo ("
  <ul><ul><table class='formbutton'>
  <TD><input type='submit' name='.submit' value='GO' /></TD>
  <TD><B><A HREF='graph_node_degree_form.php?demo=0'>RESET</A></B></TD>
  <TD><B><A HREF='graph_node_degree_form.php?demo=1'>DEMO</A></B></TD>
  </form>
  <TD><B><A HREF='help.graph_node_degree.html'>MANUAL</A></B></TD>
  <TD><B><A HREF='mailto:sylvain@scmbb.ulb.ac.be'>MAIL</A></B></TD>
  </TR></TABLE></ul></ul>");


?>
</body>
</html>