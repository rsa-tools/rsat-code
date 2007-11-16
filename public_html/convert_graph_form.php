<html>
<head>
   <title>GrA-tools - convert-graph</title>
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
    $demo_graph = $uetz;
  }

  title('convert-graph');
  echo ("<center>Convert graphs between different formats.</center>\n");
  echo ("<form method='post' action='convert_graph.php' enctype='multipart/form-data'>
  &nbsp;&nbsp;&nbsp;<a href = 'help.convert_graph.html#formats'><B>Input format</B></a>&nbsp;<select name='in_format'>
  <option selected value = 'tab'> tab-delimited format
  <option value = 'adj_matrix'> Adjacency matrix
  <option value = 'gml'> GML format
  </select><br><br>
  &nbsp;&nbsp;&nbsp;<B><a href = 'help.convert_graph.html#formats'>Output format</a></B>&nbsp;<select name='out_format'>
  <option value = 'tab'> tab-delimited format
  <option selected value = 'gml'> GML format
  <option value = 'adj_matrix'> Adjacency matrix
  </select><br>");
  if ($demo) {
    demo("This demonstration graph is the yeast two-hybrid dataset produced by <a target = 'top' href = 'http://www.ncbi.nlm.nih.gov/sites/entrez?db=pubmed&uid=10688190&cmd=showdetailview&indexed=google'>Uetz et al (2001)</a>. It consists in 865 interactions between 926 proteins.");
  }
  echo("<br>
  <textarea name='graph' rows='6' cols='65'>$demo_graph</textarea>
  <br>Upload graph from file : <br>
  <input type='file' name='graph_file' size='45' /><br>
  &nbsp;&nbsp;&nbsp;
  
  <br><a href = 'help.convert_graph.html#columns'>Column specification (only relevant for tab-delimited input)</a><br>
  <table>
  <tr><td><B><a href = 'help.convert_graph.html#scol'>Source node</a></B></td><td><input type = 'text' name='s_col' value = '$default_scol' size = 1></input></td></tr>
  <tr><td><B><a href = 'help.convert_graph.html#scol'>Target node</a></B></td><td><input type = 'text' name='t_col' value = '$default_tcol' size = 1></input></td></tr>
  <tr><td><B><a href = 'help.convert_graph.html#wcol'>Weight or label column</a></B></td><td><input type = 'text' name='w_col' size = 1></input></td></tr>
  <tr><td><B><a href = 'help.convert_graph.html#eccol'>Edge color column</a></B></td><td><input type = 'text' name='ec_col' size = 1></input></td></tr>
  <tr><td><B><a href = 'help.convert_graph.html#sccol'>Source node color column</a></B></td><td><input type = 'text' name='sc_col' size = 1></input></td></tr>
  <tr><td><B><a href = 'help.convert_graph.html#sccol'>Target node color column</a></B></td><td><input type = 'text' name='tc_col' size = 1></input></td></tr>
  </table>
  <input type='checkbox' name='layout' value='on' $default_layout/>&nbsp;<B><a href = 'help.convert_graph.html#layout'>Calculate the layout of the nodes (only relevant for GML output)</a></B><br>
  <input type='checkbox' name='undirected' value='on' />&nbsp;<B><a href = 'help.convert_graph.html#undirected'>Undirected graph (only relevant for adjacency matrix input and output)</a></B><br>
  <ul><ul><table class='formbutton'>
  <TD><input type='submit' name='.submit' value='GO' /></TD>
  <TD><B><A HREF='convert_graph_form.php?demo=0'>RESET</A></B></TD>
  <TD><B><A HREF='convert_graph_form.php?demo=1'>DEMO</A></B></TD>
  </form>
  <TD><B><A HREF='help.convert_graph.html'>MANUAL</A></B></TD>
  <TD><B><A HREF='mailto:sylvain@scmbb.ulb.ac.be'>MAIL</A></B></TD>
  </TR></TABLE></ul></ul>
 ");

?>
</body>
</html>