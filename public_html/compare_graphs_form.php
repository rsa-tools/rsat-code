<html>
<head>
   <title>GrA-tools - compare-graphs</title>
   <link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
</head>
<body class="form">
<?
  require ('functions.php');
  require ('demo_dataset.php');
  # variables definition
  $default_scolQ = 1;
  $default_tcolQ = 2;
  $default_scolR = 1;
  $default_tcolR = 2;
  # PIPE VALUES
  $pipe = $_REQUEST['pipe'];
  $query_graph_file = $_REQUEST['graph_file'];
  $query_graph_format = $_REQUEST['graph_format'];
  
  $query_scol = $_REQUEST['scol'];
  $query_tcol = $_REQUEST['tcol'];
  $query_wcol = $_REQUEST['wcol'];
  
  
  
  # demo graph
  $demo = $_REQUEST['demo'];
  if ($demo == 1) {
    $demo_graphQ = $uetz;
    $demo_graphR = $ito;
  }
  title('compare-graphs');
  echo ("<center>Computes the intersection, the union or the difference of two graphs</center>\n");
  echo ("<form method='post' action='compare_graphs.php' enctype='multipart/form-data'>");
  echo ("<table>\n");
  echo ("<tr><th>Query graph</th><th>Reference graph</th></tr>\n");
  ## QUERY INPUT FORMAT
  echo "<tr><td>";
  echo "<B>Input format</B></a>&nbsp;";
  if (!$pipe) {
    echo ("
    <select name='in_formatQ'>
    <option selected value = 'tab'> tab-delimited format
    <option value = 'adj_matrix'> Adjacency matrix
    <option value = 'gml'> GML format
    </select><br>");
  } else {
    echo ": $query_graph_format<br>";
    echo "<input type='hidden' NAME='in_formatQ' VALUE='$query_graph_format'>";
  } 
  echo "</td>";
  ## REFERENCE INPUT FORMAT
  echo ("<td>
  <B>Input format</B></a>&nbsp;<select name='in_formatR'>
  <option selected value = 'tab'> tab-delimited format
  <option value = 'adj_matrix'> Adjacency matrix
  <option value = 'gml'> GML format
  </select><br></td></tr>");
  ## QUERY INPUT GRAPH TEXTAREA
  echo "<td>";
  if (!$pipe) {
    echo ("<textarea name='graphQ' rows='6' cols='65'>$demo_graphQ</textarea></td>");
  }
  echo "</td>";
  ## REFERENCE INPUT GRAPH TEXTAREA
  echo ("<td><textarea name='graphR' rows='6' cols='65'>$demo_graphR</textarea></td>");
  echo ("<tr><td>");
  ## QUERY INPUT GRAPH FILE
  if (!$pipe) {
    echo ("Upload query graph from file : <br>
    <input type='hidden' name='MAX_FILE_SIZE' value='30000' />
    <input type='file' name='graph_fileQ' size='45' />");
  } else {
    info("Query graph uploaded from the previous treatment");
    echo "<input type='hidden' NAME='pipe_query_graph_file' VALUE='$query_graph_file'/>";
  }
  echo "</td>";
  ## REFERENCE INPUT GRAPH FILE
  echo ("<td>Upload reference graph from file : <br>
  <input type='hidden' name='MAX_FILE_SIZE' value='30000' />
  <input type='file' name='graph_fileR' size='45' /></td>");
  echo ("</tr>");
  ## QUERY INPUT GRAPH PARAMETERS
  echo ("<tr>");
  echo ("<td>");
  if (!$pipe) {
    echo("Column specification (only relevant for tab-delimited input)</a><br>
    <table>
    <tr><td><B><a href = 'help.compare_graphs.html#scol'>Source node</a></B></td><td><input type = 'text' name='s_colQ' value = '$default_scolQ' size = 1></input></td></tr>
    <tr><td><B><a href = 'help.compare_graphs.html#scol'>Target node</a></B></td><td><input type = 'text' name='t_colQ' value = '$default_tcolQ' size = 1></input></td></tr>
    <tr><td><B><a href = 'help.compare_graphs.html#wcol'>Weight or label column</a></B></td><td><input type = 'text' name='w_colQ' size = 1></input></td></tr>
    </table>");
  } else {
    if ($graph_format == 'tab') {
      echo "<input type='hidden' NAME='s_colQ' VALUE='$query_scol'>";
      echo "<input type='hidden' NAME='t_colQ' VALUE='$query_tcol'>";
      echo "<input type='hidden' NAME='w_colQ' VALUE='$query_wcol'>";
    }
  }
  echo ("</td");
  ## REFENCE INPUT GRAPH PARAMETERS
  echo ("<td>");
  echo("Column specification (only relevant for tab-delimited input)</a><br>
  <table>
  <tr><td><B><a href = 'help.compare_graphs.html#scol'>Source node</a></B></td><td><input type = 'text' name='s_colR' value = '$default_scolR' size = 1></input></td></tr>
  <tr><td><B><a href = 'help.compare_graphs.html#scol'>Target node</a></B></td><td><input type = 'text' name='t_colR' value = '$default_tcolR' size = 1></input></td></tr>
  <tr><td><B><a href = 'help.compare_graphs.html#wcol'>Weight or label column</a></B></td><td><input type = 'text' name='w_colR' size = 1></input></td></tr>
  </table>");
  echo ("</table><br><br>");
  ## OUTPUT SPECIFICATION

  echo("
  <B><a href = 'help.compare_graphs.html#return'>Output</B></a>&nbsp;<select name='return'>
  <option selected value = 'union'> union
  <option value = 'intersection'> intersection
  <option value = 'difference'> difference
  <option value = 'intersection+Q'> intersection + query
  <option value = 'intersection+R'> intersection + reference
  </select><br><br>");
  echo ("<ul>
  <B><a href = 'help.compare_graphs.html#formats'>Output format</a></B>&nbsp;<select name='out_format'>
  <option value = 'tab'> tab-delimited format
  <option selected value = 'gml'> GML format
  <option value = 'adj_matrix'> Adjacency matrix
  </select><br><br>");
  ## Weight on the edges
  echo ("
  <B><a href = 'help.compare_graphs.html#outweights'>Weight on the edges of the output graph</a></B></a>&nbsp;<select name='outweight'>
  <option value = 'Q'> weight (label) of the query
  <option value = 'R'> weight (label) of the reference
  <option value = 'sum'> sum of the weights
  <option value = 'mean'> arithemetic mean of the weights
  <option value = 'mean.g'> geometrical mean of the weights
  <option value = 'min'> min of the weights
  <option value = 'max'> max of the weights
  <option selected value = 'Q::R'> weights of query and reference
  </select><br><br>");
  ## Directed graph
  echo("
    <input type='checkbox' name='directed' value='on'/>&nbsp;<B><a href = 'help.compare_graphs.html#directed'>Graphs must be considered as directed</a></B><br><br>");
  echo("<input type='checkbox' name='self' value='on'/>&nbsp;<B><a href = 'help.compare_graphs.html#directed'>Graph admit self loops</a></B><br><br></ul>");
  echo ("<ul><ul><table class='formbutton'>
  <TD><input type='submit' name='.submit' value='GO' /></TD>
  <TD><B><A HREF='compare_graphs_form.php?demo=0'>RESET</A></B></TD>
  <TD><B><A HREF='compare_graphs_form.php?demo=1'>DEMO</A></B></TD>
  </form>
  <TD><B><A HREF='help.compare_graphs.html'>MANUAL</A></B></TD>
  <TD><B><A HREF='mailto:sylvain@scmbb.ulb.ac.be'>MAIL</A></B></TD>
  </TR></TABLE></ul></ul>");