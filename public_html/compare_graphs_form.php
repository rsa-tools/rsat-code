<html>
<head>
   <title>Network Analysis Tools - compare-graphs</title>
   <link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
</head>
<body class="form">
<?php
  require ('functions.php');
//   require ('demo_dataset.php');
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
    $demo_graphQ = storeFile("demo_files/protein_interactions_uetz.tab");
    $demo_graphR = storeFile("demo_files/protein_interactions_ito.tab");

    $demo_remark = "In this demonstration, we will compare the networks resulting from the
two first publications reporting a complete characterization of the
yeast interactome, obtained using the two-hybrid method.";
    $demo_remark .= "The first network (<a target = '_top' href = 'http://www.ncbi.nlm.nih.gov/sites/entrez?db=pubmed&uid=10688190&cmd=showdetailview&indexed=google'>Uetz et al, 2000)</a> contains 865 interactions between 926 proteins.";
    $demo_remark .= "The second network  (<a href = 'http://www.ncbi.nlm.nih.gov/sites/entrez?Db=pubmed&Cmd=ShowDetailView&TermToSearch=11283351&ordinalpos=3&itool=EntrezSystem2.PEntrez.Pubmed.Pubmed_ResultsPanel.Pubmed_RVDocSum' target = 'top'>Ito et al, 2001)</a> contains 786 interactions between 779 proteins.";

    $demo_remark .= "We will merge the two networks (i.e. compute
their union), and label each edge according to the fact that it is
found only in Ito's network, only in Uetz' network, or in both. We
will also compute the statistical significance of the intersection
between the two networks.";

  }
  title('compare-graphs');
echo ("<center>Computes the intersection, the union or the difference between two graphs.\n");
echo ("<br>This program was developed by 
<a target=_blank href=http://www.bigre.ulb.ac.be/Users/sylvain/>Sylvain Broh&eacute;e</a>,
Gilles Vanderstocken 
and 
<a target=_blank href=http://www.bigre.ulb.ac.be/Users/jvanheld/>Jacques van Helden</a>.
</center>\n");
  if($demo) {
    demo($demo_remark);
  }
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
    echo ("<textarea name='graphQ' rows='6' cols='40'>$demo_graphQ</textarea></td>");
  }
  echo "</td>";
  ## REFERENCE INPUT GRAPH TEXTAREA
  echo ("<td><textarea name='graphR' rows='6' cols='40'>$demo_graphR</textarea></td>");
  echo ("<tr><td>");
  ## QUERY INPUT GRAPH FILE
  if (!$pipe) {
    echo ("Upload query graph from file : <br>
    <input type='file' name='graph_fileQ' size='40' />");
  } else {
    info_link("Graph uploaded from the previous treatment", rsat_path_to_url($query_graph_file));
    echo "<input type='hidden' NAME='pipe_query_graph_file' VALUE='$query_graph_file'/>";
  }
  echo "</td>";
  ## REFERENCE INPUT GRAPH FILE
  echo ("<td>Upload reference graph from file : <br>
  <input type='file' name='graph_fileR' size='40' /></td>");
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
  <B><a href = 'help.compare_graphs.html#return'>Output type</B></a>&nbsp;<select name='return'>
  <option selected value = 'union'> union
  <option value = 'intersection'> intersection
  <option value = 'difference'> difference
  <option value = 'intersection+Q'> intersection + query
  <option value = 'intersection+R'> intersection + reference
  </select><br><br>");
  echo ("<ul>
  <B><a href = 'help.compare_graphs.html#formats'>Output format</a></B>&nbsp;<select name='out_format'>
  <option selected value = 'tab'> tab-delimited format
  <option value = 'gml'> GML format
  <option value = 'adj_matrix'> Adjacency matrix
  <option value = 'dot'> Dot format
  </select><br><br>");
  ## Weight/label on the edges
  echo ("
  <B><a href = 'help.compare_graphs.html#outweights'>Weight/label on the edges of the output graph</a></B></a>&nbsp;<select name='outweight'>
  <option value = 'Q'> weight / label of the query
  <option value = 'R'> weight / label of the reference
  <option value = 'sum'> sum of the weights
  <option value = 'mean'> arithmetic mean of the weights
  <option value = 'mean.g'> geometrical mean of the weights
  <option value = 'min'> min of the weights
  <option value = 'max'> max of the weights
  <option selected value = 'Q::R'> weights / label of query and reference
  </select><br><br>");
  ## Directed graph
  echo("
    <input type='checkbox' name='directed' value='on'/>&nbsp;<B><a href = 'help.compare_graphs.html#directed'>Graphs must be considered as directed</a></B><br><br>");
  echo("<input type='checkbox' name='self' value='on'/>&nbsp;<B><a href = 'help.compare_graphs.html#directed'>Graph admits self loops</a></B><br><br></ul>");
  echo ("<ul><ul><table class='formbutton'>
  <TD><input type='submit' name='.submit' value='GO' /></TD>
  <TD><B><A HREF='compare_graphs_form.php?demo=0'>RESET</A></B></TD>
  <TD><B><A HREF='compare_graphs_form.php?demo=1'>DEMO</A></B></TD>
  </form>
  <TD><B><A HREF='help.compare_graphs.html'>MANUAL</A></B></TD>
  <TD><B><A target = '_blank' HREF='".checkNeatTutorial("tutorials/neat_tutorial/Comparisons_between_network.html")."'>TUTORIAL</A></B></TD>

  <TD><B><A HREF='mailto:sylvain@bigre.ulb.ac.be' >MAIL</A></B></TD>
  </TR></TABLE></ul></ul>");
  ?>
  </body>
  </html>
