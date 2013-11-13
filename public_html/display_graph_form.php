<html>
<head>
   <title>Network Analysis Tools - display-graph</title>
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
  $default_layout = "spring";
  $default_ewidth = "";
  # PIPE VALUES
  $pipe = $_REQUEST['pipe'];
  $graph_file = $_REQUEST['graph_file'];
  $graph_format = $_REQUEST['graph_format'];
  $scol = $_REQUEST['scol'];
  $tcol = $_REQUEST['tcol'];
  $wcol = $_REQUEST['wcol'];  
  $eccol = $_REQUEST['eccol'];
  $sccol = $_REQUEST['sccol'];  
  $tccol = $_REQUEST['tccol'];
  
  # demo graph
  $demo = $_REQUEST['demo'];
  if ($demo == 1) {
    $demo_graph = storeFile("demo_files/string_coexpression_mcl_clusters.tab");
    $demo_wcol = 3;
    $demo_eccol = 4;
    $demo_ewidth = "checked";
  }

  title('display-graph');
  echo ("<center>Generate a figure for a graph\n");
   echo ("<br>This program was developed by 
	  <a target=_blank href=http://www.bigre.ulb.ac.be/Users/sylvain/>Sylvain Broh&eacute;e</a> and
          <a target=_blank href=http://www.bigre.ulb.ac.be/Users/jvanheld/>Jacques van Helden</a>.
          </center>\n");
   echo ("<b>Warning!</b> The layout and display facilities supported on the Web
	  site are restricted.<br>\n

          Several specialized programs exists for visualizing networks
	  in a dynamical way, and with a large diversity of layout
	  options. <br>\n
          
          The GML files exported by the NeAT tools can be imported in
	  <a target = '_blank' href =
	  'http://cytoscape.org/'><b>Cytoscape</b></a> or <a target =
	  '_blank' href =
	  'http://tyna.gersteinlab.org/tyna/'><b>yED</b></a>.  <br>\n

          For a flexible manipulation of the graphs, we recommend to
	  use these program. <p>\n" 
   );

  echo ("<form method='post' action='display_graph.php' enctype='multipart/form-data'>");
  echo ("&nbsp;&nbsp;&nbsp;<a href = 'help.display_graph.html#formats'><B>Input format</B></a>");
  if (!$pipe) {
    echo ("
    &nbsp;<select name='in_format'>
    <option selected value = 'tab'> tab-delimited format
    <option value = 'adj_matrix'> Adjacency matrix
    <option value = 'gml'> GML format
    </select><br><br>");
  } else {
    echo ": $graph_format<br>";
    echo "<input type='hidden' NAME='in_format' VALUE='$graph_format'>";
  }
  echo ("&nbsp;&nbsp;&nbsp;<B><a href = 'help.display_graph.html#formats'>Output format</a></B>&nbsp;<select name='out_format'>
  <option selected value = 'png'> PNG
  <option value = 'ps'> PS
  <option value = 'jpeg'> JPEG
  </select><br><br>");
  if (!$pipe) {
    if ($demo) {
      demo("This demonstration graph consists in the set of clusters extracted by the <a href = 'http://micans.org/mcl/' target = '_blank'>MCL clustering algorithm</a> from a subset of the <a href = 'http://string.embl.de/' target = '_blank'>STRING database</a> (database evidence channel) with inflation 1.8. It consists in 9991 interactions between 1231 proteins. Each edge color corresponds to a different cluster returned by the MCL algorithm. As the graph is weighted the width of the edges will be proportional to the weight.");
    }

    echo ("<textarea name='graph' rows='6' cols='65'>$demo_graph</textarea>");
    echo("<br>Upload graph from file : <br>");
    echo("<input type='file' name='graph_file' size='45' /><br>");
    echo("<br><a href = 'help.display_graph.html#columns'>Column specification (only relevant for tab-delimited input)</a><br>");
    echo("<table>");
    echo("<tr><td><B><a href = 'help.display_graph.html#scol'>Source node</a></B></td><td><input type = 'text' name='s_col' value = '$default_scol' size = 1></input></td></tr>");
    echo("<tr><td><B><a href = 'help.display_graph.html#tcol'>Target node</a></B></td><td><input type = 'text' name='t_col' value = '$default_tcol' size = 1></input></td></tr>");
    echo("<tr><td><B><a href = 'help.display_graph.html#wcol'>Weight or label </a></B></td><td><input type = 'text' name='w_col' value = '$demo_wcol' size = 1></input></td></tr>");
    echo("<tr><td><B><a href = 'help.display_graph.html#eccol'>Edge color</a></B></td><td><input type = 'text' name='ec_col' value = '$demo_eccol' size = 1></input></td></tr>");
    echo("<tr><td><B><a href = 'help.display_graph.html#sccol'>Source node color</a></B></td><td><input type = 'text' name='sc_col' size = 1></input></td></tr>");
    echo("<tr><td><B><a href = 'help.display_graph.html#tccol'>Target node color</a></B></td><td><input type = 'text' name='tc_col' size = 1></input></td></tr>");
    echo("</table>");

  } else {
    info_link("Graph uploaded from the previous treatment", rsat_path_to_url($graph_file));
    echo "<input type='hidden' NAME='pipe_graph_file' VALUE='$graph_file'>";
  }

  if ($graph_format == 'tab') {
    echo "<input type='hidden' NAME='s_col' VALUE='$scol'/>\n";
    echo "<input type='hidden' NAME='t_col' VALUE='$tcol'/>\n";
    echo "<input type='hidden' NAME='w_col' VALUE='$wcol'/>\n";
    echo "<input type='hidden' NAME='ec_col' VALUE='$eccol'/>\n";
    echo "<input type='hidden' NAME='sc_col' VALUE='$sccol'/>\n";
    echo "<input type='hidden' NAME='tc_col' VALUE='$tccol'/>\n";
  }


  echo ("<B><a href = 'help.convert_graph.html#layout'>Layout</a></B>
    <select name='layout'>
      <option value='none'> No layout
      <option selected value='spring'>Spring embedding
      <option value='random'> Random
    </select><br> "
      );

/*  echo ("
  <input type='checkbox' name='layout' value='on' $default_layout/>&nbsp;<B><a href = 'help.display_graph.html#layout'>Calculate the layout of the nodes (mandatory for all input formats except GML)</a></B><br>");*/

  echo ("<input type='checkbox' $demo_ewidth name='ewidth' value='on' />&nbsp;<B><a href = 'help.display_graph.html#ewidth'>Edge width proportional to weight</a></B><br>

  <ul><ul><table class='formbutton'>
  <TD><input type='submit' name='.submit' value='GO' /></TD>
  <TD><B><A HREF='display_graph_form.php?demo=0'>RESET</A></B></TD>
  <TD><B><A HREF='display_graph_form.php?demo=1'>DEMO</A></B></TD>
  </form>
  <TD><B><A HREF='help.display_graph.html'>MANUAL</A></B></TD>
  <TD><B><A target = '_blank' HREF='".checkNeatTutorial("tutorials/neat_tutorial/Network_visualization_forma.html")."'>TUTORIAL</A></B></TD>
  <TD><B><A HREF='mailto:sylvain@bigre.ulb.ac.be'>MAIL</A></B></TD>
  </TR></TABLE></ul></ul>
 ");

?>
</body>
</html>
