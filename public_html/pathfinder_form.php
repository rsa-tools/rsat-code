<html>
<head>
   <title>NAT - Pathfinder</title>
   <link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
</head>
<body class="form">
<?
  require ('functions.php');
  require ('Scer_biocyc.php');

  // default variables
  $default_in_format = "";
  $default_out_format = "";
  $default_sources = "";
  $default_targets = "";
  $default_graph = "";
  $default_graph_id = "";
  $default_directed = 1;
  $default_rank = 5;
  $default_weight = "unit";
  $default_store_graph = 0;

  // variables given within workflow
  $pipe = $_REQUEST['pipe'];
  $requested_graph_file = $_REQUEST['graph_file'];
  $requested_in_format = $_REQUEST['in_format'];

  # demo graph
  $demo = $_REQUEST['demo'];
  if ($demo == 1) {
     $demo_graph = $scer_biocyc;
     $demo_sources = "PAPSSULFOTRANS-RXN/SULFITE-REDUCTASE-(FERREDOXIN)-RXN";
     $demo_targets = "HOMO-SER/SER";
     $default_sources = $demo_sources;
     $default_targets = $demo_targets;
  }

  title('Pathfinder');
  echo ("<center>Do multiple-to-multiple end path finding.</center>\n");
  echo ("<form method='post' action='pathfinder.php' enctype='multipart/form-data'>
  <h2>Input/Output</h2>
  <HR SIZE='2' NOSHADE>
  <br>
  <B><a href='help.pathfinder.html#graph'>1. Input graph</a></B>
  <br>
  <br>");
  if(!$pipe){
  	echo("
  	Enter the graph here (tab-delimited or gml):
  	<br>
  	<textarea name='graph' rows='6' cols='65'>$demo_graph</textarea>
  	<br>
  	<br>
  	<b>Or</b> upload the graph from a file:
  	<br>
  	<input type='file' name='graph_file' size='45' />
  	<br>
  	<br><b>OR</b> enter a graph id (id of previously submitted graph):
  	<br>
  	<input type='text' name='graph_id' size='45' />
  	<br>
  	<br>");
  } else {
  	info("Graph has been stored from previous treatment.");
  	echo "<input type='hidden' NAME='pipe_graph_file' VALUE='$requested_graph_file'>";
  }
   echo("
   <B><a href='help.pathfinder.html#terminals'>2. Seed nodes</a></B>
   <br>
   <br>
   Enter source and target nodes:
  <br>
  <br>
  <table>
     <tr><td>Source nodes&nbsp;&nbsp;</td>        <td><input type='text' name='sources' value='$default_sources' size=100></input></td></tr>
     <tr><td>Target nodes</td>          <td><input type= 'text' name='targets' value='$default_targets' size=100></input></td></tr>
  </table>
  <br>
  <b>Or</b> upload a batch file:
   <br>
   <br>
  <input type='file' name='batch_file' size='45' />
  <br>
  <br>
   <B><a href='help.pathfinder.html#formats'>3. Input/Output graph formats</a></B>
  <br>");
  if(!$pipe){
    echo("<br>
       Input format&nbsp;
       <select name='in_format'>
       <option selected value = 'flat'> tab-delimited format
       <option value = 'gml'> GML format
       </select>
       <br>");
  }else{
  		info("Input graph format: ".$requested_in_format);
  		 echo "<input type='hidden' NAME='pipe_in_format' VALUE='$requested_in_format'>";
  }
  echo("<br><br>
  Output format&nbsp;
   <select name='out_format'>
  <option selected value = 'flat'> tab-delimited format
  <option value = 'GML'> GML format
  </select>
  <br>
  <br>
  <B><a href = 'help.pathfinder.html#output'>4. Output type</a></B>
  <br>
  <br>
   <select name='outputType'>
  <option selected value = 'pathsTable'>table of paths
  <option value = 'pathsMarked'>input graph with paths highlighted
  <option value = 'pathsGraphs'>paths separated in one graph
  <option value = 'pathsUnion'>paths unified in one graph
  </select>
  <br>
  <br>
  <h2>Options</h2>
  <HR SIZE='2' NOSHADE>
  <br>
  <table>
     <tr><td><B><a href = 'help.pathfinder.html#rank'>Rank</a></B></td>                   <td><input type = 'text' name='rank' value = '$default_rank' size = 10></input></td></tr>
     ");
     if(!$pipe){
     echo("<tr><td><B><a href = 'help.pathfinder.html#weight'>Weighting scheme</a></B></td>     <td><select name='weight'>
                <option value = 'unit'>unit weight
                <option value = 'none'>as given in input graph
                <option selected value = 'con'>degree of nodes as weight
     	</select></td></tr>");
     }else{
     	info("No weight selection available for already submitted graphs");
     }
     echo("<tr><td><B><a href = 'help.pathfinder.html#directed'>Directed</a></B></td> <td><input type='checkbox' name='directed' value='on'></input></td></tr>
     <tr><td><B><a href = 'help.pathfinder.html#server'>Store graph on server</a></B></td> <td><input type='checkbox' name='store_graph' value='on'></input></td></tr>
  </table>
  <h2>Launch</h2>
  <HR SIZE='2' NOSHADE>
  <br>
  <table class='formbutton'>
  <TD><input type='submit' name='.submit' value='GO' /></TD>
  <TD><B><A HREF='pathfinder_form.php?demo=0'>RESET</A></B></TD>
  <TD><B><A HREF='pathfinder_form.php?demo=1'>DEMO</A></B></TD>
  </form>
  <TD><B><A HREF='help.pathfinder.html'>MANUAL</A></B></TD>
  <TD><B><A HREF='mailto:kfaust@ulb.ac.be'>MAIL</A></B></TD>
   <TD><B><A HREF='help.pathfinder.html#credits'>CREDITS</A></B></TD>
  </TR></TABLE></ul></ul>
 ");
?>
</body>
</html>