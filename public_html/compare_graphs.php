<html>
<head>
   <title>Network Analysis Tools - compare-graphs</title>
   <link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
      <style type="text/css">
    <!--
    div.hourglass{position: absolute; top: 80px; left: 400px }
    div.hide{position: absolute; top: 80px; left: 400px }
   -->
</style>
</head>
<body class="results">
<?php 
  require ('functions.php');
  # log file update
  UpdateLogFile("neat","","");

  # File to store the commands
  $cmd_file = getTempFileName('compgraph_cmd', '.txt');
  $cmd_handle = fopen($cmd_file, 'a');

  title('compare-graphs - results');
  # Error status
  $error = 0;
  # Get parameters
  $in_formatQ = $_REQUEST['in_formatQ'];
  $in_formatR = $_REQUEST['in_formatR'];
  
  
  $out_format = $_REQUEST['out_format'];
  if ($_FILES['graph_fileQ']['name'] != "") {
    $graph_fileQ = uploadFile('graph_fileQ');
  } else if ($_REQUEST['pipe_query_graph_file'] != "") {
    $graph_fileQ = $_REQUEST['pipe_query_graph_file'];
  }
  if ($_FILES['graph_fileR']['name'] != "") {
    $graph_fileR = uploadFile('graph_fileR');
  }
  $now = date("Ymd_His");
  $graphQ = $_REQUEST['graphQ'];
  $graphR = $_REQUEST['graphR'];
  $directed = $_REQUEST['directed'];
  $s_colQ = $_REQUEST['s_colQ'];
  $t_colQ = $_REQUEST['t_colQ'];
  $w_colQ = $_REQUEST['w_colQ'];
  $s_colR = $_REQUEST['s_colR'];
  $t_colR = $_REQUEST['t_colR'];
  $w_colR = $_REQUEST['w_colR'];
  if ($directed == "on") {
    $directed = 1;
  }
  $directed = 0;
  $self = $_REQUEST['self'];  
  if ($self == "on") {
    $self = 1;
  }

  $return =  $_REQUEST['return'];
  $outweight =  $_REQUEST['outweight']; 
  ## If a query graph file and a query graph are submitted -> error
  if ($graphQ != "" && $graph_fileQ!= "") {
    $error = 1;
    error("You cannot submit both a query graph and a query graph file");
  }
  ## If a query graph file and a query graph are submitted -> error
  if ($graphR != "" && $graph_fileR!= "") {
    $error = 1;
    error("You cannot submit both a reference graph and a reference graph file");
  }
  ## No specification of the source and target columns
  if ($in_format_Q == "tab" && $s_colQ == "" && $t_colQ == "") {
    warning("Default value for source and target columns for query graph in tab-delimited input format are 1 and 2");
  }
  ## No specification of the source and target columns
  if ($in_format_R == "tab" && $s_colR == "" && $t_colR == "") {
    warning("Default value for source and target columns for reference graph in tab-delimited input format are 1 and 2");
  }
  ## put the content of the file $graph_file in $graph
  if ($graph_fileQ != "" && $graphQ == "") {
    $graphQ = storeFile($graph_fileQ);
  } 
  ## put the content of the file $graph_file in $graph
  if ($graph_fileR != "" && $graphR == "") {
    $graphR = storeFile($graph_fileR);
  }
  ## If no graph are submitted -> error
  if ($graphQ == "" && $graph_fileQ == "") {
    $error = 1;
    error("You must submit a query input graph");
  }
  ## If no graph are submitted -> error
  if ($graphR == "" && $graph_fileR == "") {
    $error = 1;
    error("You must submit a reference input graph");
  }
  
   if (!$error) { 
     $graphQ = trim_text($graphQ);
     $graphR = trim_text($graphR);
     $parameters = array( 
       "request" => array (
         "Qinformat"=>$in_formatQ,
         "Rinformat"=>$in_formatR,
         "outformat"=>$out_format,
         "outweight"=>$outweight,
         "Rinputgraph"=>$graphR,
         "Qinputgraph"=>$graphQ,
         "Qwcol"=>$w_colQ,
         "Qscol"=>$s_colQ,
         "Qtcol"=>$t_colQ,
         "Rwcol"=>$w_colR,
         "Rscol"=>$s_colR,
         "Rtcol"=>$t_colR,
         "directed"=>$directed,
         "return"=>$return,
         "self"=>$self
       )
     );

     // report parameters if echo requested
     if ($properties['rsat_echo'] >= 2) {
       info("rsat_main\t".$rsat_main);
       info("WWW_RSA\t".$WWW_RSA);
       info("neat_wsdl\t".$neat_wsdl);
       print_r($parameters);
       if ($properties['rsat_echo'] >= 3) {
	 phpinfo();
       }
     }

    # Info message
    info("Results will appear below");
    echo"<hr>\n";
    hourglass("on");

    # Open the SOAP client
    $client = new SoapClient(
			     $neat_wsdl,
			     array(
				   'trace' => 1,
				   'soap_version' => SOAP_1_1,
				   'style' => SOAP_DOCUMENT,
				   'encoding' => SOAP_LITERAL
				   )
			     );
    flush();
    # Execute the command
#    echo "<pre>";
    $echoed = $client->compare_graphs($parameters);
#    echo "</pre>"; 
    $response =  $echoed->response;
    $command = $response->command;
    $server = $response->server;
    $client = $response->client;
    store_command($command, "graph comparison", $cmd_handle);
    $URL['Result'] = rsat_path_to_url($server);

    # The comment file has the same name as the
    # result file with ".comments" at the end of the string.
    $comments_temp_file = $server.".comments";
    $comments = storeFile($comments_temp_file);

    hourglass("off");

    ## Close command handle
    fclose($cmd_handle);
    if ($properties['rsat_echo'] >= 1) {
      $URL['Server commands'] = rsat_path_to_url($cmd_file);
    }

    info("server commands\t".$URL['Server commands']);

    ## Print the comparison statistics
    echo "<pre>";
    echo "$comments";
    echo "</pre><hr>";

    ## DISPLAY THE RESULT
    print_url_table($URL);
    
    echo "<hr>\n";
  echo "
  <TABLE CLASS = 'nextstep'>
    <TR>
      <Th colspan = 3>
        Next step
      </Th>
    </TR>
    <TR>
      <TD>
        <FORM METHOD='POST' ACTION='${neat_www_root}/display_graph_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$server'>
          <input type='hidden' NAME='graph_format' VALUE='$out_format'>";
          if ($out_format == 'tab') {
            echo "
            <input type='hidden' NAME='scol' VALUE='1'>
            <input type='hidden' NAME='tcol' VALUE='2'>
            <input type='hidden' NAME='wcol' VALUE='3'>
            <input type='hidden' NAME='eccol' VALUE='4'>";
          }
          echo "
          <INPUT type='submit' value='Display the graph'>
        </form>
      </td>
      <TD>
        <FORM METHOD='POST' ACTION='${neat_www_root}/convert_graph_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$server'>
          <input type='hidden' NAME='graph_format' VALUE='$out_format'>";
          if ($out_format == 'tab') {
            echo "
            <input type='hidden' NAME='scol' VALUE='1'>
            <input type='hidden' NAME='tcol' VALUE='2'>
            <input type='hidden' NAME='wcol' VALUE='3'>
            <input type='hidden' NAME='eccol' VALUE='4'>";
          }
          echo "
          <INPUT type='submit' value='Convert $out_format to another format'>
        </form>
      </td>
      <TD>
        <FORM METHOD='POST' ACTION='${neat_www_root}/random_graph_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$server'>
          <input type='hidden' NAME='graph_format' VALUE='$out_format'>";
          if ($out_format == 'tab') {
            echo "
            <input type='hidden' NAME='scol' VALUE='1'>
            <input type='hidden' NAME='tcol' VALUE='2'>
            <input type='hidden' NAME='wcol' VALUE='3'>";
          }
          echo "
          <INPUT type='submit' value='Randomize this graph'>
        </form>
      </td>
    </tr>
    <tr>
      <TD>
        <FORM METHOD='POST' ACTION='${neat_www_root}/graph_get_clusters_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$server'>
          <input type='hidden' NAME='graph_format' VALUE='$out_format'>";
          if ($out_format == 'tab') {
            echo "
            <input type='hidden' NAME='scol' VALUE='1'>
            <input type='hidden' NAME='tcol' VALUE='2'>
            <input type='hidden' NAME='wcol' VALUE='3'>";
          }
          echo "
          <INPUT type='submit' value='Map clusters or extract a subnetwork'>
        </form>
      </td>
      <TD>
        <FORM METHOD='POST' ACTION='${neat_www_root}/graph_topology_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$server'>
          <input type='hidden' NAME='graph_format' VALUE='$out_format'>";
          if ($out_format == 'tab') {
            echo "
            <input type='hidden' NAME='scol' VALUE='1'>
            <input type='hidden' NAME='tcol' VALUE='2'>
            <input type='hidden' NAME='wcol' VALUE='3'>";
          }
          echo "
          <INPUT type='submit' value='Nodes topology statistics'>
        </form>
      </td>
      <TD>
        <FORM METHOD='POST' ACTION='${neat_www_root}/graph_neighbours_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$server'>
          <input type='hidden' NAME='graph_format' VALUE='$out_format'>";
          if ($out_format == 'tab') {
            echo "
            <input type='hidden' NAME='scol' VALUE='1'>
            <input type='hidden' NAME='tcol' VALUE='2'>
            <input type='hidden' NAME='wcol' VALUE='3'>";
          }
          echo "
          <INPUT type='submit' value='Neighbourhood analysis'>
       </form>
      </td>    
    </tr>
    <TR>
      <TD>
        <FORM METHOD='POST' ACTION='${neat_www_root}/mcl_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$server'>
          <input type='hidden' NAME='graph_format' VALUE='$out_format'>";
          if ($out_format == 'tab') {
            echo "
            <input type='hidden' NAME='scol' VALUE='1'>
            <input type='hidden' NAME='tcol' VALUE='2'>
            <input type='hidden' NAME='wcol' VALUE='3'>";
          }
          echo "
          <INPUT type='submit' value='MCL Graph clustering'>
        </form>
      </td>
      <TD>
        <FORM METHOD='POST' ACTION='${neat_www_root}/alter_graph_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$server'>
          <input type='hidden' NAME='graph_format' VALUE='$out_format'>";
          if ($out_format == 'tab') {
            echo "
            <input type='hidden' NAME='scol' VALUE='1'>
            <input type='hidden' NAME='tcol' VALUE='2'>
            <input type='hidden' NAME='wcol' VALUE='3'>";
          }
          echo "
          <INPUT type='submit' value='Graph alteration'>
        </form>
      </td>
      <TD>
        <FORM METHOD='POST' ACTION='${neat_www_root}/pathfinder_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$server'>
          <input type='hidden' NAME='in_format' VALUE='$out_format'>";
          echo "
          <INPUT type='submit' value='Path Finding'>
        </form>
      </td>
      </tr>
      <tr>
        <TD>
          <FORM METHOD='POST' ACTION='${neat_www_root}/visant.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='visant_graph_file' VALUE='$server'>
          <input type='hidden' NAME='visant_graph_format' VALUE='$out_format'>
          <input type='hidden' NAME='visant_directed' VALUE='$directed'>
          <input type='hidden' NAME='tab_java' VALUE='0'>";
          echo "
          <INPUT type='submit' value='Load in VisANT'>
          </form>
        </td>";
	  if (($out_format == 'tab')&&($outweight != "Q::R")){
	    echo "
      <TD>
        <FORM METHOD='POST' ACTION='${neat_www_root}/roc-stats_form.cgi'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='roc-stats_graph_file' VALUE='$resultURL'>
          <input type='hidden' NAME='sc_col' VALUE='3'>
          <input type='hidden' NAME='status_col' VALUE='5'>
          <input type='hidden' NAME='pos' VALUE=\"Q.and.R\nR.not.Q\">
	  <input type='hidden' NAME='neg' VALUE=\"Q.not.R\">
	  <input type='hidden' NAME='null' VALUE=\"\<NULL\>\">";
	    echo "
          <INPUT type='submit' value='ROC stats'>
        </form>
      </td>";
	  }
	  echo "
      </tr>
  </table>";
   }
?>
</HTML>
