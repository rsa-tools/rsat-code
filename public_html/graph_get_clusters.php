<html>
<head>
   <title>Network Analysis Tools - graph-get-clusters</title>
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
  $cmd_file = getTempFileName('commands_graph_getclusters');
  $cmd_handle = fopen($cmd_file, 'a');

  title('graph-get-clusters - results');
  # Error status
  $error = 0;
  # Get parameters
  $in_format = $_REQUEST['in_format'];
  $out_format = $_REQUEST['out_format'];
  if ($_FILES['graph_file']['name'] != "") {
    $graph_file = uploadFile('graph_file');
  } else if ($_REQUEST['pipe_graph_file'] != "")  {
    $graph_file = $_REQUEST['pipe_graph_file'];
  }
  if ($_FILES['clusters_file']['name'] != "") {
    $clusters_file = uploadFile('clusters_file');
  } 
  $now = date("Ymd_His");
  $graph = $_REQUEST['graph'];
  $return =  $_REQUEST['return'];
  $clusters = $_REQUEST['clusters'];
  if ($clusters == "" && $_REQUEST['pipe_clusters_file']) {
    $pipe_clusters_file = $_REQUEST['pipe_clusters_file'];
    $clusters = storeFile($pipe_clusters_file);
  }
  $distinct = $_REQUEST['distinct'];
  $induced = $_REQUEST['induced'];
  if ($induced == 'on') {
    $induced = 1;
  } else {
    $induced = 0;
  }
  $s_col = $_REQUEST['s_col'];
  $t_col = $_REQUEST['t_col'];
  $w_col = $_REQUEST['w_col'];
  
  
  
  ## If a file and a graph are submitted -> error
  if ($graph != "" && $graph_file != "") {
    $error = 1;
    error("You must not submit both a graph and a graph file");
  }
  ## If clusters and a cluster file are submitted -> error
  if ($clusters != "" && $cluster_file != "") {
    $error = 1;
    error("You must not submit both clusters and a cluster file");
  }
  ## No specification of the source and target columns
  if ($in_format == "tab" && $s_col == "" && $t_col == "") {
    warning("Default value for source and target columns for tab-delimited input format are 1 and 2 respectively");
  }
  ## put the content of the file $graph_file in $graph
  if ($graph_file != "" && $graph == "") {
    $graph = storeFile($graph_file);
  }
  ## put the content of the file $clusters_file in $clusters
  if ($clusters_file != "" && $clusters == "") {
    $clusters = storeFile($clusters_file);
  }
  ## If no graph are submitted -> error
  if ($graph == "" && $graph_file == "") {
    $error = 1;
    error("You must submit an input graph");
  }
  ## If no clusters are submitted -> error
  if ($clusters == "") {
    $error = 1;
    error("You must submit input clusters");
  }  
  ## If no clusters are submitted -> error
  if ($induced && $return == 'graph') {
    $error = 1;
    error("You must not ask for an annotated graph when inducing the graph");
  }    
  if (!$error) { 
    $graph = trim_text($graph);
    $clusters = trim_text($clusters);
    ## Load the parameters of the program in to an array
    $parameters = array( 
      "request" => array(
        "informat"=>$in_format,
        "outformat"=>$out_format,
        "return"=>$return,
        "clusters"=>$clusters,
        "inputgraph"=>$graph,
        "scol"=>$s_col,
        "tcol"=>$t_col,
        "wcol"=>$w_col,
        "distinct"=>$distinct,
        "induced"=>$induced
      )
    );
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
    # Execute the command and catch the errors
    try {
      $echoed = $client->graph_get_clusters($parameters);
      $soap_error = 0;
    } catch (Exception $soap_exception) {
      echo ("<pre>");
      echo "Error : \n\t",  $soap_exception->getMessage(), "\n";
      echo ("</pre>");
      $soap_error = 1;
      exit(1);
    }  


    $response =  $echoed->response;
    $command = $response->command;
    $server = $response->server;
    echo ("<p><b>Command :</b> $command</p>");
    $client = $response->client;
#    $server = rtrim ($server);
#    $temp_file = explode('/',$server);
#    $temp_file = end($temp_file);
#    $resultURL = $WWW_RSA."/tmp/".$temp_file;
    store_command($command, "graph comparison", $cmd_handle);
    $URL['Result'] = rsat_path_to_url($server);

    hourglass("off");


    ## Close command handle
    fclose($cmd_handle);
    $URL['Server commands'] = rsat_path_to_url($cmd_file);

    ## DISPLAY THE RESULT
    print_url_table($URL);

    # Display the results
    echo "The results is available at the following URL ";
    echo "<a href = '$resultURL'>$resultURL</a>"; 
    echo "<hr>\n";
    if ($return == 'clusters' || $return == 'graph') {
      echo "
        <TABLE CLASS = 'nextstep'>
      <TR>
        <Th colspan = 3>
          Next step
        </Th>
      </TR>
      <TR>
        <TD>
          <FORM METHOD='POST' ACTION='display_graph_form.php'>
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
          <FORM METHOD='POST' ACTION='compare_graphs_form.php'>
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
            <INPUT type='submit' value='Compare this graph to another one'>
          </form>
        </td>
        <TD>
          <FORM METHOD='POST' ACTION='random_graph_form.php'>
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
        <FORM METHOD='POST' ACTION='convert_graph_form.php'>
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
        <FORM METHOD='POST' ACTION='graph_topology_form.php'>
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
        <FORM METHOD='POST' ACTION='graph_neighbours_form.php'>
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
        <FORM METHOD='POST' ACTION='mcl_form.php'>
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
        <FORM METHOD='POST' ACTION='alter_graph_form.php'>
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
        <FORM METHOD='POST' ACTION='pathfinder_form.php'>
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
          <FORM METHOD='POST' ACTION='visant.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='visant_graph_file' VALUE='$server'>
          <input type='hidden' NAME='visant_graph_format' VALUE='$out_format'>
          <input type='hidden' NAME='visant_directed' VALUE='$directed'>
          <input type='hidden' NAME='tab_java' VALUE='0'>";
          echo "
          <INPUT type='submit' value='Load in VisANT'>
          </form>
        </td>
      </tr>
    </table>";
    }

  }
?>
