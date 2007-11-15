<html>
<head>
   <title>GrA-tools - graph-get-clusters</title>
   <link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
</head>
<body class="results">
<?php 
  require ('functions.php');
  title('graph-get-clusters - results');
  # Error status
  $error = 0;
  # Get parameters
  $in_format = $_REQUEST['in_format'];
  $out_format = $_REQUEST['out_format'];
  if ($_FILES['graph_file']['name'] != "") {
    $graph_file = uploadFile('graph_file');
  }
  if ($_FILES['clusters_file']['name'] != "") {
    $clusters_file = uploadFile('clusters_file');
  }
  $now = date("Ymd_His");
  $graph = $_REQUEST['graph'];
  $return =  $_REQUEST['return'];
  $clusters = $_REQUEST['clusters'];
  $distinct = $_REQUEST['distinct'];
  if ($layout == 'on') {
    $layout = 1;
  }
  $induced = $_REQUEST['induced'];
  if ($layout == 'on') {
    $induced = 1;
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
  if ($clusters == "" && $clusters == "") {
    $error = 1;
    error("You must submit an input graph");
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
        "distinct"=>$layout,
        "induced"=>$undirected
      )
    );
    # Info message
    info("Results will appear below");
    echo"<hr>\n";
  
    # Open the SOAP client
    $client = new SoapClient(
                       'http://rsat.scmbb.ulb.ac.be/rsat/web_services/RSATWS.wsdl',
                           array(
                                 'trace' => 1,
                                 'soap_version' => SOAP_1_1,
                                 'style' => SOAP_DOCUMENT,
                                 'encoding' => SOAP_LITERAL
                                 )
                           );
    # Execute the command
    $echoed = $client->graph_get_clusters($parameters);

    $response =  $echoed->response;
    $command = $response->command;
    $server = $response->server;
    $client = $response->client;
    $server = rtrim ($server);
    $temp_file = explode('/',$server);
    $temp_file = end($temp_file);
    $resultURL = $WWW_RSA."/tmp/".$temp_file;
    # Display the results
    echo "The results is available at the following URL ";
    echo "<a href = '$resultURL'>$resultURL</a>"; 
    echo "<hr>\n";
    echo("
<TABLE CLASS = 'nextstep'>
    <TR>
      <Th colspan = 3>
        Next step
      </Th>
    </TR>
    <TR>
      <td>
      </td>
      <TD>
        <FORM METHOD='POST' ACTION='display_graph_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$server'>
          <input type='hidden' NAME='graph_format' VALUE='$out_format'>");
          if ($out_format == 'tab') {
            echo "
            <input type='hidden' NAME='scol' VALUE='1'>
            <input type='hidden' NAME='tcol' VALUE='2'>
            <input type='hidden' NAME='wcol' VALUE='3'>";
          }
          echo "
          <INPUT type='submit' value='Display the graph'>
        </form>
      </td>
    </tr>
  </table>";
  }
?>
