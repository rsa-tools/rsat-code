<html>
<head>
   <title>GrA-tools - graph-node-degree</title>
   <link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
</head>
<body class="results">
<?php 
  require ('functions.php');
  title('graph-node-degree - results');
  # Error status
  $error = 0;
  # Get parameters
  $in_format = $_REQUEST['in_format'];
  if ($_FILES['graph_file']['name'] != "") {
    $graph_file = uploadFile('graph_file');
  }
  if ($_FILES['nodes_file']['name'] != "") {
    $nodes_file = uploadFile('nodes_file');
  }
  $now = date("Ymd_His");
  $graph = $_REQUEST['graph'];
  $nodes = $_REQUEST['nodes'];
  $s_col = $_REQUEST['s_col'];
  $t_col = $_REQUEST['t_col'];
  $all = 0;
  ## If a file and a graph are submitted -> error
  if ($graph != "" && $graph_file != "") {
    $error = 1;
    error("You must not submit both a graph and a graph file");
  }
  ## If nodes and a node file are submitted -> error
  if ($nodes != "" && $nodes_file != "") {
    $error = 1;
    error("You must not submit both nodes and a nodes file");
  }
  ## No specification of the source and target columns
  if ($in_format == "tab" && $s_col == "" && $t_col == "") {
    warning("Default value for source and target columns for tab-delimited input format are 1 and 2 respectively");
  }
  ## put the content of the file $graph_file in $graph
  if ($graph_file != "" && $graph == "") {
    $graph = storeFile($graph_file);
  }
  ## put the content of the file $nodes_file in $nodes
  if ($nodes_file != "" && $nodes == "") {
    $nodes = storeFile($nodes_file);
  }
  ## If no graph are submitted -> error
  if ($graph == "" && $graph_file == "") {
    $error = 1;
    error("You must submit an input graph");
  }
  ## If no nodes are submitted -> info : all nodes will be computed
  if ($nodes == "") {
   info("You did not submit any nodes. The degree of all nodes will be computed");
   $all = 1;
  }  
  
  if (!$error) { 
    $graph = trim_text($graph);
    $nodes = trim_text($nodes);
    ## Load the parameters of the program in to an array
    $parameters = array( 
      "request" => array(
        "informat"=>$in_format,
        "nodefile"=>$nodes,
        "inputgraph"=>$graph,
        "scol"=>$s_col,
        "tcol"=>$t_col,
        "all"=>$all
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
    echo "<pre>";
    $echoed = $client->graph_node_degree($parameters);

    $response =  $echoed->response;
    $command = $response->command;
    $server = $response->server;
    $client = $response->client;
    echo "</pre>";
    $server = rtrim ($server);
    $temp_file = explode('/',$server);
    $temp_file = end($temp_file);
    $resultURL = $WWW_RSA."/tmp/".$temp_file;
    # Display the results
    echo "The results is available at the following URL ";
    echo "<a href = '$resultURL'>$resultURL</a>"; 
    echo "<hr>\n";
   
  }
?>
