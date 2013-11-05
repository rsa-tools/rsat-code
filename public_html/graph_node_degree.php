<html>
<head>
   <title>NeAT - graph-node-degree</title>
   <link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
</head>
<body class="results">
<?php 
  require ('functions.php');
  title('graph-node-degree - results');
  # log file update
  UpdateLogFile("neat","","");
  # Error status
  $error = 0;
  # Get parameters
  $in_format = $_REQUEST['in_format'];
  if ($_FILES['graph_file']['name'] != "") {
    $graph_file = uploadFile('graph_file');
  } else if ($_REQUEST['pipe_graph_file'] != "")  {
    $graph_file = $_REQUEST['pipe_graph_file'];
  }
  if ($_FILES['nodes_file']['name'] != "") {
    $nodes_file = uploadFile('nodes_file');
  }
  $now = date("Ymd_His");
  $graph = $_REQUEST['graph'];
  $nodes = $_REQUEST['nodes'];
  $s_col = $_REQUEST['s_col'];
  $t_col = $_REQUEST['t_col'];
  
  $all = $_REQUEST['allnodes'];
  if ($all == 'all') {
    $all = 1;
  } else {
    $all = 0;
  }
  
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
    $soap_client = new SoapClient(
                       $neat_wsdl,
                           array(
                                 'trace' => 1,
                                 'soap_version' => SOAP_1_1,
                                 'style' => SOAP_DOCUMENT,
                                 'encoding' => SOAP_LITERAL
                                 )
                           );
    # Execute the command
    echo "<pre>";
    $echoed = $soap_client->graph_node_degree($parameters);

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
     
   $graph_node_degree_result = storeFile($server);
   ### CLASSFREQ + XY-GRAPH (all nodes)
   # classfreq
   $cf_all_parameters = array( 
      "request" => array(
        "inputFile"=>$graph_node_degree_result,
        "col"=>4,
        "classinterval"=>1
      )
    );
    echo "<pre>";
    $cf_all_echoed = $soap_client->classfreq($cf_all_parameters);

    $cf_all_response =  $cf_all_echoed->response;
    $cf_all_command = $cf_all_response->command;
    $cf_all_server = $cf_all_response->server;
    $cf_all_client = $cf_all_response->client;
    $cf_all_server = rtrim ($cf_all_server);
    $cf_all_temp_file = explode('/',$cf_all_server);
    $cf_all_temp_file = end($cf_all_temp_file);
    $cf_all_resultURL = "tmp/".$cf_all_temp_file;
    echo "</pre>";
    $cf_all_server = rtrim ($cf_all_server);
    
    # xy graph
    $cf_all_results = storeFile($cf_all_server);
    $xy_all_parameters = array( 
       "request" => array(
         "inputFile"=>$cf_all_results,
         "xcol"=>"2",
         "ycol"=>"4,6",
         "xmin"=>0,
         "format"=>"png",
         "lines"=>1,
         "xleg1"=>"Degree",
         "yleg1"=>"Number of nodes",
         "title1"=>"Degree distribution",
         "legend"=>1,
         "header"=>1
       )
     );
     echo "<pre>";
     $xy_all_echoed = $soap_client->xygraph($xy_all_parameters);
     $xy_all_response =  $xy_all_echoed->response;
     $xy_all_command = $xy_all_response->command;
     $xy_all_server = $xy_all_response->server;
     $xy_all_client = $xy_all_response->client;

     echo "</pre>";
     $xy_all_server = rtrim ($xy_all_server);
         $xy_all_temp_file = explode('/',$xy_all_server);
    $xy_all_temp_file = end($xy_all_temp_file);
    $xy_all_resultURL = "tmp/".$xy_all_temp_file;
    
    
    # xy graph (log)
    $xy_all_parameters_log = array( 
       "request" => array(
         "inputFile"=>$cf_all_results,
         "xcol"=>"2",
         "ycol"=>"4,6",
         "xmin"=>0,
         "format"=>"png",
         "lines"=>1,
         "xleg1"=>"Degree",
         "yleg1"=>"Number of nodes",
         "title1"=>"Degree distribution",
         "legend"=>1,
         "header"=>1,
         "xlog"=>10,
         "ylog"=>10
       )
     );
     echo "<pre>";
     $xy_all_echoed_log = $soap_client->xygraph($xy_all_parameters_log);
     $xy_all_response_log =  $xy_all_echoed_log->response;
     $xy_all_command_log = $xy_all_response_log->command;
     $xy_all_server_log = $xy_all_response_log->server;
     $xy_all_client_log = $xy_all_response_log->client;

     echo "</pre>";
     $xy_all_server_log = rtrim ($xy_all_server_log);
     $xy_all_temp_file_log = explode('/',$xy_all_server_log);
     $xy_all_temp_file_log = end($xy_all_temp_file_log);
     $xy_all_resultURL_log = "tmp/".$xy_all_temp_file_log;

     
   ### CLASSFREQ + XY-GRAPH (intra nodes degree)
   # classfreq
   $cf_in_parameters = array( 
      "request" => array(
        "inputFile"=>$graph_node_degree_result,
        "col"=>2,
        "classinterval"=>1
      )
    );
    echo "<pre>";
    $cf_in_echoed = $soap_client->classfreq($cf_in_parameters);

    $cf_in_response =  $cf_in_echoed->response;
    $cf_in_command = $cf_in_response->command;
    $cf_in_server = $cf_in_response->server;
    $cf_in_client = $cf_in_response->client;
    $cf_in_server = rtrim ($cf_in_server);
    $cf_in_temp_file = explode('/',$cf_in_server);
    $cf_in_temp_file = end($cf_in_temp_file);
    $cf_in_resultURL = "tmp/".$cf_in_temp_file;
    echo "</pre>";
    $cf_in_server = rtrim ($cf_in_server);
    
    # xy graph 
    $cf_in_results = storeFile($cf_in_server);
    $xy_in_parameters = array( 
       "request" => array(
         "inputFile"=>$cf_in_results,
         "xcol"=>"2",
         "ycol"=>"4,6",
         "format"=>"png",
         "lines"=>1,
         "xmin"=>0,
         "title1"=>"In-Degree distribution",
         "xleg1"=>"Degree",
         "yleg1"=>"Number of nodes",
         "legend"=>1,
         "header"=>1
       )
     );
     echo "<pre>";
     $xy_in_echoed = $soap_client->xygraph($xy_in_parameters);
     $xy_in_response =  $xy_in_echoed->response;
     $xy_in_command = $xy_in_response->command;
     $xy_in_server = $xy_in_response->server;
     $xy_in_client = $xy_in_response->client;

     echo "</pre>";
     $xy_in_server = rtrim ($xy_in_server);
         $xy_in_temp_file = explode('/',$xy_in_server);
    $xy_in_temp_file = end($xy_in_temp_file);
    $xy_in_resultURL = "tmp/".$xy_in_temp_file;
  
    # xy graph (log)
    $xy_in_parameters_log = array( 
       "request" => array(
         "inputFile"=>$cf_in_results,
         "xcol"=>"2",
         "ycol"=>"4,6",
         "format"=>"png",
         "lines"=>1,
         "xmin"=>0,
         "title1"=>"In-Degree distribution",
         "xleg1"=>"Degree",
         "yleg1"=>"Number of nodes",
         "legend"=>1,
         "header"=>1,
         "xlog"=>10,
         "ylog"=>10
       )
     );
     echo "<pre>";
     $xy_in_echoed_log = $soap_client->xygraph($xy_in_parameters_log);
     $xy_in_response_log =  $xy_in_echoed_log->response;
     $xy_in_command_log = $xy_in_response_log->command;
     $xy_in_server_log = $xy_in_response_log->server;
     $xy_in_client_log = $xy_in_response_log->client;

     echo "</pre>";
     $xy_in_server_log = rtrim ($xy_in_server_log);
     $xy_in_temp_file_log = explode('/',$xy_in_server_log);
     $xy_in_temp_file_log = end($xy_in_temp_file_log);
     $xy_in_resultURL_log = "tmp/".$xy_in_temp_file_log;  
  
  
  
   ### CLASSFREQ + XY-GRAPH (extra nodes degree)
   # classfreq
   $cf_out_parameters = array( 
      "request" => array(
        "inputFile"=>$graph_node_degree_result,
        "col"=>3,
        "classinterval"=>1
      )
    );
    echo "<pre>";
    $cf_out_echoed = $soap_client->classfreq($cf_out_parameters);

    $cf_out_response =  $cf_out_echoed->response;
    $cf_out_command = $cf_out_response->command;
    $cf_out_server = $cf_out_response->server;
    $cf_out_client = $cf_out_response->client;
    $cf_out_server = rtrim ($cf_out_server);
    $cf_out_temp_file = explode('/',$cf_out_server);
    $cf_out_temp_file = end($cf_out_temp_file);
    $cf_out_resultURL = "tmp/".$cf_out_temp_file;
    echo "</pre>";
    $cf_out_server = rtrim ($cf_out_server);
    
    # xy graph
    $cf_out_results = storeFile($cf_out_server);
    $xy_out_parameters = array( 
       "request" => array(
         "inputFile"=>$cf_out_results,
         "xcol"=>"2",
         "ycol"=>"4,6",
         "format"=>"png",
         "lines"=>1,
         "xmin"=>0,
         "title1"=>"Out-Degree distribution",
         "xleg1"=>"Degree",
         "yleg1"=>"Number of nodes",
         "legend"=>1,
         "header"=>1
       )
     );
     echo "<pre>";
     $xy_out_echoed = $soap_client->xygraph($xy_out_parameters);
     $xy_out_response =  $xy_out_echoed->response;
     $xy_out_command = $xy_out_response->command;
     $xy_out_server = $xy_out_response->server;
     $xy_out_client = $xy_out_response->client;

     echo "</pre>";
     $xy_out_server = rtrim ($xy_out_server);
     $xy_out_temp_file = explode('/',$xy_out_server);
     $xy_out_temp_file = end($xy_out_temp_file);
     $xy_out_resultURL = "tmp/".$xy_out_temp_file;
   
    # xy graph (log)
    $xy_out_parameters_log = array( 
       "request" => array(
         "inputFile"=>$cf_out_results,
         "xcol"=>"2",
         "ycol"=>"4,6",
         "format"=>"png",
         "lines"=>1,
         "xmin"=>0,
         "title1"=>"Out-Degree distribution",
         "xleg1"=>"Degree",
         "yleg1"=>"Number of nodes",
         "legend"=>1,
         "header"=>1,
         "xlog"=>10,
         "ylog"=>10
       )
     );
     echo "<pre>";
     $xy_out_echoed_log = $soap_client->xygraph($xy_out_parameters_log);
     $xy_out_response_log =  $xy_out_echoed_log->response;
     $xy_out_command_log = $xy_out_response_log->command;
     $xy_out_server_log = $xy_out_response_log->server;
     $xy_out_client_log = $xy_out_response_log->client;

     echo "</pre>";
     $xy_out_server_log = rtrim ($xy_out_server_log);
     $xy_out_temp_file_log = explode('/',$xy_out_server_log);
     $xy_out_temp_file_log = end($xy_out_temp_file_log);
     $xy_out_resultURL_log = "tmp/".$xy_out_temp_file_log;   
   
   
     echo "<table>
       <th align = 'center' colspan = 4><b>Global, in- and out- degree distributions</b></th>
       <tr>
         
         <td><a href = '$cf_all_resultURL'>Global degree distribution</a></td>
         <td><a href = '$cf_in_resultURL'>In degree distribution</a></td>
         <td><a href = '$cf_out_resultURL'>Out degree distribution</a></td>
       </tr>
       <tr>
         
         <td><a href = '$xy_all_resultURL'><img src='$xy_all_resultURL' width = '100%'></a></td>
         <td><a href = '$xy_in_resultURL'><img src='$xy_in_resultURL' width = '100%'></a></td>
         <td><a href = '$xy_out_resultURL'><img src='$xy_out_resultURL' width = '100%'></a></td>
       </tr>
       <tr>
         <td><a href = '$xy_all_resultURL_log'><img src='$xy_all_resultURL_log' width = '100%'></a></td>
         <td><a href = '$xy_in_resultURL_log'><img src='$xy_in_resultURL_log' width = '100%'></a></td>
         <td><a href = '$xy_out_resultURL_log'><img src='$xy_out_resultURL_log' width = '100%'></a></td>
       </tr>       
       
     </table>";
  }
?>
