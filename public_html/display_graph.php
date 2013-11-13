<html>
<head>
   <title>GrA-tools - display-graph</title>
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
  
  
  title('display-graph - results');

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
  $now = date("Ymd_His");  ini_set('default_socket_timeout', 3000);
  

  
  
  $graph = $_REQUEST['graph'];
  $layout = $_REQUEST['layout'];
  if ($layout == 'spring') {
    $layout = 1;
  }
  $s_col = $_REQUEST['s_col'];
  $t_col = $_REQUEST['t_col'];
  $w_col = $_REQUEST['w_col'];
  $ec_col = $_REQUEST['ec_col'];
  $sc_col = $_REQUEST['sc_col'];
  $tc_col = $_REQUEST['tc_col'];
  $ewidth = $_REQUEST['ewidth'];
  if ($ewidth == 'on') {
    $ewidth = 1;
  } else {
    $ewidth = 0;
  }  
  ## If a file and a graph are submitted -> error
  if ($graph != "" && $graph_file != "") {
    $error = 1;
    error("You must not submit both a graph and a graph file");
  }

  ## No specification of the source and target columns
  if ($in_format == "tab" && $s_col == "" && $t_col == "") {
    warning("Default value for source and target columns for tab-delimited input format are 1 and 2");
  }

  ## put the content of the file $graph_file in $graph
  if ($graph_file != "" && $graph == "") {
    $graph = storeFile($graph_file);
  }
  ## If no graph are submitted -> error
  if ($graph == "" && $graph_file == "") {
    $error = 1;
    error("You must submit an input graph");
  }
  ## If display is not selected and the format is not GML -> error
  if ($in_format != "gml" && $layout != "1") {
    $error = 1;
    error("Layout is only compatible with the GML format. To obtain a layout,
          please first convert your graph in GML using the tool convert-graph,
  an then send the result to display-graph.");
  }
  
  if (!$error) { 
  
    $graph = trim_text($graph);
    ## Load the parameters of the program in to an array

//     echo "<pre> 
//    $graph
//     </pre>";
    $parameters = array( 
      "request" => array(
      "informat"=>$in_format,
        "outformat"=>$out_format,
        "inputgraph"=>$graph,
        "ewidth"=>$ewidth,
        "scol"=>$s_col,
        "tcol"=>$t_col,
        "wcol"=>$w_col,
        "eccol"=>$ec_col,
        "sccol"=>$sc_col,
        "tccol"=>$tc_col,
        "layout"=>$layout,
        "tccol"=>$tc_col,
        "sccol"=>$sc_col,
        "eccol"=>$ec_col
        
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
    # Execute the command
    # Work with exception catch
    flush();
    try {
      $echoed = $client->display_graph($parameters);
      $soap_error = 0;
    } catch (Exception $soap_exception) {
      echo ("<pre>");
      echo "Error : \n",  $soap_exception->getMessage(), "\n";
      echo ("</pre>");
      $soap_error = 1;
    } 
    if (!$soap_error) {
      $response =  $echoed->response;
      $command = $response->command;
      //       echo "COMMAND $command";
      $client = $response->client;
      $server = $response->server;
      $resultURL = rsat_path_to_url($server);
#      $temp_file = explode('/',$server);
#      $temp_file = end($temp_file);
#      $resultURL = $WWW_RSA."/tmp/".$temp_file;
#      $URL['Display (png)'] = $resultURL;
      hourglass("off");
      
      # Display the results 
      echo "The results is available at the following URL ";
      echo "<a href = '$resultURL'>$resultURL</a><br>";
      echo "<a href = '$resultURL'><img src = '$resultURL' width = '50%'></a><br>Click on the icon to display the high-resolution image";
   }
 }

?>
