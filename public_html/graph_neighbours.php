<html>
<head>
   <title>NeA-tools - graph-neighbours</title>
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
  title('graph-neighbours - results');
  # Error status
  $error = 0;
  # Get parameters
  $in_format = $_REQUEST['in_format'];
  if ($_FILES['graph_file']['name'] != "") {
    $graph_file = uploadFile('graph_file');
  } else if ($_REQUEST['pipe_graph_file'] != "")  {
    $graph_file = $_REQUEST['pipe_graph_file'];
  }
  if ($_FILES['seeds_file']['name'] != "") {
    $seeds_file = uploadFile('seeds_file');
  }
  $now = date("Ymd_His");
  $graph = $_REQUEST['graph'];
  $seeds = $_REQUEST['seeds'];
  $s_col = $_REQUEST['s_col'];
  $t_col = $_REQUEST['t_col'];
  $w_col = $_REQUEST['w_col'];
  $stats = $_REQUEST['stats'];
  $steps = $_REQUEST['steps'];
  $all = $_REQUEST['allseeds'];
  $self = $_REQUEST['self'];
  if ($all == 'all') {
    $all = 1;
  } else {
    $all = 0;
  }
  if ($self == 'on') {
    $self = 1;
  } else {
    $self = 0;
  }
    
  ## If a file and a graph are submitted -> error
  if ($graph != "" && $graph_file != "") {
    $error = 1;
    error("You must not submit both a graph and a graph file");
  }
  ## If seeds and a node file are submitted -> error
  if ($seeds != "" && $seeds_file != "") {
    $error = 1;
    error("You must not submit both seeds and a seeds file");
  }
  ## If stats and self -> error
  if ($stats == 1 && $self == 1) {
    $error = 1;
    error("The 'one line per seed format' cannot be computed with seed node inclusion");
  }
  ## No specification of the source and target columns
  if ($in_format == "tab" && $s_col == "" && $t_col == "") {
    warning("Default value for source and target columns for tab-delimited input format are 1 and 2 respectively");
  }
  ## No specification of the weight column and stats output -> error
  if ($stats && $w_col == "") {
    error("One line per seed node output is only possible for weighted graph");
    $error = 1;
  }
  ## No specification of the weight column and stats output -> error
  if ($stats && $steps > 1) {
    error("One line per seed node output is only possible for distance 1 from seed nodes");
    $error = 1;
  }
  ## put the content of the file $graph_file in $graph
  if ($graph_file != "" && $graph == "") {
    $graph = storeFile($graph_file);
  }
  ## put the content of the file $seeds_file in $seeds
  if ($seeds_file != "" && $seeds == "") {
    $seeds = storeFile($seeds_file);
  }
  ## If no graph are submitted -> error
  if ($graph == "" && $graph_file == "") {
    $error = 1;
    error("You must submit an input graph");
  }

  if (!$error) { 
    $graph = trim_text($graph);
    $seeds = trim_text($seeds);
     
    ## Load the parameters of the program in to an array
    $gn_parameters = array( 
      "request" => array(
        "informat"=>$in_format,
        "seedfile"=>$seeds,
        "inputgraph"=>$graph,
        "scol"=>$s_col,
        "tcol"=>$t_col,
        "wcol"=>$w_col,
        "all"=>$all,
        "stats"=>$stats,
        "self"=>$self,
        "steps"=>$steps
      )
    );
    # Info message
    info("Results will appear below");
    echo"<hr>\n";
    hourglass("on");
  
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
    $gn_echoed = $soap_client->graph_neighbours($gn_parameters);

    $gn_response =  $gn_echoed->response;
    $gn_command = $gn_response->command;
    $gn_server = $gn_response->server;
    $gn_client = $gn_response->client;
    echo "</pre>";
    $gn_server = rtrim ($gn_server);
    $gn_temp_file = explode('/',$gn_server);
    $gn_temp_file = end($gn_temp_file);
    $gn_resultURL = $WWW_RSA."/tmp/".$gn_temp_file;
    # Text-to-html
    $gn_file = storeFile($gn_server);
    $tth_parameters = array( 
      "request" => array(
        "inputfile"=>$gn_file,
        "chunk"=>1000,
      )
    );
    
    $tth_echoed = $soap_client->text_to_html($tth_parameters);

    $tth_response =  $tth_echoed->response;
    $tth_command = $tth_response->command;
    $tth_server = $tth_response->server;
    $tth_client = $tth_response->client;
    echo "</pre>";
    hourglass("off");
    $tth_server = rtrim ($tth_server);
    $tth_temp_file = explode('/',$tth_server);
    $tth_temp_file = end($tth_temp_file);
    $tth_resultURL = $WWW_RSA."/tmp/".$tth_temp_file;    
    
    
    
    
    
    
    # Display the results
    echo "The results is available as text file at the following URL ";
    echo "<a href = '$gn_resultURL'>$gn_resultURL</a><br>"; 
    echo "The results is available as HTML page at the following URL ";
    echo "<a href = '$tth_resultURL'>$tth_resultURL</a><br>"; 
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
        <FORM METHOD='POST' ACTION='compare_classes_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='class_file' VALUE='$gn_server'>
          <INPUT type='submit' value='Compare the groups of neighbours'>
        </form>
      </td> 
   ";  
  }
?>
