<html>
<head>
   <title>Network Analysis Tools - contingency-stats</title>
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
  $cmd_file = getTempFileName('commands_conting_stats', '.txt');
  $cmd_handle = fopen($cmd_file, 'a');

  title('contingency-stats - results');
  # Error status
  $error = 0;
  # Get parameters
  if ($_FILES['matrix_file']['name'] != "") {
    $matrix_file = uploadFile('matrix_file');
  } else if ($_REQUEST['pipe_matrix_file'] != "")  {
    $matrix_file = $_REQUEST['pipe_matrix_file'];
  }
  $now = date("Ymd_His");
  $matrix = $_REQUEST['matrix'];
  $row_stats = $_REQUEST['rowstats'];
  $col_stats = $_REQUEST['colstats'];
  $stats = $_REQUEST['stats'];
  $tables = $_REQUEST['tables'];
  $margins = $_REQUEST['margins'];
  $return = array();
 
  if ($row_stats == 'on') {
    array_push($return, 'rowstats');
  }
  if ($col_stats == 'on') {
    array_push($return,'colstats');
  }
  if ($stats == 'on') {
    array_push($return,'stats');
  }
  if ($margins == 'on') {
    array_push($return,'margins');
  }
  if ($tables == 'on') {
    array_push($return,'tables');
  }
  $return = join(",", $return);

  ## If both a file upload and a matrix are submitted -> error
  if ($matrix != "" && $matrix_file != "") {
    $error = 1;
    error("You must not submit both a table and a table file");
  }

  if ($matrix_file != "" && $matrix == "") {
    $matrix = storeFile($matrix_file);
  }

  ## If no matrix is submitted -> error
  if ($matrix == "" && $matrix_file == "") {
    $error = 1;
    error("You must submit a contingency table");
  }
  
  if (!$error) { 
  
    $matrix = trim_text($matrix);

    ## Load the parameters of the program in to an array
    $parameters = array( 
      "request" => array(
        "inputfile"=>$matrix,
        "decimals"=>4,
        "return"=>$return
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
    $echoed = $soap_client->contingency_stats($parameters);
    $response =  $echoed->response;
    $command = $response->command;
    $server = $response->server;

#    $client = $response->client;
    store_command($command, "graph comparison", $cmd_handle);
    $URL['Contingency stats (tab)'] = rsat_path_to_url($server);
    hourglass("off");
    ## Text-to-html
    $server = rtrim ($server);
    $file = storeFile($server);
    $tth_parameters = array( 
      "request" => array(
        "inputfile"=>$file,
        "chunk"=>1000,
      )
    );
    $tth_echoed = $soap_client->text_to_html($tth_parameters);
    $tth_response =  $tth_echoed->response;
    $tth_command = $tth_response->command;
    $tth_server = $tth_response->server;
    $tth_client = $tth_response->client;
    store_command($tth_command, "text-to-html", $cmd_handle);
    $URL['Contingency stats (html)'] = rsat_path_to_url($tth_server);

    ## Close command handle
    fclose($cmd_handle);
    $URL['Server commands'] = rsat_path_to_url($cmd_file);

    ## DISPLAY THE RESULT
    print_url_table($URL);
    
  }
?>
