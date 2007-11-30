<html>
<head>
   <title>NeA-tools - contingency-stats</title>
   <link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
</head>
<body class="results">
<?php 
  require ('functions.php');
  title('contigency-stats - results');
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
  ## If a file and a graph are submitted -> error
  if ($matrix != "" && $matrix_file != "") {
    $error = 1;
    error("You must not submit both a table and a table file");
  }

  if ($matrix_file != "" && $matrix == "") {
    $matrix = storeFile($matrix_file);
  }
  ## If no graph are submitted -> error
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
    echo ("<pre>");
    $echoed = $client->contingency_stats($parameters);

    $response =  $echoed->response;
    $command = $response->command;
    $server = $response->server;
    $client = $response->client;
    $server = rtrim ($server);
    $temp_file = explode('/',$server);
    $temp_file = end($temp_file);
    $resultURL = $WWW_RSA."/tmp/".$temp_file;
    echo ("</pre>");
    # Display the results
    echo "The results is available at the following URL ";
    echo "<a href = '$resultURL'>$resultURL</a>"; 
    echo "<hr>\n";
    
  }
?>
