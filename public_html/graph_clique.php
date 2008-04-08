<html>
<head>
   <title>NeA-tools - graph-clique</title>
   <link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
</head>
<body class="results">
<?php 
  require ('functions.php');
  # log file update
  UpdateLogFile("neat","","");
  title('graph-clique - results');
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
  $s_col = $_REQUEST['s_col'];
  $t_col = $_REQUEST['t_col'];
  $min_size = $_REQUEST['min_size'];
  $max_size = $_REQUEST['max_size'];
    
  ## If a file and a graph are submitted -> error
  if ($graph != "" && $graph_file != "") {
    $error = 1;
    error("You must not submit both a graph and a graph file");
  }
  ## No specification of the source and target columns
  if ($in_format == "tab" && $s_col == "" && $t_col == "") {
    warning("Default value for source and target columns for tab-delimited input format are 1 and 2 respectively");
  }
  if (preg_match("/\d/", $thr)) {
    $error = 1;
    error("$min_size is not a valid minimal threshold for the size of clique");
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

  if (!$error) { 
    $graph = trim_text($graph);
     
    ## Load the parameters of the program in to an array
    $gc_parameters = array( 
      "request" => array(
        "informat"=>$in_format,
        "inputgraph"=>$graph,
        "scol"=>$s_col,
        "tcol"=>$t_col,
        "min_size"=>$min_size,
        "max_size"=>$max_size
        
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
    $gc_echoed = $soap_client->graph_clique($gc_parameters);

    $gc_response = $gc_echoed->response;
    $gc_command = $gc_response->command;
    $gc_server = $gc_response->server;
    $gc_client = $gc_response->client;
    echo "</pre>";
    $gc_server = rtrim ($gc_server);
    $gc_temp_file = explode('/',$gc_server);
    $gc_temp_file = end($gc_temp_file);
    $gc_resultURL = $WWW_RSA."/tmp/".$gc_temp_file;
    # Text-to-html
    $gc_file = storeFile($gc_server);
    $tth_parameters = array( 
      "request" => array(
        "inputfile"=>$gc_file,
        "chunk"=>1000,
      )
    );
    
    $tth_echoed = $soap_client->text_to_html($tth_parameters);

    $tth_response =  $tth_echoed->response;
    $tth_command = $tth_response->command;
    $tth_server = $tth_response->server;
    $tth_client = $tth_response->client;
    echo "</pre>";
    $tth_server = rtrim ($tth_server);
    $tth_temp_file = explode('/',$tth_server);
    $tth_temp_file = end($tth_temp_file);
    $tth_resultURL = $WWW_RSA."/tmp/".$tth_temp_file;    
    
    
    
    
    
    
    # Display the results
    echo "The results is available as text file at the following URL ";
    echo "<a href = '$gc_resultURL'>$gc_resultURL</a><br>"; 
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
          <input type='hidden' NAME='class_file' VALUE='$gc_server'>
          <INPUT type='submit' value='Compare the cliques'>
        </form>
      </td> 
   ";  
  }
?>
