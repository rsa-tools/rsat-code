<html>
<head>
   <title>Network Analysis Tools - convert-graph</title>
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
  $cmd_file = getTempFileName('drawheatmap_cmd');
  $cmd_handle = fopen($cmd_file, 'a');

  title('draw-heatmap - results');
  # Error status
  $error = 0;
  # Get parameters
  $out_format = $_REQUEST['out_format'];
  if ($_FILES['table_file']['name'] != "") {
    $classes_file = uploadFile('table_file');
  } 
  $now = date("Ymd_His");
  $classes = $_REQUEST['table'];

  $row_h = $_REQUEST['rowh'];
  $col_w = $_REQUEST['colw'];
  $min = $_REQUEST['min'];
  $max = $_REQUEST['max'];
  $gradient = $_REQUEST['gradient'];
  $no_text = $_REQUEST['no_text'];
  $row_names = $_REQUEST['rownames'];
  
  ## If a file and a graph are submitted -> error
  if ($table != "" && $table_file != "") {
    $error = 1;
    error("You must not submit both table and a table file");
  }


  ## put the content of the file classes_file in $classes
  if ($table_file != "" && $table == "") {
    $table = storeFile($table_file);
  }
  ## If no graph are submitted -> error
  if ($classes == "" && $classes_file == "") {
    $error = 1;
    error("You must submit a table");
  }

  
  if (!$error) { 
//       echo ("<pre>$classes</pre>");
    $table = trim_text($table);

    

    ## Load the parameters of the program in to an array
    $parameters= array( 
      "request" => array(
        "inputfile"=>$table,
        "outformat"=>$out_format,
        "col_width"=>$col_w,
        "row_height"=>$row_h,
        "min"=>$min,
        "max"=>$max,
        "gradient"=>$gradient,
        "no_text"=>$no_text,
        "row_names"=>$row_names
        
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
//     echo ("<pre>");
//     $echoed = $client->convert_graph($parameters);
//     echo ("</pre>");
    # Work with exception catch
    try {
      $echoed = $client->draw_heatmap($parameters);
      $soap_error = 0;
    } catch (Exception $soap_exception) {
      echo ("<pre>");
      echo "Error : \n\t",  $soap_exception->getMessage(), "\n";
      echo ("</pre>");
      $soap_error = 1;
    }  
    if (!$soap_error) {
      $response =  $echoed->response;
//       print_r ($response);
      $command = $response->command;
#      echo("<p><b>Convert-graph command:</b> $command</p>\n");
      $server = $response->server;
      $client = $response->client;
#      $server = rtrim ($server);
#      $temp_file = explode('/',$server);
#      $temp_file = end($temp_file);
       store_command($command, "class conversion", $cmd_handle);
       $URL['Result'] = rsat_path_to_url($server);

      hourglass("off");


    ## Close command handle
    fclose($cmd_handle);
    $URL['Server commands'] = rsat_path_to_url($cmd_file);

    ## DISPLAY THE RESULT
    print_url_table($URL);

      ## Send result to next step
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
          <input type='hidden' NAME='class_file' VALUE='$server'>
          <INPUT type='submit' value='Compare these clusters with other clusters'>
        </form>
      </td>
</tr>

   ";  
  }}
?>
