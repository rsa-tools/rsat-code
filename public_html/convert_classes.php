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
  $cmd_file = getTempFileName('commands_cclasses');
  $cmd_handle = fopen($cmd_file, 'a');

  title('convert-classes - results');
  # Error status
  $error = 0;
  # Get parameters
  $in_format = $_REQUEST['in_format'];
  $out_format = $_REQUEST['out_format'];
  if ($_FILES['classes_file']['name'] != "") {
    $classes_file = uploadFile('classes_file');
  } 
  $now = date("Ymd_His");
  $classes = $_REQUEST['classes'];

  $s_col = $_REQUEST['s_col'];
  $m_col = $_REQUEST['m_col'];
  $w_col = $_REQUEST['w_col'];
  
  ## If a file and a graph are submitted -> error
  if ($classes != "" && $classes_file != "") {
    $error = 1;
    error("You must not submit both a set of classes and a classes file");
  }

  ## No specification of the source and target columns
  if ($in_format == "tab" && $m_col == "" && $c_col == "") {
    warning("Default value for member and classes columns for tab-delimited input format are 1 and 2 respectively");
  }
  ## put the content of the file classes_file in $classes
  if ($classes_file != "" && $classes == "") {
    $classes = storeFile($classes_file);
  }
  ## If no graph are submitted -> error
  if ($classes == "" && $classes_file == "") {
    $error = 1;
    error("You must submit input classes");
  }

  
  if (!$error) { 
//       echo ("<pre>$classes</pre>");
    $classes = trim_text($classes);

    

    ## Load the parameters of the program in to an array
    $parameters= array( 
      "request" => array(
        "inputclasses"=>$classes,
        "informat"=>$in_format,
        "outformat"=>$out_format,
        "member_col"=>$m_col,
        "class_col"=>$c_col,
        "score_col"=>$w_col
      )
    );
    # Info message
    info("Results will appear below");
    echo"<hr>\n";
    hourglass("on");
    # Open the SOAP client
    $client = new SoapClient(
                       $neat_wsdl,
// "http://rsat.ulb.ac.be/rsat/web_services/RSATWS2.wsdl",
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
      $echoed = $client->convert_classes($parameters);
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
       echo ($command);
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
        <FORM METHOD='POST' ACTION='graph_get_clusters_form.php'>
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
          <INPUT type='submit' value='Map clusters or extract a subnetwork'>
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
          <input type='hidden' NAME='visant_directed' VALUE='$directed'>";
          echo "
          <INPUT type='submit' value='Load in VisANT'>
          </form>
        </td>
      </tr>

    </table>";
    }
  }
?>
