<html>
<head>
   <title>Network Analysis Tools - random-graph</title>
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
  $cmd_file = getTempFileName('commands_rdmgraph', '.txt');
  $cmd_handle = fopen($cmd_file, 'a');

  title('random-graph - results');
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

  $now = date("Ymd_His");
  $graph = $_REQUEST['graph'];

  $mean = $_REQUEST['mean'];
  $sd = $_REQUEST['stddev'];
  $nodes = $_REQUEST['nodes'];
  $edges = $_REQUEST['edges'];
  $self = $_REQUEST['self'];
  $normal = $_REQUEST['normal'];
  $duplicate = $_REQUEST['duplicate'];
  $col_conservation = $_REQUEST['col_conservation'];
  $no_single = $_REQUEST['no_single'];
  $directed = $_REQUEST['directed'];
  if ($duplicate == 'on') {
    $duplicate = 1;
  }
  if ($col_conservation == 'on') {
    $col_conservation = 1;
  }
  if ($self == 'on') {
    $self = 1;
  }
  if ($no_single == 'on') {
    $no_single = 1;
  }
  if ($directed == 'on') {
    $directed = 1;
  }
  if ($normal == 'on') {
    $normal = 1;
  }  
  
  $s_col = $_REQUEST['s_col'];
  $t_col = $_REQUEST['t_col'];
  $w_col = $_REQUEST['w_col'];
  $random_type = $_REQUEST['random_type'];
  
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
  ## If no graph are submitted and random_type = scratch -> error
  if ($graph == "" && $graph_file == "" &&  $random_type !=  'scratch') {
    $error = 1;
    error("You must submit an input graph with randomization type : $random_type");
  }
  ## no_single option is only valid for ER random type
  if ($random_type !=  'ER' && $no_single) {
    $error = 1;
    error("Option -no_single cannot be used with randomization type $random_type");
  }
  ## scratch and no edge or node requirement -> error
  if ($random_type ==  'scratch' && ($nodes == "" || $edges == "")) {
    $error = 1;
    error("You did not specify any node or edge number for the 'from scratch' randomization type");
  }
  ## If no graph are submitted and random_type != ER ou scratch -> error
  if (($random_type != 'ER' &&  $random_type != 'scratch')  && ($mean != "" && $sd != "")) {
    $error = 1;
    error("You must not specify randomization type $random_type	with -sd and -mean option");
  }  
  if (!$error) { 
  
    $graph = trim_text($graph);
    ## Load the parameters of the program in to an array
    $parameters = array( 
      "request" => array(
        "informat"=>$in_format,
        "outformat"=>$out_format,
        "inputgraph"=>$graph,
        "scol"=>$s_col,
        "tcol"=>$t_col,
        "wcol"=>$w_col,
        "directed"=>$directed,
        "normal"=>$normal,
        "mean"=>$mean,
        "sd"=>$sd,
        "no_single"=>$no_single,
        "duplicate"=>$duplicate,
        "nodes"=>$nodes,
        "edges"=>$edges,
        "random_type"=>$random_type,
	"self"=>$self
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
#    echo "<pre>";

    try {
      $echoed = $client->random_graph($parameters);
      $soap_error = 0;
    } catch (Exception $soap_exception) {
      echo ("<pre>");
      echo "Error : \n",  $soap_exception->getMessage(), "\n";
      echo ("</pre>");
      $soap_error = 1;
    } 



    
#    echo "</pre>";
    $response =  $echoed->response;
    $command = $response->command;
    $server = $response->server;
    $client = $response->client;
#    $temp_file = explode('/',$server);
#    $temp_file = end($temp_file);
#    $resultURL = $WWW_RSA."/tmp/".$temp_file;
    store_command($command, "random-graph", $cmd_handle);
    $URL['Random graph'] = rsat_path_to_url($server);
     
    $server = rtrim ($server);
    hourglass("off");



//     # Display the results
//     echo "The results is available at the following URL ";
//     echo "<a href = '$resultURL'>$resultURL</a>";
//     echo "<hr>";


    ## Close command handle
    fclose($cmd_handle);
    $URL['Server commands'] = rsat_path_to_url($cmd_file);

    ## DISPLAY THE RESULT
    print_url_table($URL);
    

    ## Display the "Next step" table
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
            <input type='hidden' NAME='wcol' VALUE='3'>";
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
        <FORM METHOD='POST' ACTION='convert_graph_form.php'>
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
          <INPUT type='submit' value='Convert $out_format to another format'>
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
          <input type='hidden' NAME='visant_directed' VALUE='$directed'>
          <input type='hidden' NAME='tab_java' VALUE='0'>";
          echo "
          <INPUT type='submit' value='Load in VisANT'>
          </form>
        </td>
      </tr>
  </table>";
 }

?>
