<html>
<head>
   <title>Network Analysis Tools - convert-graph</title>
   <link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
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
  title('string dataset download - results');
  # Error status
  $error = 0;
  # Get parameters
  $now = date("Ymd_His");
  $gene_list = $_REQUEST['genes'];
  $gene_list = explode("-sep-", $gene_list);
  $genes = array();
  for ($i = 0; $i <= count($gene_list); $i++) {
    $gene = $gene_list[$i];
    $gene = trim($gene);
    if ($_REQUEST[$gene] == 1) {
      $gene = str_replace("_", ".", $gene);
      array_push ($genes, $gene);
    } 
  }
  $organism = $_REQUEST['organism'];
  $channels = $_REQUEST['channels'];
  $uth = $_REQUEST['uth'];
  $lth = $_REQUEST['lth'];
  if ($uth == "none") {
    $uth = "";
  } else if (!is_double($uth) || $uth > 1) { 
    warning("Upper threshold value is not a valid value. Threshold $uth will be ignored.");
    $uth = "";
  } 
  if ($lth == "none") {
    $lth = "";
  } else if (!is_double($lth) || $lth > 1) { 
    warning("Upper threshold value is not a valid value. Threshold $lth will be ignored.");
    $lth = "";
  } 
  $file_content_tab = "";
  # Info message
  info("Results will appear below");
  echo"<hr>\n";
  hourglass("on");
  # WGET COMMAND
  for ($i = 0; $i < count($genes); $i++) {
    $genes[$i] = rtrim ($genes[$i]);
    if ($genes[$i] == "") {
      continue;
    }
    $temp_file = writeTempFile("String", "");
    $wget_command = "wget 'http://stitch.embl.de/api/psi-mi/interactions?identifier=$genes[$i]&species=$organism' -O $temp_file";
    exec($wget_command);
    $file_content_xml = file_get_contents ($temp_file);
    ## Load the parameters of the program in to an array
    $parameters = array( 
      "request" => array(
        "inputfile"=>$file_content_xml,
        "channels"=>$channels,
        "uth"=>$uth,
        "lth"=>$lth,
        "interactor_type"=>'protein'
      )
    );

    # Open the SOAP client
    $client = new SoapClient(
                       $neat_wsdl,
// "http://rsat.ulb.ac.be/rsat/web_services/RSATWS.wsdl",
                           array(
                                 'trace' => 1,
                                 'soap_version' => SOAP_1_1,
                                 'style' => SOAP_DOCUMENT,
                                 'encoding' => SOAP_LITERAL
                                 )
                           );
    try {
      $echoed = $client->parse_psi_xml($parameters);
      $soap_error = 0;
    } catch (Exception $soap_exception) {
      echo ("<pre>");
      echo "Error : \n",  $soap_exception->getMessage(), "\n";
      echo ("</pre>");
      $soap_error = 1;
    }
    $response = $echoed->response;
    $server = $response->server;
    $server = rtrim($server);
    $command = $response->command;
    $file_content_tab_add = storeFile($server);
    $file_content_tab = $file_content_tab.$file_content_tab_add;
  }
  
  $temp_file = writeTempFile("string", $file_content_tab);
     $resultURL = $WWW_RSA.$temp_file;
hourglass("off");
      # Display the results
      echo "The results is available at the following URL ";
      echo "<a href = '$resultURL'>$resultURL</a>"; 
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
          <FORM METHOD='POST' ACTION='display_graph_form.php'>
            <input type='hidden' NAME='pipe' VALUE='1'>
            <input type='hidden' NAME='graph_file' VALUE='$temp_file'>
            <input type='hidden' NAME='graph_format' VALUE='tab'>

              <input type='hidden' NAME='scol' VALUE='1'>
              <input type='hidden' NAME='tcol' VALUE='2'>
              <input type='hidden' NAME='wcol' VALUE='3'>
              <input type='hidden' NAME='eccol' VALUE='4'>
            <INPUT type='submit' value='Display the graph'>
          </form>
        </td>
        <TD>
          <FORM METHOD='POST' ACTION='compare_graphs_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$temp_file'>
          <input type='hidden' NAME='graph_format' VALUE='tab'>
            <input type='hidden' NAME='scol' VALUE='1'>
            <input type='hidden' NAME='tcol' VALUE='2'>
            <input type='hidden' NAME='wcol' VALUE='3'>

          <INPUT type='submit' value='Compare this graph to another one'>
        </form>
      </td>
      <TD>
        <FORM METHOD='POST' ACTION='random_graph_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$temp_file'>
          <input type='hidden' NAME='graph_format' VALUE='tab'>
            <input type='hidden' NAME='scol' VALUE='1'>
            <input type='hidden' NAME='tcol' VALUE='2'>
            <input type='hidden' NAME='wcol' VALUE='3'>
          <INPUT type='submit' value='Randomize this graph'>
        </form>
      </td>
    </tr>
    <tr>
      <TD>
        <FORM METHOD='POST' ACTION='graph_get_clusters_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$temp_file'>
          <input type='hidden' NAME='graph_format' VALUE='tab'>
            <input type='hidden' NAME='scol' VALUE='1'>
            <input type='hidden' NAME='tcol' VALUE='2'>
            <input type='hidden' NAME='wcol' VALUE='3'>
          <INPUT type='submit' value='Map clusters or extract a subnetwork'>
        </form>
      </td>
      <TD>
        <FORM METHOD='POST' ACTION='graph_topology_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$temp_file'>
          <input type='hidden' NAME='graph_format' VALUE='tab'>
            <input type='hidden' NAME='scol' VALUE='1'>
            <input type='hidden' NAME='tcol' VALUE='2'>
            <input type='hidden' NAME='wcol' VALUE='3'>
          <INPUT type='submit' value='Nodes topology statistics'>
        </form>
      </td>
      <TD>
        <FORM METHOD='POST' ACTION='graph_neighbours_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$temp_file'>
          <input type='hidden' NAME='graph_format' VALUE='tab'>
            <input type='hidden' NAME='scol' VALUE='1'>
            <input type='hidden' NAME='tcol' VALUE='2'>
            <input type='hidden' NAME='wcol' VALUE='3'>
          <INPUT type='submit' value='Neighbourhood analysis'>
        </form>
      </td>    
    </tr>
    <TR>
      <TD>
        <FORM METHOD='POST' ACTION='mcl_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$temp_file'>
          <input type='hidden' NAME='graph_format' VALUE='tab'>
            <input type='hidden' NAME='scol' VALUE='1'>
            <input type='hidden' NAME='tcol' VALUE='2'>
            <input type='hidden' NAME='wcol' VALUE='3'>
          <INPUT type='submit' value='MCL Graph clustering'>
        </form>
      </td>
      <TD>
        <FORM METHOD='POST' ACTION='alter_graph_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$temp_file'>
          <input type='hidden' NAME='graph_format' VALUE='tab'>
            <input type='hidden' NAME='scol' VALUE='1'>
            <input type='hidden' NAME='tcol' VALUE='2'>
            <input type='hidden' NAME='wcol' VALUE='3'>
          <INPUT type='submit' value='Graph alteration'>
        </form>
      </td>
      <TD>
        <FORM METHOD='POST' ACTION='pathfinder_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$temp_file'>
          <input type='hidden' NAME='in_format' VALUE='tab'>";
          echo "
          <INPUT type='submit' value='Path Finding'>
        </form>
      </td>
      </tr>
        <tr>
        <TD>
          <FORM METHOD='POST' ACTION='visant.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='visant_graph_file' VALUE='$temp_file'>
          <input type='hidden' NAME='visant_graph_format' VALUE='tab'>
          <input type='hidden' NAME='visant_directed' VALUE='$directed'>";
          echo "
          <INPUT type='submit' value='Load in VisANT'>
          </form>
        </td>
      </tr>

    </table>";
?>
