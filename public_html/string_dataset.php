<html>
<head>
   <title>NeA-tools - convert-graph</title>
   <link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
</head>
<body class="results">
<?php 
  require ('functions.php');
  # log file update
  UpdateLogFile("neat","","");
  title('string dataset download- results');
  # Error status
  $error = 0;
  # Get parameters
  $now = date("Ymd_His");
  $genes = $_REQUEST['genes'];
  $genes = trim_text($genes);
  $organism = $_REQUEST['organism'];
  $channels = array();
  $file_content_tab = "";
  if ($_REQUEST['automated_textmining'] != "") {
    array_push ($channels, 'automated_textmining');
  }
  if ($_REQUEST['experimental_interaction_data'] != "") {
    array_push ($channels, 'experimental_interaction_data');
  }
  if ($_REQUEST['gene_cooccurence'] != "") {
    array_push ($channels, 'gene_cooccurence');
  }
  if ($_REQUEST['gene_fusion_events'] != "") {
    array_push ($channels, 'gene_fusion_events');
  }
  if ($_REQUEST['gene_coexpression'] != "") {
    array_push ($channels, 'gene_coexpression');
  }
  if ($_REQUEST['combined_confidence'] != "") {
    array_push ($channels, 'combined_confidence');
  }
  

  
  # WGET COMMAND
  $gene_list = explode("\n", $genes); 
  echo "<pre>";
  print_r ($genes);
  print_r ($gene_list);
  echo "</pre>";
  for ($i = 0; $i < count($gene_list); $i++) {
    $gene_list[$i] = rtrim ($gene_list[$i]);
    if ($gene_list[$i] == "") {
      next;
    }
    $temp_file;
    $temp_command = "mktemp tmp/string.XXXXXX";
    exec($temp_command, $temp_file);
    trim($temp_file[0]);
    $temp_file = $temp_file[0];
    $wget_command = "wget 'http://stitch.embl.de/api/psi-mi/interactions?identifier=$gene_list[$i]&species=$organism' -O $temp_file";
    exec($wget_command);
    $file_content_xml = file_get_contents ($temp_file);
    ## Load the parameters of the program in to an array
    $parameters = array( 
      "request" => array(
        "inputfile"=>$file_content_xml,
        "channels"=>implode(",", $channels)
      )
    );

    # Open the SOAP client
    $client = new SoapClient(
//                        $neat_wsdl,
"http://rsat.scmbb.ulb.ac.be/rsat/web_services/RSATWS-test.wsdl",
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
    $file_content_tab_add = storeFile($server);
    $file_content_tab = $file_content_tab.$file_content_tab_add;
  }
  
  $temp_file = writeTempFile("string", $file_content_tab);
     $resultURL = $WWW_RSA.$temp_file;
      # Info message
      info("Results will appear below");
      echo"<hr>\n";
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
        <FORM METHOD='POST' ACTION='graph_node_degree_form.php'>
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
          <INPUT type='submit' value='Nodes degrees computation'>
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
?>
