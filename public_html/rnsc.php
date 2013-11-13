<html>
<head>
   <title>Network Analysis Tools - RNSC</title>
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
  $cmd_file = getTempFileName('commands_rnsc', '.txt');
  $cmd_handle = fopen($cmd_file, 'a');

  title('RNSC - results');
  # Error status
  $error = 0;
  # Get parameters
  $in_format = $_REQUEST['in_format'];
  if ($_FILES['graph_file']['name'] != "") {
    $graph_file = uploadFile('graph_file');
  } else if ($_REQUEST['pipe_graph_file'] != "")  {
    $graph_file = $_REQUEST['pipe_graph_file'];
  }
  $now = date("Ymd_His");
  $graph = $_REQUEST['graph'];
  $s_col = $_REQUEST['s_col'];
  $t_col = $_REQUEST['t_col'];
  $max_clust = $_REQUEST['max_clust'];
  $tabulength = $_REQUEST['tabulength'];
  $tabulist = $_REQUEST['tabulist'];
  $naive_stop = $_REQUEST['naive_stop'];
  $scale_stop = $_REQUEST['scale_stop'];
  $div_freq = $_REQUEST['div_freq'];
  $shf_div_len = $_REQUEST['shf_div_len'];
  $exp_nb = $_REQUEST['exp_nb'];
  
  
  
  
  ## If a file and a graph are submitted -> error
  if ($graph != "" && $graph_file != "") {
    $error = 1;
    error("The graph can be entered either in the 'Graph' box, or uploaded, but you cannot fill both options simultaneously");
  }

  ## No specification of the source and target columns
  if ($in_format == "tab" && $s_col == "" && $t_col == "") {
    warning("Default value for source and target columns for tab-delimited input format are 1 and 2 respectively. ");
    info("The graph is considered as unweighted.");
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

    # Info message
    info("Results will appear below");
    echo"<hr>\n";
  hourglass("on");
    # Open the SOAP client
    $client = new SoapClient(
                       $neat_wsdl,
// "http://localhost/rsat/web_services/RSATWS2.wsdl",
                           array(
                                 'trace' => 1,
                                 'soap_version' => SOAP_1_1,
                                 'style' => SOAP_DOCUMENT,
                                 'encoding' => SOAP_LITERAL
                                 )
                           );



    ## Convert-graph
    $cg_parameters = array( 
      "request" => array(
        "inputgraph"=>$graph,
        "informat"=>$in_format,
        "outformat"=>"rnsc",
        "scol"=>$s_col,
        "tcol"=>$t_col
      )
    );
//     print_r ($cg_parameters);

    # Execute the command and catch the errors
    try {
      $cg_echoed = $client->convert_graph($cg_parameters);
      $soap_error = 0;
    } catch (Exception $soap_exception) {
      echo ("<pre>");
      echo "Error : \n\t",  $soap_exception->getMessage(), "\n";
      echo ("</pre>");
      $soap_error = 1;
      exit(1);
    }  
    $cg_response = $cg_echoed->response;
    $cg_command = $cg_response->command;
    store_command($cg_command, "graph conversion",$cmd_handle);
    $cg_server = $cg_response->server;
    $cg_client = $cg_response->client;
    $cg_server = rtrim ($cg_server);
    $cg_temp_file = explode('/',$cg_server);
    $cg_temp_file = end($cg_temp_file);
    $rnsc_graph = $cg_server.".rnsc";
    $rnsc_nodes = $cg_server."_node_names.rnsc";
    
    ## RNSC                           
    ## Load the parameters of the program in to an array
    $cg_graph = storeFile($rnsc_graph);
    $cg_nodes = storeFile($rnsc_nodes);
    $cg_graph = $cg_graph."\n";
    $rnsc_parameters = array( 
      "request" => array(
        "inputgraph"=>$cg_graph,
        "max_clust"=>$max_clust,
        "tabulength"=>$tabulength,
        "tabulist"=>$tabulist,
        "naive_stop"=>$aive_stop,
        "scale_stop"=>$scale_stop,
        "div_freq"=>$div_freq,
        "shf_div_len"=>$shf_div_len,
        "exp_nb"=>$exp_nb
      )
    );

    # Execute the command and catch the errors
    try {
      $rnsc_echoed = $client->rnsc($rnsc_parameters);
      $soap_error = 0;
    } catch (Exception $soap_exception) {
      echo ("<pre>");
      echo "Error : \n\t",  $soap_exception->getMessage(), "\n";
      echo ("</pre>");
      $soap_error = 1;
      exit(1);
    }  
    $rnsc_response = $rnsc_echoed->response;
    $rnsc_command = $rnsc_response->command;
    store_command($rnsc_command, "rnsc",$cmd_handle);
    $rnsc_server = $rnsc_response->server;
    $rnsc_client = $rnsc_response->client;
//     $rnsc_server = rtrim ($rnsc_server);
//     $rnsc_temp_file = explode('/',$rnsc_server);
//     $rnsc_temp_file = end($rnsc_temp_file);
    //    $URL['Clusters (RNSC format) old'] = $WWW_RSA."/tmp/".$rnsc_temp_file;
    $URL['Clusters (RNSC format)'] = rsat_path_to_url($rnsc_server);


    # Convert-classes 
    ## Load the parameters of the program in to an array
    $input_classes = storeFile($rnsc_server);
    
    $cc_parameters = array( 
      "request" => array(
        "inputclasses"=>$input_classes,
        "informat"=>"rnsc",
        "outformat"=>"tab",
        "names"=>$cg_nodes
      )
    );
    # Execute the command
    $cc_echoed = $client->convert_classes($cc_parameters);

    $cc_response = $cc_echoed->response;
    $cc_command = $cc_response->command;
    store_command($cc_command, "cluster formatting",$cmd_handle);
    $cc_server = $cc_response->server;
    $cc_client = $cc_response->client;
    //     $cc_server = rtrim ($cc_server);
    //     $cc_temp_file = explode('/',$cc_server);
    //     $cc_temp_file = end($cc_temp_file);
    $URL['Clusters (tab format)'] = rsat_path_to_url($cc_server);
    
    # contingency-table
    ## Load the parameters of the program into an array
    $cc_input_file = storeFile($cc_server);
    $ct_parameters = array(
       "request" => array(
        "inputfile" => $cc_input_file,
        "col1" => 2,
        "col2" => 1,
        "margin" => 1,
          )
      );
    $ct_echoed = $client->contingency_table($ct_parameters);
    $ct_response = $ct_echoed->response;
    $ct_command = $ct_response->command;
    store_command($ct_command, "cluster sizes",$cmd_handle);
    $ct_server = $ct_response->server;
    $ct_client = $ct_response->client;
    //     $ct_server = rtrim ($ct_server);
    //     $ct_temp_file = explode('/',$ct_server);
    //     $ct_temp_file = end($ct_temp_file);
    $URL['CLass/members table'] = rsat_path_to_url($ct_server);

   # Compute cluster size distribution with classfreq
    $cf_inputfile =  storeFile($ct_server);
//     $cf_echoed = $client->contingency_table($cf_parameters);
    $cf_parameters = array(
      "request" => array(
        "inputFile" => $cf_inputfile,
        "col" => 2,
        "classinterval"=>1
        )
      );    
    $cf_echoed = $client->classfreq($cf_parameters);
    $cf_response = $cf_echoed->response;
    $cf_command = $cf_response->command;
    store_command($cf_command, "cluster size distribution",$cmd_handle);
    $cf_server = $cf_response->server;
    $cf_client = $cf_response->client;
//     $cf_server = rtrim ($cf_server);
//     $cf_temp_file = explode('/',$cf_server);
//     $cf_temp_file = end($cf_temp_file);
    $URL['Cluster size distrib'] = rsat_path_to_url($cf_server);
    # XYgraph
    $xy_inputfile =  storeFile($cf_server);
    $xy_parameters = array( 
       "request" => array(
         "inputFile"=>$xy_inputfile,
         "xcol"=>"2",
         "ycol"=>"4",
         "format"=>"png",
         "lines"=>1,
         "xmin"=>0,
         "title1"=>"Cluster size distribution",
         "xleg1"=>"Cluster size",
         "yleg1"=>"Number of clusters",
       )
     );
    $xy_echoed = $client->xygraph($xy_parameters);
    $xy_response = $xy_echoed->response;
    $xy_command = $xy_response->command;
    store_command($xy_command, "cluster size distribution graph",$cmd_handle);
    $xy_server = $xy_response->server;
    $xy_client = $xy_response->client;
//     $xy_server = rtrim ($xy_server);
//     $xy_temp_file = explode('/',$xy_server);
//     $xy_temp_file = end($xy_temp_file);
    $xy_resultURL = rsat_path_to_url($xy_server);
    $URL['Cluster size distrib graph'] = $xy_resultURL;
    
    echo ("<a href = '$xy_resultURL'><img align = 'center' src='$xy_resultURL'></a><br>");
    echo "<br><hr>\n";
    
    # Input file
    $input_graph_file = writeTempFile("rnsc_input", $graph);
    hourglass("off");

    ## Close command handle
    fclose($cmd_handle);
//     $cmd_basename = explode('/',$cmd_file);
//     $cmd_basename = end($cmd_basename);
    $URL['Server commands'] = rsat_path_to_url($cmd_file);    

    ## DISPLAY THE RESULT
    print_url_table($URL);
     
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
          <input type='hidden' NAME='class_file' VALUE='$cc_server'>
          <INPUT type='submit' value='Compare these clusters with other clusters'>
        </form>
      </td>

      <TD>
        <FORM METHOD='POST' ACTION='graph_get_clusters_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'\>
          <input type='hidden' NAME='graph_file' VALUE='$input_graph_file'\>
          <input type='hidden' NAME='cluster_file' VALUE='$cc_server'\>
          <input type='hidden' NAME='graph_format' value = '$in_format'\>
          <input type='hidden' NAME='scol' VALUE='1'\>
          <input type='hidden' NAME='tcol' VALUE='2'\>
          
          <INPUT type='submit' value='Map those clusters on the network'\>
        </form>
      </td>
      <TD>
        <FORM METHOD='POST' ACTION='graph_cluster_membership_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'\>
          <input type='hidden' NAME='graph_file' VALUE='$input_graph_file'\>
          <input type='hidden' NAME='cluster_file' VALUE='$cc_server'\>
          <input type='hidden' NAME='graph_format' value = '$in_format'\>
          <input type='hidden' NAME='scol' VALUE='1'\>
          <input type='hidden' NAME='tcol' VALUE='2'\>
          
          <INPUT type='submit' value='Cluster membership'\>
        </form>
      </td> 
   ";  
  }
  
  
  
  
?>
