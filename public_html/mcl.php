<html>
<head>
   <title>NeA-tools - MCL</title>
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
  title('MCL - results');
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
  $w_col = $_REQUEST['w_col'];
  $inflation = $_REQUEST['inflation'];
  
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
//                        $neat_wsdl,
"http://rsat.scmbb.ulb.ac.be/rsat/web_services/RSATWS2.wsdl",
                           array(
                                 'trace' => 1,
                                 'soap_version' => SOAP_1_1,
                                 'style' => SOAP_DOCUMENT,
                                 'encoding' => SOAP_LITERAL
                                 )
                           );
                           
    ## Convert-graph
    $cg_parameters =  array( 
      "request" => array(
        "inputgraph"=>$graph,
        "informat"=>$in_format,
        "outformat"=>"tab",
        "scol"=>$s_col,
        "tcol"=>$t_col,
        "wcol"=>$w_col
      )
    );
    # Execute the command 
    $cg_echoed = $client->convert_graph($cg_parameters);

    $cg_response = $cg_echoed->response;
    $cg_command = $cg_response->command;
//     echo ("$cg_command");
    $cg_server = $cg_response->server;
    $cg_client = $cg_response->client;
    $cg_server = rtrim ($cg_server);
    $cg_temp_file = explode('/',$cg_server);
    $cg_temp_file = end($cg_temp_file);
    $cg_resultURL = $WWW_RSA."/tmp/".$cg_temp_file;    
    
    ## MCL                           
    ## Load the parameters of the program in to an array
    $cg_graph = storeFile($cg_server);
    $mcl_parameters = array( 
      "request" => array(
        "inputgraph"=>$cg_graph,
        "inflation"=>$inflation
      )
    );
    # Execute the command
    $mcl_echoed = $client->mcl($mcl_parameters);

    $mcl_response = $mcl_echoed->response;
    $mcl_command = $mcl_response->command;
//     echo ("$mcl_command");
    $mcl_server = $mcl_response->server;
    $mcl_client = $mcl_response->client;
    $mcl_server = rtrim ($mcl_server);
    $mcl_temp_file = explode('/',$mcl_server);
    $mcl_temp_file = end($mcl_temp_file);
    $mcl_resultURL = $WWW_RSA."/tmp/".$mcl_temp_file;
    # Convert-classes 
    ## Load the parameters of the program in to an array
    $input_classes = storeFile($mcl_server);
    
    $cc_parameters = array( 
      "request" => array(
        "inputclasses"=>$input_classes,
        "informat"=>"mcl",
        "outformat"=>"tab"
      )
    );
    # Execute the command
    $cc_echoed = $client->convert_classes($cc_parameters);

    $cc_response = $cc_echoed->response;
    $cc_command = $cc_response->command;
    $cc_server = $cc_response->server;
    $cc_client = $cc_response->client;
    $cc_server = rtrim ($cc_server);
    $cc_temp_file = explode('/',$cc_server);
    $cc_temp_file = end($cc_temp_file);
    $cc_resultURL = $WWW_RSA."/tmp/".$cc_temp_file;    
    
    
    # contingency-table
    ## Load the parameters of the program into an array
    $cc_input_file = storeFile($cc_server);
    $ct_parameters = array(
      "request" => array(
        "inputfile" => $cc_input_file,
        "col1" => 2,
        "col2" => 1,
        )
      );
    echo"<pre>";
    $ct_echoed = $client->contingency_table($ct_parameters);
    $ct_response = $ct_echoed->response;
    $ct_command = $ct_response->command;
    $ct_server = $ct_response->server;
    $ct_client = $ct_response->client;
    $ct_server = rtrim ($ct_server);
    $ct_temp_file = explode('/',$ct_server);
    $ct_temp_file = end($ct_temp_file);
    $ct_resultURL = $WWW_RSA."/tmp/".$ct_temp_file;
    echo "</pre>";
    # classfreq 
    $cf_inputfile =  storeFile($ct_server);
    $cf_echoed = $client->contingency_table($cf_parameters);
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
    $cf_server = $cf_response->server;
    $cf_client = $cf_response->client;
    $cf_server = rtrim ($cf_server);
    $cf_temp_file = explode('/',$cf_server);
    $cf_temp_file = end($cf_temp_file);
    $cf_resultURL = $WWW_RSA."/tmp/".$cf_temp_file;    
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
    $xy_server = $xy_response->server;
    $xy_client = $xy_response->client;
    $xy_server = rtrim ($xy_server);
    $xy_temp_file = explode('/',$xy_server);
    $xy_temp_file = end($xy_temp_file);
    $xy_resultURL = $WWW_RSA."/tmp/".$xy_temp_file;
    
    hourglass("off");
    echo ("<a href = '$xy_resultURL'><img align = 'center' src='$xy_resultURL' width ='50%'></a><br>");
    echo "<br><hr>\n";
    
    
    # Display the results
    echo "The results is available at the following URL ";
    echo "<a href = '$cc_resultURL'>$cc_resultURL</a>"; 
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
          <input type='hidden' NAME='class_file' VALUE='$cc_server'>
          <INPUT type='submit' value='Compare these clusters with other clusters'>
        </form>
      </td>

      <TD>
        <FORM METHOD='POST' ACTION='graph_get_clusters_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'\>
          <input type='hidden' NAME='graph_file' VALUE='$cg_server'\>
          <input type='hidden' NAME='cluster_file' VALUE='$cc_server'\>
          <input type='hidden' NAME='graph_format' value = 'tab'\>
          <input type='hidden' NAME='scol' VALUE='1'\>
          <input type='hidden' NAME='tcol' VALUE='2'\>
          <input type='hidden' NAME='wcol' VALUE='3'\>
          
          <INPUT type='submit' value='Map those clusters on the network'\>
        </form>
      </td>
      <TD>
        <FORM METHOD='POST' ACTION='graph_cluster_membership_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'\>
          <input type='hidden' NAME='graph_file' VALUE='$cg_server'\>
          <input type='hidden' NAME='cluster_file' VALUE='$cc_server'\>
          <input type='hidden' NAME='graph_format' value = 'tab'\>
          <input type='hidden' NAME='scol' VALUE='1'\>
          <input type='hidden' NAME='tcol' VALUE='2'\>
          <input type='hidden' NAME='wcol' VALUE='3'\>
          
          <INPUT type='submit' value='Cluster membership'\>
        </form>
      </td> 
   ";  
  }
  
  
  
  
?>
