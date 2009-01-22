<html>
<head>
   <title>NeA-tools - graph-cluster-membership</title>
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
  title('graph-cluster-membership - results');
  # Error status
  $error = 0;
  # Get parameters
  $in_format = $_REQUEST['in_format'];
  if ($_FILES['graph_file']['name'] != "") {
    $graph_file = uploadFile('graph_file');
  } else if ($_REQUEST['pipe_graph_file'] != "")  {
    $graph_file = $_REQUEST['pipe_graph_file'];
  }
  if ($_FILES['clusters_file']['name'] != "") {
    $clusters_file = uploadFile('clusters_file');
  } 
  $now = date("Ymd_His");
  $graph = $_REQUEST['graph'];
  $clusters = $_REQUEST['clusters'];
  if ($clusters == "" && $_REQUEST['pipe_clusters_file']) {
    $pipe_clusters_file = $_REQUEST['pipe_clusters_file'];
    $clusters = storeFile($pipe_clusters_file);
  }
  $s_col = $_REQUEST['s_col'];
  $t_col = $_REQUEST['t_col'];
  $w_col = $_REQUEST['w_col'];
  $decimals = $_REQUEST['decimals'];
  $stat = $_REQUEST['stat'];
  
  ## If a file and a graph are submitted -> error
  if ($graph != "" && $graph_file != "") {
    $error = 1;
    error("You must not submit both a graph and a graph file");
  }
  ## If clusters and a cluster file are submitted -> error
  if ($clusters != "" && $cluster_file != "") {
    $error = 1;
    error("You must not submit both clusters and a cluster file");
  }
  ## No specification of the source and target columns
  if ($in_format == "tab" && $s_col == "" && $t_col == "" && $w_col = "") {
    warning("Default value for source and target columns for tab-delimited input format are 1 and 2, respectively");
  }
  ## put the content of the file $graph_file in $graph
  if ($graph_file != "" && $graph == "") {
    $graph = storeFile($graph_file);
  }
  
  ## put the content of the file $clusters_file in $clusters
  if ($clusters_file != "" && $clusters == "") {
    $clusters = storeFile($clusters_file);
  }
  ## If no graph is submitted -> error
  if ($graph == "" && $graph_file == "") {
    $error = 1;
    error("You must submit an input graph");
  }
  ## If no clusters are submitted -> error
  if ($clusters == "") {
    $error = 1;
    error("You must submit input clusters");
  }  
  if (!$error) { 
    $graph = trim_text($graph);
    $clusters = trim_text($clusters);
    ## Load the parameters of the program in to an array
    $parameters = array( 
      "request" => array(
        "informat"=>$in_format,
        "clusters"=>$clusters,
        "inputgraph"=>$graph,
        "scol"=>$s_col,
        "tcol"=>$t_col,
        "wcol"=>$w_col,
	"stat"=>$stat,
	"decimals"=>$decimals,
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
    try {
      $echoed = $client->graph_cluster_membership($parameters);
      $soap_error = 0;
    } catch (Exception $soap_exception) {
      echo ("<pre>");
      echo "Error : \n",  $soap_exception->getMessage(), "\n";
      echo ("</pre>");
      $soap_error = 1;
    } 

if ($soap_error!=1) {
    $response =  $echoed->response;
    $command = $response->command;
    $server = $response->server;
    $client_result = $response->client;
    $server = rtrim ($server);
    $temp_file = explode('/',$server);
    $temp_file = end($temp_file);
    $resultURL = $WWW_RSA."/tmp/".$temp_file;
    #comment_file 
    $comment_file = $server.".comments";
    $comments = storeFile($comment_file);
    if ($comments != "") {
      warning($comments);
      echo "<hr>";
    }
    # Text-to-html
    $file = storeFile($server);
    # echo "</pre>"; 
    #$tth_parameters = array( 
    #  "request" => array(
     #  "inputfile"=>$file,
     #  "chunk"=>1000,
     #  "no_sort"=>1
    #  )
   # );
  #  $tth_echoed = $client->text_to_html($tth_parameters);
    
    #$tth_response =  $tth_echoed->response;
    #$tth_command = $tth_response->command;
   # $tth_server = $tth_response->server;
   # $tth_client = $tth_response->client;
    #echo "</pre>";
    #$tth_server = rtrim ($tth_server);
   # $tth_temp_file = explode('/',$tth_server);
   # $tth_temp_file = end($tth_temp_file);
   # $tth_resultURL = $WWW_RSA."/tmp/".$tth_temp_file;
  
    # Draw heatmap (low resolution)
     echo "</pre>"; 
    $ldh_parameters = array( 
      "request" => array(
       "row_names"=>1,
       "inputfile"=>$file,
       "outformat"=>"png",
       "col_width"=>"10",
       "row_height"=>"10",
//        "html"=>1,
       "header"=>1
      )
    );
    $ldh_echoed = $client->draw_heatmap($ldh_parameters);
    $ldh_response = $ldh_echoed->response;
    $ldh_command = $ldh_response->command;
    $ldh_server = $ldh_response->server;
    $ldh_client = $ldh_response->client;
    echo "</pre>";
    $ldh_server = rtrim ($ldh_server);
    $ldh_temp_file = explode('/',$ldh_server);
    $ldh_temp_file = end($ldh_temp_file);
    $ldh_resultURL = $WWW_RSA."/tmp/".$ldh_temp_file;
    $ldh_resultURL_map = $WWW_RSA."/tmp/".$ldh_temp_file.".html";
    
    # Draw heatmap (high resolution)
     echo "</pre>"; 
    $hdh_parameters = array( 
      "request" => array(
       "row_names"=>1,
       "inputfile"=>$file,
       "outformat"=>"png",
       "header"=>1
      )
    );
    $hdh_echoed = $client->draw_heatmap($hdh_parameters);
    $hdh_response = $hdh_echoed->response;
    $hdh_command = $hdh_response->command;
    $hdh_server = $hdh_response->server;
    $hdh_client = $hdh_response->client;
    echo "</pre>";
    $hdh_server = rtrim ($hdh_server);
    $hdh_temp_file = explode('/',$hdh_server);
    $hdh_temp_file = end($hdh_temp_file);
    $hdh_resultURL = $WWW_RSA."/tmp/".$hdh_temp_file;  
  
  
    hourglass("off");
    $lowres_size = urlfilesize($ldh_resultURL, "mb");
    $highres_size = urlfilesize($hdh_resultURL, "mb");
    
    
    # Display the results
    echo "The membership table (text file) is available at the following URL";
    echo "<a href = '$resultURL'>$resultURL</a><br>";
    #echo "The results are available as HTML page at the following URL ";
    #echo "<a href = '$tth_resultURL'>$tth_resultURL</a><br>"; 
    echo "<hr>";
    info ("Depending on the web browser, and of the number of clusters and nodes, images might not be viewable as the dimension of the figure are too large. In this case, right click (save as) and open it with an other visualization tool");
    echo "The heat map is available at the folllowing URL<br>";
    echo "&nbsp;&nbsp;low resolution ($lowres_size Mb) <a href = '$ldh_resultURL'>$ldh_resultURL</a>.<br>"; 
    echo "&nbsp;&nbsp;high resolution ($highres_size Mb) <a href = '$hdh_resultURL'>$hdh_resultURL</a>.<br>"; 
    echo "<hr>\n";
 
  	}
  	
  }
?>
