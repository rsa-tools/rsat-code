<html>
<head>
   <title>NAT - VisANT</title>
   <link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
   <style type="text/css">
    <!--
    div.hourglass{position: absolute; top: 80px; left: 500px }
    div.hide{position: absolute; top: 80px; left: 500px }
   -->
</style>
</head>
<body class="results">
<?php
require ('functions.php');
# log file update thanks to Sylvain
UpdateLogFile("neat","","");

title('Display graph in VisANT');
# default variables
$out_format="VisML";
$directed = 0;
$nodeAttribs = "";
$edgeAttribs = "";
$noNodes = 0;
$return_type="server";
$error = 0;
$tmpGraphFile = "";
# server-related params
$result_location = $tmp.'/';
$html_location = $WWW_RSA.'/tmp/';
# visant address
$address = "http://visant.bu.edu:8080/vserver/DAI?command=link&location=";

############## prepare parameters ##########################

# read in params
$graph_location = $_REQUEST['visant_graph_file'];
$graph_format = $_REQUEST['visant_graph_format'];
$directed = $_REQUEST['visant_directed'];
$tab_java = $_REQUEST['tab_java'];

# check parameters
if($graph_format == ""){
    error("You did not specify the format of the graph to be converted into visML");
}else if($graph_format == 'tab'){
    $graph_format = 'flat';
}else if($graph_format == 'gml'){
    $graph_format = 'GML';
}
# check whether the graph comes from NeAT java tools ($tab_java = 1) or NeAT perl tools ($tab_java = 0)
if ($tab_java == "") {
  $tab_java = 1;
}
# read in graph
if ($graph_location != "") {
    $graph = storeFile($graph_location);
}else{
    error("You did not specify the location of the graph to be converted into visML!");
}

# if tab_java = 0 and $graph_format = flat ... conversion to the tab_java format (use of convert-graph)
if($tab_java == 0 && $graph_format == 'flat'){
$parameters = array(
      "request" => array(
        "informat"=>"tab",
        "outformat"=>"tab_java",
        "inputgraph"=>$graph,
        "scol"=>1,
        "tcol"=>2,
        "wcol"=>3,
        "eccol"=>4
      )
    );


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
    # Work with exception catch
    try {
      $echoed = $client->convert_graph($parameters);
      $soap_error = 0;
    } catch (Exception $soap_exception) {
      echo ("<pre>");
      echo "Error : \n",  $soap_exception->getMessage(), "\n";
      echo ("</pre>");
      $soap_error = 1;
    }
    if (!$soap_error) {
      $response =  $echoed->response;
      $command = $response->command;
      $server = $response->server;
      $client = $response->client;
      $server = rtrim ($server);
      $temp_file = explode('/',$server);
      $temp_file = end($temp_file);
      $resultURL = $WWW_RSA."/tmp/".$temp_file;
}

$graph = storeFile($server);
$graph = rtrim($graph);
}

############## call conversion web service #################

  $parameters = array(
      "request" => array(
        'graphString'=>$graph,
        'tmpInGraphFile'=>$tmpGraphFile,
        'inFormat'=>$graph_format,
        'outFormat'=>$out_format,
    	'directed'=>$directed,
        'nodeAttribs'=>$nodeAttribs,
        'edgeAttribs'=>$edgeAttribs,
        'noNodes'=>$noNodes,
    	'returnType'=>$return_type
      )
    );
    info_link("Graph uploaded from the previous treatment is converted into visML", rsat_path_to_url($graph_location));
    info("Results will appear below");
    echo"<hr>\n";
    echo("<div id='hourglass' class='hourglass'><img src='images/animated_hourglass.gif' height='50' border='1'></div>");
    flush();
    $client = new SoapClient(
                      $neat_java_wsdl,
                           array(
                                 'trace' => 1,
                                 'soap_version' => SOAP_1_1,
                                 'style' => SOAP_DOCUMENT,
                                 'encoding' => SOAP_LITERAL
                                 )
                           );
    $functions = $client->__getFunctions();
    try{
        $echoed = $client->graphconversion($parameters);
    }catch(SoapFault $fault){
        echo("The following error occurred:");
        error($fault);
        $error = 1;
    }
    echo("<div id='hide' class='hide'><img src='images/hide_hourglass.jpg' height='60' border='0'></div>");
################ process results ##########################

    $response =  $echoed->response;
    # result processing
    $command = $response->command;
    $server = $response->server;
    $client = $response->client;

    if($error == 0){
         # Display the results
    	echo("<align='left'>Click on the image below to display graph via visANT java web start:<br><br>
    	<a href='$address$html_location$server'><img src='images/visant_icon.png' alt='VisANT icon' name'VisANT icon' width='150' border='1'></a></align>
    	<br><br>Download graph in VisML format:<br><br>
    	<a href='$html_location$server'>$html_location$server</a>
    	<br><br><br><br>
    	More about VisANT: see <a href='http://visant.bu.edu/' target='_blank'>VisANT homepage</a>");
    }else{
    	echo("An error occurred. Could not display graph in VisANT.");
    }

?>
</body>
</html>