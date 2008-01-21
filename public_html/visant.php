<html>
<head>
   <title>NAT - Pathfinder</title>
   <link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
</head>
<body class="results">
<?php
require ('functions.php');
# log file update thanks to Sylvain
UpdateLogFile("neat","","");
# default variables
$out_format="VisML";
$directed = 0;
$nodeAttribs = "";
$edgeAttribs = "";
$noNodes = 0;
$return_type="server";
$error = 0;

############## prepare parameters ##########################

# read in params
$graph_location = $_REQUEST['visant_graph_file'];
$graph_format = $_REQUEST['visant_graph_format'];
$directed = $_REQUEST['visant_directed'];

# check parameters
if($graph_format == ""){
    error("You did not specify the format of the graph to be converted into visML");
}else if($graph_format == 'tab'){
    $graph_format = 'flat';
}else if($graph_format == 'gml'){
    $graph_format = 'GML';
}

# read in graph
if ($graph_location != "") {
    $graph = storeFile($graph_location);
}else{
    error("You did not specify the location of the graph to be converted into visML!");
}

############## call conversion web service #################

  $parameters = array(
      "request" => array(
        'graphString'=>$graph,
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
     $client = new SoapClient(
                      'http://rsat.scmbb.ulb.ac.be/be.ac.ulb.bigre.graphtools.server/wsdl/GraphAlgorithms.wsdl',
                           array(
                                 'trace' => 1,
                                 'soap_version' => SOAP_1_1,
                                 'style' => SOAP_DOCUMENT,
                                 'encoding' => SOAP_LITERAL
                                 )
                           );
    $functions = $client->__getFunctions();
    #info(print_r($parameters));
    #info(print_r($functions));
    try{
        $echoed = $client->graphconversion($parameters);
    }catch(SoapFault $fault){
        echo("The following error occurred:");
        error($fault);
        $error = 1;
    }

################ process results ##########################

    $response =  $echoed->response;
    # result processing
    $command = $response->command;
    $server = $response->server;
    $client = $response->client;

    if($error == 0){
        $address = " http://visant.bu.edu:8080/vserver/DAI?command=link&location=";
         # Display the results
    	echo("<align='left'>Click on the image below to display graph via visANT java web start:<br><br>
    	<a href='$address$server'><img src='images/visant_icon.png' alt='VisANT icon' name'VisANT icon' width='150' border='1'></a></align>
    	<br><br>Download graph in VisML format:<br><br>
    	<a href='$server'>$server</a>
    	<br><br><br><br>
    	More about VisANT: see <a href='http://visant.bu.edu/' target='_blank'>VisANT homepage</a>");
    }

?>
</body>
</html>