<html>
<head>
   <title>NAT - Pathfinder</title>
   <link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
</head>
<body class="results">
<?

# default variables
$out_format="VisML";
$nodeAttribs = "";
$edgeAttribs = "";
$noNodes = 0;
$return_type="server";

# read in params
$graph_location = $_REQUEST['visant_graph_file'];
$graph_format = $_REQUEST['visant_graph_format'];
$directed = $_REQUEST['visant_directed'];

# check parameters
if($graph_format == ""){
    error("You did not specify the format of the graph to be converted into visML");
}

# read in graph
if ($graph_location != "") {
    $graph = storeFile($graph_location);
}else{
    error("You did not specify the location of the graph to be converted into visML!");
}




# call conversion web service
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


?>
</body>
</html>