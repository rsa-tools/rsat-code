<html>
<head>
   <title>NAT Pathfinder</title>
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
  # log file update thanks to Sylvain
  UpdateLogFile("neat","","");
  # File to store the commands
  $cmd_file = getTempFileName('commands_pathfinder', '.txt'); 
  $cmd_handle = fopen($cmd_file, 'a');
  title('Results Pathfinder');
  # Error status
  $error = 0;
  # default params
  $default_returnType = "server";
  $default_attribs = "";
  $algorithm = "rea";
  $tab_java = 1;
  $noHTMLTable = 1;
  $tmpGraphFile = "";
  # server-related params
  $result_location = $tmp.'/';
  $html_location = $WWW_RSA.'/tmp/';
  $sylvain_input_format = "tab";
  $sylvain_input_graph = "";

  # Get parameters
  $piped_graph = $_REQUEST['pipe_graph_file'];
  $piped_in_format = $_REQUEST['pipe_in_format'];
  $in_format = $_REQUEST['in_format'];
  $out_format = $_REQUEST['out_format'];
  $sources = $_REQUEST['sources'];
  $targets = $_REQUEST['targets'];
  $to_exclude = $_REQUEST['exclude'];
  $to_include = $_REQUEST['include'];
  if ($_FILES['batch_file']['name'] != "") {
    $batch_file = uploadFile('batch_file');
  }
  if ($_FILES['graph_file']['name'] != "") {
    $graph_file = uploadFile('graph_file');
  }
  $graph = $_REQUEST['graph'];
  $graph_id = $_REQUEST['graph_id'];
  $directed = $_REQUEST['directed'];
  $weight = $_REQUEST['weight'];
  $rank = $_REQUEST['rank'];
  $outputchoice = $_REQUEST['outputchoice'];
  $outputType =   $_REQUEST['outputType'];
  $store_graph = $_REQUEST['store_graph'];
  $return_type = $default_returnType;
  # Get advanced parameters
  $algorithm = $_REQUEST['algorithm'];
  $maxWeight = $_REQUEST['maxWeight'];
  $maxLength = $_REQUEST['maxLength'];
  $minLength = $_REQUEST['minLength'];
  $exAttrib = $_REQUEST['exAttrib'];
  $metabolic = $_REQUEST['metabolic'];
  $email = $_REQUEST['email'];
  # $neat_java_wsdl = 'http://localhost:8080/be.ac.ulb.bigre.graphtools.server/services/GraphAlgorithms';

  ############ Check input ########################

  if($piped_graph != ""){
		$graph_file = $piped_graph;
		$in_format = $piped_in_format;
  }

   if ($store_graph == 'on') {
    $store_graph = 1;
  } else {
    $store_graph = 0;
  }

  if ($directed == 'on') {
    $directed = 1;
  } else {
    $directed = 0;
  }

  if($metabolic == 'on'){
    $metabolic = 1;
  }else{
    $metabolic = 0;
  }

  ## check outputType
  if($outputchoice == 'pathsTable'){
    $outputType = $outputchoice;
  }

  ## convert format names
  if($out_format == 'GML'){
  	$sylvain_input_format = 'gml';
  }
  if($out_format == 'flat'){
  	$sylvain_input_format = 'tab';
  }

  ## forbidden to set both batchfile and source and target nodes
  if($batch_file != "" && $sources != "" && $targets != ""){
  	$error = 1;
  	error("Set either a path finding batch file or source and target nodes, but not both.");
  }

  ## set batch file
  if($batch_file != ""){
  	$sources =  storeFile($batch_file);
  	$targets = "";
  }

  ## If a file and a graph are submitted -> error
  if ($graph != "" && $graph_file != "") {
    $error = 1;
    error("You must not submit both a graph and a graph file");
  }

  ## If a graph id and a graph file are submitted -> error
  if ($graph_id != "" && $graph_file != "") {
    $error = 1;
    error("You must not submit both a graph id and a graph file");
  }

   ## If a graph id and a graph  are submitted -> error
  if ($graph_id != "" && $graph != "") {
    $error = 1;
    error("You must not submit both a graph id and a graph.");
  }

  ## No specification of the source and target nodes or of a batch file
  if ($source == "" && $targets == "" && $batch_file == "") {
    error("You need to specify source and target nodes or a batch file.");
  }

  ## put the content of the file $graph_file in $graph
  if ($graph_file != "" && $graph == "") {
    $graph = storeFile($graph_file);
  }

  ## put the content of the graph id into graph
  if($graph_id != "" && $graph == ""){
  	$graph = $graph_id;
  }

  ## If no graph are submitted -> error
  if ($graph == "" && $graph_file == "") {
    $error = 1;
    error("You must submit an input graph");
  }

  if (!$error) {
    # convert two spaces in a row into a tab delimiter
    if(strcmp($out_format,'flat') == 0){
        $graph = spaces_to_tab($graph,2);
    }

   ########## Launch the client ###############

    $parameters = array(
      "request" => array(
     	'source'=>$sources,
     	'target'=>$targets,
     	'graphString'=>$graph,
     	'tmpInGraphFile'=>$tmpGraphFile,
        'inFormat'=>$in_format,
        'outFormat'=>$out_format,
    	'directed'=>$directed,
    	'metabolic'=>$metabolic,
    	'exclusionAttr'=>$exAttrib,
    	'weight'=>$weight,
    	'algorithm'=>$algorithm,
    	'rank'=>$rank,
    	'nodesPresent'=>$to_include,
    	'nodesAbsent'=>$to_exclude,
    	'attribs'=>$default_attribs,
    	'maxWeight'=>$maxWeight,
    	'maxLength'=>$maxLength,
    	'minLength'=>$minLength,
    	'outputType'=>$outputType,
    	'storeInputGraph'=>$store_graph,
    	'returnType'=>$return_type
      )
    );

    if($email == ""){
        ## announce the results
        echo("<br><br><br><br>");
        info("Results will appear below");
        echo"<hr>\n";
        echo("<div id='hourglass' class='hourglass'><img src='images/animated_hourglass.gif' height='50' border='1'></div>");
        flush();

         # Open the SOAP client
        # echo($neat_java_wsdl);
        $client = new SoapClient(
                      $neat_java_wsdl,
                           array(
                                 'trace' => 1,
                                 'soap_version' => SOAP_1_1,
                                 'style' => SOAP_DOCUMENT,
                                 'encoding' => SOAP_LITERAL
                                 )
                           );
        # Execute the command
        $functions = $client->__getFunctions();
        $types = $client->__getTypes();
        #info(print_r($parameters));
        #info(print_r($functions));
        #info(print_r($types));
 	  try{
        $echoed = $client->pathfinding($parameters);
        }catch(SoapFault $fault){
        echo("The following error occurred:");
        error($fault);
        $error = 1;
        }

    echo("<div id='hide' class='hide'><img src='images/hide_hourglass.jpg' height='60' border='0'></div>");
    ########## Process results ###############

    $response =  $echoed->response;
    # result processing
    $command = $response->command;
    $server = $response->server;
    $client = $response->client;
    $graphid = $response->graphid;
    if(ereg('PATHFINDER ERROR',$server)){
    	$error = 1;
    	error("$server");
    }
    if($error == 0){
   		# location of result file on server (absolute path)
    	$file_location = $result_location . $server;
        # content of result file
        $fileContent = storeFile($file_location);
        # Display warning if required
        if(ereg('maximal path number',$client)){
        	echo("<font color='red'>Warning: Path enumeration had to be interrupted, because too many paths were found. Paths might have been missed!</font>");
        }else if(ereg('time out',$client)){
        	echo("<font color='red'>Warning: Path enumeration was interrupted by a time out. Paths might have been missed!</font>");
        }
        # Display the results
    	# echo($command);
    	# $resultURL = rsat_path_to_url($server); results in error when clicking result link
    	$resultURL = $html_location.$server;
        $URL['result'] = $resultURL;
			
    	# store command in a file
        store_command("$command", "k-shortest path finding", $cmd_handle);

    	# Text-to-html web service (for table of paths only)
    	if(strcmp($outputType,'pathsTable') == 0 && $noHTMLTable == 0){
    	 $rsat_client = new SoapClient(
                          $neat_wsdl,
                           array(
                                 'trace' => 1,
                                 'soap_version' => SOAP_1_1,
                                 'style' => SOAP_DOCUMENT,
                                 'encoding' => SOAP_LITERAL
                                 )
                           );
        $tth_parameters = array(
          "request" => array(
          "inputfile"=>$fileContent,
          "chunk"=>1000,
         	)
        );

        $tth_echoed = $rsat_client->text_to_html($tth_parameters);
        $tth_response =  $tth_echoed->response;
        $tth_command = $tth_response->command;
        $tth_server = $tth_response->server;
       	$tth_server = rtrim ($tth_server);
        $tth_temp_file = explode('/',$tth_server);
   	    $tth_temp_file = end($tth_temp_file);
    	$tth_resultURL = $WWW_RSA."/tmp/".$tth_temp_file;
    	
    	store_command("$tth_command", "Text to html", $cmd_handle);
    	
    	$tth_resultURL = rsat_path_to_url($tth_server);
        $URL['HTML'] = $tth_resultURL;
    	
    	} # end pathsTable
    	
    	        
         ## Close command handle
         fclose($cmd_handle);
         $URL['Server commands'] = rsat_path_to_url($cmd_file);
         ## DISPLAY THE RESULT
         print_url_table($URL);
    	
    	# in case of tab-format, truncate nodes to make it readable by Sylvain Brohee's tools
    	if(strcmp($out_format,'flat') == 0){
    		if(ereg(';ARCS',$fileContent)){
    		    # get string without nodes
    			$fileContent = end(explode(';ARCS	rgb_color	color',$fileContent));
    			$sylvain_input_graph = $fileContent;
    		}
    	}else{
    		$sylvain_input_graph = $fileContent;
   	    }
    	# remove leading or trailing white spaces or end of lines
    	$sylvain_input_graph = ltrim($sylvain_input_graph,"\n");
    	$sylvain_input_graph = rtrim($sylvain_input_graph,"\n");
    	$sylvain_input_graph = ltrim($sylvain_input_graph);
    	$sylvain_input_graph = rtrim($sylvain_input_graph);

	   # generate temp file
	   $tempFileName = tempnam($result_location,"Pathfinder_tmpGraph_");
	   $fh = fopen($tempFileName, 'w') or die("Can't open file $tempFileName");
	   fwrite($fh, $sylvain_input_graph);
	   fclose($fh);

        if($store_graph) {
   		   echo "<br><align='left'>Your stored input graph has the id:<br> $graphid<br>
   		   Submit this id to speed up other path finding jobs on this input graph.</align>";
        }
    echo "<hr>\n";
    if(strcmp($outputType,'pathsTable') == 0){
         echo "
     	To process your result with another tool, click one of the buttons listed below.
     	<br>
     	<br>
 	  <TABLE CLASS = 'nextstep'>
  		<TR>
      	<Th colspan = 3>Next steps</Th>
    	</TR>
    	<TR>
      	<TD>
       <FORM METHOD='POST' ACTION='convert_graph_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$file_location'>
          <input type='hidden' NAME='graph_format' VALUE='path'>
          <input type='hidden' NAME='pathcol' VALUE='7'>
          <INPUT type='submit' value='Convert the paths table into a graph with paths unified'>
         </form>
        </td>
        <td>
         <FORM METHOD='POST' ACTION='convert_graph_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$file_location'>
          <input type='hidden' NAME='graph_format' VALUE='path'>
          <input type='hidden' NAME='distinct_path' VALUE='on'>
          <input type='hidden' NAME='pathcol' VALUE='7'>
          <INPUT type='submit' value='Convert the paths table into a graph with paths separated'>
         </form>
        </td>
        </tr>
        </table>";

    }
    if(strcmp($outputType,'pathsUnion') == 0 || strcmp($outputType,'pathsMarked') == 0 || strcmp($outputType,'pathsGraphs') == 0){
     echo "
     	To process your result with another tool, click one of the buttons listed below.
     	<br>
     	<br>
 	  <TABLE CLASS = 'nextstep'>
  		<TR>
      	<Th colspan = 3>Next steps</Th>
    	</TR>
    	<TR>
      	<TD>
        <FORM METHOD='POST' ACTION='display_graph_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$tempFileName'>
          <input type='hidden' NAME='graph_format' VALUE='$sylvain_input_format'>";
          if ($sylvain_input_format == 'tab') {
            echo "
            <input type='hidden' NAME='scol' VALUE='1'>
            <input type='hidden' NAME='tcol' VALUE='2'>
            <input type='hidden' NAME='eccol' VALUE='3'>";
          }
          echo "
          <INPUT type='submit' value='Display the graph'>
         </form>
        </td>
        <TD>
        <FORM METHOD='POST' ACTION='visant.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='visant_graph_file' VALUE='$file_location'>
          <input type='hidden' NAME='visant_graph_format' VALUE='$sylvain_input_format'>
          <input type='hidden' NAME='visant_directed' VALUE='$directed'>
          <input type='hidden' NAME='tab_java' VALUE='$tab_java'>
          <INPUT type='submit' value='Display the graph with VisANT'>
         </form>
        </td>
       <TD>
         <FORM METHOD='POST' ACTION='compare_graphs_form.php'>
           <input type='hidden' NAME='pipe' VALUE='1'>
           <input type='hidden' NAME='graph_file' VALUE='$tempFileName'>
           <input type='hidden' NAME='graph_format' VALUE='$sylvain_input_format'>";
          if ($sylvain_input_format == 'tab') {
            echo "
             <input type='hidden' NAME='scol' VALUE='1'>
             <input type='hidden' NAME='tcol' VALUE='2'>
             <input type='hidden' NAME='wcol' VALUE='3'>";
          }
          echo "
          <INPUT type='submit' value='Compare this graph to another one'>
         </form>
       </td>
     </tr>
     <tr>
     	<TD>
        <FORM METHOD='POST' ACTION='pathfinder_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$file_location'>
          <input type='hidden' NAME='in_format' VALUE='$out_format'>";
          echo "
          <INPUT type='submit' value='Do path finding on this graph'>
         </form>
        </td>
       <TD>
         <FORM METHOD='POST' ACTION='graph_get_clusters_form.php'>
           <input type='hidden' NAME='pipe' VALUE='1'>
           <input type='hidden' NAME='graph_file' VALUE='$tempFileName'>
           <input type='hidden' NAME='graph_format' VALUE='$sylvain_input_format'>";
          if ($sylvain_input_format == 'tab') {
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
          <input type='hidden' NAME='graph_file' VALUE='$tempFileName'>
          <input type='hidden' NAME='graph_format' VALUE='$sylvain_input_format'>";
          if ($sylvain_input_format == 'tab') {
            echo "
             <input type='hidden' NAME='scol' VALUE='1'>
             <input type='hidden' NAME='tcol' VALUE='2'>
             <input type='hidden' NAME='wcol' VALUE='3'>";
           }
           echo "
          <INPUT type='submit' value='Node degree computation'>
         </form>
       </td>
     </tr>
     <tr>
      <TD>
        <FORM METHOD='POST' ACTION='convert_graph_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$tempFileName'>
          <input type='hidden' NAME='graph_format' VALUE='$sylvain_input_format'>";
          if ($out_format == 'tab') {
            echo "
            <input type='hidden' NAME='scol' VALUE='1'>
            <input type='hidden' NAME='tcol' VALUE='2'>";
          }
          echo "
          <INPUT type='submit' value='Convert $sylvain_input_format to another format'>
        </form>
      </td>
       <TD>
        <FORM METHOD='POST' ACTION='graph_neighbours_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$tempFileName'>
          <input type='hidden' NAME='graph_format' VALUE='$sylvain_input_format'>";
          if ($out_format == 'tab') {
            echo "
            <input type='hidden' NAME='scol' VALUE='1'>
            <input type='hidden' NAME='tcol' VALUE='2'>";
          }
          echo "
          <INPUT type='submit' value='Neighbourhood analysis'>
        </form>
      </td>
       <TD>
        <FORM METHOD='POST' ACTION='mcl_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$tempFileName'>
          <input type='hidden' NAME='graph_format' VALUE='$sylvain_input_format'>";
          if ($out_format == 'tab') {
            echo "
            <input type='hidden' NAME='scol' VALUE='1'>
            <input type='hidden' NAME='tcol' VALUE='2'>";
          }
          echo "
          <INPUT type='submit' value='MCL Graph clustering'>
        </form>
      </td>
     </tr>
     <tr>
      <TD>
        <FORM METHOD='POST' ACTION='random_graph_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$tempFileName'>
          <input type='hidden' NAME='graph_format' VALUE='$sylvain_input_format'>";
          if ($sylvain_input_format == 'tab') {
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
   </table>";
  		} # if output type is a graph
  	  } # path finding error
  	  # send result via email
      }else{

        $pathfinderParams = array(
     	'source'=>$sources,
     	'target'=>$targets,
     	'graphString'=>$graph,
     	'tmpInGraphFile'=>$tmpGraphFile,
        'inFormat'=>$in_format,
        'outFormat'=>$out_format,
    	'directed'=>$directed,
    	'metabolic'=>$metabolic,
    	'exclusionAttr'=>$exAttrib,
    	'attribs'=>$default_attribs,
    	'weight'=>$weight,
    	'algorithm'=>$algorithm,
    	'rank'=>$rank,
    	'nodesPresent'=>$to_include,
    	'nodesAbsent'=>$to_exclude,
    	'maxWeight'=>$maxWeight,
    	'maxLength'=>$maxLength,
    	'minLength'=>$minLength,
    	'outputType'=>$outputType,
    	'storeInputGraph'=>$store_graph,
    	'returnType'=>$return_type
         );

        $mixedRequest = array("PathfinderRequest"=>$pathfinderParams,"GraphConverterRequest"=>NULL, "MetabolicGraphConstructorRequest"=>NULL,"PathwayinferenceRequest"=>NULL);
        $requestArray = array(0=>$mixedRequest);
        $emailParams = array("request" => array(
                'email'=>$email,
                'requestArray'=>$requestArray
            )
        );
        # Open the SOAP client
        $emailclient = new SoapClient(
                           $neat_java_wsdl,
                           array(
                                 'trace' => 1,
                                 'soap_version' => SOAP_1_1,
                                 'style' => SOAP_DOCUMENT,
                                 'encoding' => SOAP_LITERAL
                                 )
                           );
 	   try{
 	      $functions = $emailclient->__getFunctions();
          $types = $emailclient->__getTypes();
 	      # info(print_r($emailParams));
 	      # info(print_r($types));
          $echoed = $emailclient->workflow($emailParams);
        }catch(SoapFault $fault){
            echo("The following error occurred:");
            error($fault);
            $error = 1;
        }
        if($error == 0){
            info("After computation has finished, the result will be sent to the specified email address: ".$email);
        }
    } # end send by email
} # input param error
?>
</body>
</html>
