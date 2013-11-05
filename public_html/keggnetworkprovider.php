<html>
<head>
   <title>NAT - KEGG network provider</title>
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
# File to store the commands
$cmd_file = getTempFileName('commands_keggnetworkprovider', '.txt');
$cmd_handle = fopen($cmd_file, 'a');
# Sylvain's upload function modified to accept upload location
Function uploadFileToGivenLocation($file, $location) {
    $nomDestination = $_FILES[$file]["name"];
    $now = date("Ymd_His");
    $nomDestination = $nomDestination.$now;

    if (is_uploaded_file($_FILES[$file]['tmp_name'])) {
        if (rename($_FILES[$file]['tmp_name'], $location.$nomDestination)) {
        } else {
            echo "Could not move $_FILES[$file]['tmp_name']"." check that $location exists<br>";
        }
    } else {
       echo "File ".$_FILES[$file]['tmp_name']." could not be uploaded<br>";
    }
    return $location.$nomDestination;
}

title('Results KEGG network provider');
# default variables
$organisms = "";
$organisms_file = "";
$reactions = "";
$reactions_file = "";
$excludedcompounds = "";
$excludedreactions = "";
$excludedrpairclasses = "";
$directed = 0;
$rpairs = 0;
$keepIrre = 0;
$attribs = "";
$out_format = "";
$return_type="server";
$email = "";
$error = 0;
$tab_java = 1;
# helper variables
$selected_attribs = "";
$fileContent = "";
$sylvain_input_graph = "";
$graph_type = "";
$customgraphdb = "";
$symmetric = "true";

# server-related params
$result_location = $tmp.'/';
$html_location = $WWW_RSA.'/tmp/';
# $host= parse_url($WWW_RSA,PHP_URL_HOST);
$metabolicpathfinder_location = $neat_java_host.'/metabolicpathfinding/metabolicPathfinder_form.jsp';
$pathwayinference_location = $neat_java_host.'/metabolicpathfinding/pathwayinference_form.jsp';

############## prepare parameters ##########################

# read in params
$organisms = $_REQUEST['organisms'];
if ($_FILES['organisms_file']['name'] != "") {
    $organisms_file = uploadFileToGivenLocation('organisms_file',$result_location);
    if($organisms_file != ""){
        # file is readable for everyone, but writable for owner only
    	chmod($organisms_file,644);
    }
}
$reactions = $_REQUEST['reactions'];
if ($_FILES['reactions_file']['name'] != "") {
    $reactions_file = uploadFileToGivenLocation('reactions_file',$result_location);
     if($reactions_file != ""){
    	chmod($reactions_file,644);
    }
}
$excludedcompounds = $_REQUEST['excludedcompounds'];
$excludedreactions = $_REQUEST['excludedreactions'];
$excludedrpairclasses = $_REQUEST['excludedrpairclasses'];
$directed = $_REQUEST['directed'];
$keepIrre = $_REQUEST['keepIrre'];
$rpairs = $_REQUEST['rpair'];
# name of parameter in form is compoundattributes[] and reactionattributes[]
$compoundattribs = $_REQUEST['compoundattributes'];
$reactionattribs = $_REQUEST['reactionattributes'];
$out_format = $_REQUEST['outFormat'];
$email = $_REQUEST['email'];

# check parameters

if ($directed == 'on') {
    $directed = 1;
} else {
    $directed = 0;
}

if ($keepIrre == 'true') {
    $keepIrre = 1;
} else {
    $keepIrre = 0;
}

  if($rpairs == 'on'){
    $rpairs = 's';
  }else{
    $rpairs = 'r';
  }

# separate file names from path
$organisms_file = basename($organisms_file);
$reactions_file = basename($reactions_file);

# load attribs from multiple select form
if($compoundattribs){
 foreach ($compoundattribs as $t){
 	$selected_attribs .= $t."/";
 }
}
if($reactionattribs){
 foreach ($reactionattribs as $t){
 	$selected_attribs .= $t."/";
 }
}

############## call metabolicgraphconstruction web service #################

  $parameters = array(
      "request" => array(
        'organismNames'=>$organisms,
        'organismFile'=>$organisms_file,
        'reactionIds'=>$reactions,
        'reactionFile'=>$reactions_file,
    	'directed'=>$directed,
    	'keepIrreversible'=>$keepIrre,
    	'graphType'=>$rpairs,
    	'attributes'=>$selected_attribs,
    	'outFormat'=>$out_format,
    	'excludeCompounds'=>$excludedcompounds,
    	'excludeReactions'=>$excludedreactions,
    	'excludeRPairClasses'=>$excludedrpairclasses,
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
    	# info(print_r($parameters));
    	# info(print_r($functions));
    	try{
      	  $echoed = $client->metabolicgraphconstruction($parameters);
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

    	if(ereg('KEGG organisms or custom reactions or both have to be specified',$server)){
    		$error = 1;
    		error("$server");
    	}

    	if($error == 0){

			# compact result display
			# $resultURL = rsat_path_to_url($server); results in error when clicking result link
			$resultURL = $html_location.$server;
			$URL['tab'] = $resultURL;
			
			# no html format for kegg network provider
			
			# store command in a file
            store_command("$command", "KEGG network generation", $cmd_handle);
            # Close command handle
            fclose($cmd_handle);
            $URL['Server commands'] = rsat_path_to_url($cmd_file);
			
			 ## DISPLAY THE RESULT
             print_url_table($URL);
			
    	# next step panel only available for graphs in gml or tab-delimited format
    	if(strcmp($out_format,'tab') == 0 || strcmp($out_format,'gml') == 0){
        	# location of result file
        	$file_location = $result_location . $server;

            # content of result file
        	$fileContent = storeFile($file_location);

        	# in case of tab-format, truncate nodes to make it readable by Sylvain Brohee's tools
    		if(strcmp($out_format,'tab') == 0){
    			if(ereg(';ARCS',$fileContent)){
    		    	# get string without nodes
    		    	if(ereg(';ARCS	RPAIRS.Linkage.Type',$fileContent)){
    					$fileContent = end(explode(';ARCS	RPAIRS.Linkage.Type',$fileContent));
    				}else{
    					$fileContent = end(explode(';ARCS',$fileContent));
    				}
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
	   		$tempFileName = tempnam($result_location,"Graphconstruction_tmpGraph_");
	   		$fh = fopen($tempFileName, 'w') or die("Can't open file $tempFileName");
	   		fwrite($fh, $sylvain_input_graph);
	   		fclose($fh);

			# Display the next steps panel
    		 echo("
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
          <input type='hidden' NAME='graph_format' VALUE='$out_format'>");
          if ($out_format == 'tab') {
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
          <input type='hidden' NAME='visant_graph_format' VALUE='$out_format'>
          <input type='hidden' NAME='visant_directed' VALUE='$directed'>
          <input type='hidden' NAME='tab_java' VALUE='$tab_java'>
          <INPUT type='submit' value='Display the graph with VisANT'>
         </form>
        </td>
       <TD>
         <FORM METHOD='POST' ACTION='compare_graphs_form.php'>
           <input type='hidden' NAME='pipe' VALUE='1'>
           <input type='hidden' NAME='graph_file' VALUE='$tempFileName'>
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
     </tr>
     <tr>
       <td>";
       if($rpairs == 's'){
       	$graph_type="subreactiongraph";
       	$customgraphdb = "keggrpair";
       }else{
       	$graph_type="reactiongraph";
       	$customgraphdb = "keggreaction";
       }
       if($directed == 1){
          $directed = "true";
       }else{
          $directed = "false";
       }
       if($keepIrre == 1){
       	  $symmetric = "false";
       }
       echo "
        <FORM METHOD='POST' ACTION='$metabolicpathfinder_location'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$file_location'>
          <input type='hidden' NAME='graph_format' VALUE='$out_format'>
          <input type='hidden' NAME='graph_type' VALUE='$graph_type'>
          <input type='hidden' NAME='directed' VALUE='$directed'>
          <input type='hidden' NAME='irreversible' VALUE='$keepIrre'>
          <input type='hidden' NAME='organisms' VALUE='$organisms'>
          <INPUT type='submit' value='Find metabolic paths in this graph'>
         </form>
       </td>
       <td>
       <FORM METHOD='POST' ACTION='$pathwayinference_location'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$file_location'>
          <input type='hidden' NAME='graph_format' VALUE='$out_format'>
          <input type='hidden' NAME='directed' VALUE='$directed'>
          <input type='hidden' NAME='metabolic' VALUE='true'>
          <input type='hidden' NAME='symmetric' VALUE='$symmetric'>
          <input type='hidden' NAME='irreversible' VALUE='$keepIrre'>
           <input type='hidden' NAME='customgraphdb' VALUE='$customgraphdb'>
          <input type='hidden' NAME='organisms' VALUE='$organisms'>
          <INPUT type='submit' value='Infer a pathway from this graph'>
         </form>
       </td>
       <TD>
         <FORM METHOD='POST' ACTION='graph_get_clusters_form.php'>
           <input type='hidden' NAME='pipe' VALUE='1'>
           <input type='hidden' NAME='graph_file' VALUE='$tempFileName'>
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
     </tr>
     <tr>
      <TD>
        <FORM METHOD='POST' ACTION='convert_graph_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$tempFileName'>
          <input type='hidden' NAME='graph_format' VALUE='$out_format'>";
          if ($out_format == 'tab') {
            echo "
            <input type='hidden' NAME='scol' VALUE='1'>
            <input type='hidden' NAME='tcol' VALUE='2'>";
          }
          echo "
          <INPUT type='submit' value='Convert $out_format to another format'>
        </form>
      </td>
       <TD>
        <FORM METHOD='POST' ACTION='graph_neighbours_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$tempFileName'>
          <input type='hidden' NAME='graph_format' VALUE='$out_format'>";
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
          <input type='hidden' NAME='graph_format' VALUE='$out_format'>";
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
       <FORM METHOD='POST' ACTION='graph_node_degree_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$tempFileName'>
          <input type='hidden' NAME='graph_format' VALUE='$out_format'>";
          if ($out_format == 'tab') {
            echo "
             <input type='hidden' NAME='scol' VALUE='1'>
             <input type='hidden' NAME='tcol' VALUE='2'>
             <input type='hidden' NAME='wcol' VALUE='3'>";
           }
           echo "
          <INPUT type='submit' value='Node degree computation'>
         </form>
       </td>
      <TD>
        <FORM METHOD='POST' ACTION='random_graph_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$tempFileName'>
          <input type='hidden' NAME='graph_format' VALUE='$out_format'>";
          if ($out_format == 'tab') {
            echo "
             <input type='hidden' NAME='scol' VALUE='1'>
             <input type='hidden' NAME='tcol' VALUE='2'>
             <input type='hidden' NAME='wcol' VALUE='3'>";
          }
          echo("
          <INPUT type='submit' value='Randomize this graph'>
         </form>
       </td>
       <TD>
        <FORM METHOD='POST' ACTION='pathfinder_form.php'>
          <input type='hidden' NAME='pipe' VALUE='1'>
          <input type='hidden' NAME='graph_file' VALUE='$file_location'>
          <input type='hidden' NAME='in_format' VALUE='$out_format'>
          <INPUT type='submit' value='Do path finding on this graph'>
         </form>
        </td>
     </tr>
   		</table>");
		}else{
		   # not the right format for next step panel
		}
      }else{
    		echo("An error occurred. Could not construct the requested KEGG network.");
      }
    # send results by email
    }else{
        $parameters = array(
        'organismNames'=>$organisms,
        'organismFile'=>$organisms_file,
        'reactionIds'=>$reactions,
        'reactionFile'=>$reactions_file,
    	'directed'=>$directed,
    	'keepIrreversible'=>$keepIrre,
    	'graphType'=>$rpairs,
    	'attributes'=>$selected_attribs,
    	'outFormat'=>$out_format,
    	'excludeCompounds'=>$excludedcompounds,
    	'excludeReactions'=>$excludedreactions,
    	'excludeRPairClasses'=>$excludedrpairclasses,
    	'returnType'=>$return_type
    	);
	   $mixedRequest = array("MetabolicGraphConstructorRequest"=>$parameters,"GraphConverterRequest"=>NULL, "PathfinderRequest"=>NULL,"PathwayinferenceRequest"=>NULL);
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
            info("After network construction has finished, the result will be sent to the specified email address: ".$email);
        }
    }

?>
</body>
</html>
