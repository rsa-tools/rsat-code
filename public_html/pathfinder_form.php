<html>
<head>
   <title>Network Analysis Tools - Pathfinder</title>
   <link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
</head>
<body class="form">
<?php
  require ('functions.php');

  // default variables
  $default_in_format = "";
  $default_out_format = "";
  $default_sources = "";
  $default_targets = "";
  $default_exclude = "";
  $default_include = "";
  $default_graph = "";
  $default_graph_id = "";
  $default_rank = 5;
  $default_weight = "unit";
  // default advanced params
  $default_maxWeight = 1000000;
  $default_maxLength = 1000000;
  $default_minLength = 0;
  $default_exAttrib = "";
  $default_algorithm = "rea";
  $default_metabolic = 0;
  $default_email='';
  $demo2_graph_id = 'Pathfinder_tmpGraph_d597cd86-3095-475f-99e7-d70a542d072a.tab';

  # variables given within workflow
  $pipe = $_REQUEST['pipe'];
  $requested_graph_file = $_REQUEST['graph_file'];
  $requested_in_format = $_REQUEST['in_format'];

  # advanced options
  $advanced = $_REQUEST['advanced'];

  # demo graph
  $demo = $_REQUEST['demo'];
  if ($demo == 1) {
     $demo_graph = storeFile("demo_files/Scer_biocyc.tab");
     $demo_sources = "GLY";
     $demo_targets = "PROTOHEME";
     $default_sources = $demo_sources;
     $default_targets = $demo_targets;
  }

  # demo 2 graph
  $demo2 = $_REQUEST['demo2'];
  if($demo2 == 1){
    $default_weight = "none";
    $default_sources = "MID2";
    $default_targets = "RLM1/SWI4/SWI6";
    $default_rank = 5;
  }

  # $host = parse_url($WWW_RSA,PHP_URL_HOST);
  $metabolic_pathfinder_location = $neat_java_host.'/metabolicpathfinding/metabolicPathfinder_form.jsp';

  title('Pathfinder');
  echo ("<center>Do multiple-to-multiple end path finding. Click on <a href='help.pathfinder.html'><img src='images/question-mark-button.jpg' width='15' alt='help'></a> for help.<br>
  Web service and interface by <a href='mailto:kfaust@ulb.ac.be'>Karoline Faust</a>.
  This program calls REA (Jimenez and Marzal, 1999). See
  <a href='help.pathfinder.html#credits'>Credits</a> for other contributors.</center>\n");
  echo ("<form method='post' action='pathfinder.php' enctype='multipart/form-data'>
  <hr>
  <h2>1. Network <a href='help.pathfinder.html#graph'><img src='images/question-mark-button.jpg' width='15' alt='help'></a>
  </h2>
  <br>
  <br>");
  if(!$pipe){
  	if(!$demo2){
  	echo("
  	Enter the network here (tab-delimited* or gml):
  	<br>
  	<textarea name='graph' rows='6' cols='65'>$demo_graph</textarea>
  	<br>
  	*The tab-delimiter can be replaced by two white spaces.
  	<br>
  	<br>
  	<b>Or</b> upload the network from a file:
  	<br>
  	<input type='file' name='graph_file' size='45' />
  	<br>
  	<br><b>OR</b> enter a network id (id of previously submitted network):
  	<br>
  	<input type='text' name='graph_id' size='45' />
  	<br>
  	<br>");
  	}else{
  	 echo("The network has been already uploaded. Its network identifier is: $demo2_graph_id<br>
  	 <input type='hidden' name='graph_id' value=$demo2_graph_id></input>
  	 ");
  	}
  } else {
  	  info_link("Network uploaded from the previous treatment", rsat_path_to_url($requested_graph_file));
  	echo "<input type='hidden' NAME='pipe_graph_file' VALUE='$requested_graph_file'>";
  }
  if($demo == 1){
  		demo("The demo network is the union of all paths annotated for <i>S. cerevisiae</i> in <a href='http://www.biocyc.org/' target='_blank'>BioCyc</a> release 10.6. It is an undirected graph consisting
  		of 2,656 edges.<br>
  		The seed nodes are the start and end compound of the <a href='http://biocyc.org/YEAST/NEW-IMAGE?object=HEME-BIOSYNTHESIS-II' target='_blank'>heme biosynthesis II pathway</a> as annotated in BioCyc.<br>
  		This pathway is one of the study cases described in Croes et al., J. Mol. Biol. 356: 222-236 (see our <a href='neat_publications.html '>list of publications</a>.)<br>
  		Note that for this demo, path finding is done on a smaller, undirected graph without differentially weighting compounds and reactions (both receive a weight according to their degree). To find paths in directed metabolic networks with mutually exclusive reaction directions, check the advanced options of the Pathfinder tool.
  		Use the <a href='$metabolic_pathfinder_location'> Metabolic path finding tool</a> to find paths in the complete KEGG network.<br>
  		The path of first rank corresponds to the annotated heme biosynthesis II pathway. To see the influence of the weighting scheme, you can set the weighting scheme to unit weight and the rank to 1 (for quicker computation). You will obtain an entirely different result.<br><br>");
  }else if($demo2 == 1){
  demo("In this demo, the network has been already pre-loaded. We use the weighted STRING database network, a network of protein-protein interactions, as demo network. Its weights, which describe edge reliability, have been inverted to describe
  edge costs. The demo network is part of our <a href='$WWW_RSA/data/published_data/nature_protocols/network_analysis/'>sample data collection</a>. We want to recover a signal
  transduction pathway known to regulate cell wall integrity in <i>S. cerevisiae</i> (Saito and Tatebayashi, J. Biochem, 2004). This pathway is reported to start with WSC1 or MID2 and to end with RLM1, SWI4 or SWI6.
  We will look for the lightest paths connecting MID2 to RLM1/SWI4/SWI6, WSC1 being absent in the demo network. The paths of first rank miss one step (ROM2) and replace MPK1 by KSS1, otherwise they are correct.
  ROM2 is missing in the demo network. With improved quality of the protein-protein interaction network used, path finding accuracy is expected to increase.");
  }
   echo("
   <br>
   ");
  	if(!$pipe){
  	 if(!$demo2){
        echo("<br>
        <a href='help.pathfinder.html#formats'><b>Input format</b></a>&nbsp;
        <select name='in_format'>
        <option selected value = 'flat'> tab-delimited format</option>
        <option value = 'gml'> GML format</option>
        </select>
        <br>");
        }else{
            echo("<input type='hidden' name='in_format' value='flat'></input>");
        }
  	}else{
  		info("Input graph format: ".$requested_in_format);
  		echo "<input type='hidden' NAME='pipe_in_format' VALUE='$requested_in_format'>";
  }
 if(!$demo2){
  echo("
   	<br>
   <table>
   	<tr><td><B><a href = 'help.pathfinder.html#directed'>Directed network</a></B></td> <td><input type='checkbox' name='directed' value='on'></input></td></tr>
    <tr><td><B><a href = 'help.pathfinder.html#server'>Store network on server</a></B></td> <td><input type='checkbox' name='store_graph' value='on'></input></td></tr>
   </table>
   <br>");
 }else{
    echo("<input type='hidden' name='directed' value='off'></input><input type='hidden' name='store_graph' value='off'></input>");
 }
 echo("
   <br>
   <hr>
   <h2>2. Seed nodes
   <a href='help.pathfinder.html#terminals'><img src='images/question-mark-button.jpg' width='15' alt='help'></a>
   </h2>
   <br>
   <br>");
   if($pipe && ($requested_in_format == 'gml' || $requested_in_format == 'GML')){
    echo("Please make sure to enter the identifiers (numbers) and not the labels of the nodes. If you are unsure about the node identifiers, check the input graph by clicking on 'Graph uploaded from the previous treatment'.");
   }
   echo("
  <font color='#006400'>Info: Matching of provided identifiers to identifiers of nodes in the network is case-insensitive.</font>
   <br>
   <br>
   Enter source and target nodes:
  <br>");
  echo("
  <br>
  <table>
     <tr><td>Source nodes&nbsp;&nbsp;</td> <td><input type='text' name='sources' value='$default_sources' size=100></input></td></tr>
     <tr><td>Target nodes</td><td><input type= 'text' name='targets' value='$default_targets' size=100></input></td></tr>
  </table>
  <br>
  <b>Or</b> upload a batch file:
   <br>
   <br>
  <input type='file' name='batch_file' size='45' />
  <br>
  <br>
   <hr>
  <h2>3. Path finding options
  <a href='help.pathfinder.html#options'><img src='images/question-mark-button.jpg' width='15' alt='help'></a>
  </h2>
  <br>
  <br>");
  if(!$demo2){
    echo("
    <table>
     <tr><td><B><a href = 'help.pathfinder.html#rank'>Rank</a></B></td>                   <td><input type = 'text' name='rank' value = '$default_rank' size = 10></input></td></tr>
     <tr><td><font color='#CC6600'>Warning: The edge weight is seen as edge cost, not as edge strength!</font></td></tr>
	 <tr><td><B><a href = 'help.pathfinder.html#weight'>Weighting scheme</a></B></td>     <td><select name='weight'>
                <option value = 'unit'>unit weight</option>
                <option value = 'none'>as given in input graph</option>
                <option selected value = 'con'>degree of nodes as weight</option>
               </select>
        </td></tr>
         <tr><td><B><a href = 'help.pathfinder.html#constraints'>Identifiers of nodes to exclude</a></B></td> <td><input type='text' NAME='exclude' VALUE='$default_exclude'></input></td></tr>
     <tr><td><B><a href = 'help.pathfinder.html#constraints'>Identifiers of nodes to include</a></B></td> <td><input type='text' NAME='include' VALUE='$default_include'></input></td></tr>
         <tr><td><B><a href = 'help.pathfinder.html#constraints'>Maximal path weight</a></B></td> <td><input type='text' NAME='maxWeight' VALUE='$default_maxWeight'></input></td></tr>
         <tr><td><B><a href = 'help.pathfinder.html#constraints'>Maximal path length</a></B></td> <td><input type='text' name='maxLength' value='$default_maxLength'></input></td></tr>
         <tr><td><B><a href = 'help.pathfinder.html#sconstraints'>Minimal path length</a></B></td> <td><input type='text' name='minLength' value='$default_minLength'></input></td></tr>
        </table>
    <br>");
  }else{
   echo("
    The weights have been already set.
    <br>
    <table>
     <tr><td><B><a href = 'help.pathfinder.html#rank'>Rank</a></B></td><td><input type = 'text' name='rank' value = '$default_rank' size = 10></input></td></tr>
     <tr><td><B><a href = 'help.pathfinder.html#constraints'>Identifiers of nodes to exclude</a></B></td> <td><input type='text' NAME='exclude' VALUE='$default_exclude'></input></td></tr>
     <tr><td><B><a href = 'help.pathfinder.html#constraints'>Identifiers of nodes to include</a></B></td> <td><input type='text' NAME='include' VALUE='$default_include'></input></td></tr>
   	 <tr><td><B><a href = 'help.pathfinder.html#constraints'>Maximal path weight</a></B></td> <td><input type='text' NAME='maxWeight' VALUE='$default_maxWeight'></input></td></tr>
     <tr><td><B><a href = 'help.pathfinder.html#constraints'>Maximal path length</a></B></td> <td><input type='text' name='maxLength' value='$default_maxLength'></input></td></tr>
     <tr><td><B><a href = 'help.pathfinder.html#sconstraints'>Minimal path length</a></B></td> <td><input type='text' name='minLength' value='$default_minLength'></input></td></tr>
   </table>
   <input type='hidden' name='weight' value=$default_weight></input>
    <br>
   ");
  }
  if($advanced){
    echo("
    <h4>Advanced Options
    <a href='help.pathfinder.html#advanced'><img src='images/question-mark-button.jpg' width='15' alt='help'></a>
    </h4>
    <br>
    <table>
          <tr><td><B><a href = 'help.pathfinder.html#exclusion'>Exclusion attribute</a></B></td> <td><input type='text' name='exAttrib' value='$default_exAttrib'></input></td></tr>
           <tr><td><B><a href = 'help.pathfinder.html#metabolic'>Metabolic</a></B></td> <td><input type='checkbox' name='metabolic' value='on'></input></td></tr>
         <tr><td><font color='#CC6600'>Warning: Backtracking is only available for metabolic graphs!</font></td></tr>
          <tr><td><B><a href = 'help.pathfinder.html#algorithm'>k shortest paths algorithm</a></B></td> <td>
          <select name='algorithm'>
            <option selected value = 'rea'>REA
            <option value = 'backtrack'>Backtracking
        </select>
  </td></tr>
    </table>
    <br>
    ");
  }else{
    echo "<input type='hidden' NAME='exAttrib' VALUE='$default_exAttrib'>";
    echo "<input type='hidden' NAME='metabolic' VALUE='$default_metabolic'>";
    echo "<input type='hidden' NAME='algorithm' VALUE='$default_algorithm'>";
  }
  echo("
  <br>
  <br>
  <hr>
  <h2>4. Output
  <a href='help.pathfinder.html#output'><img src='images/question-mark-button.jpg' width='15' alt='help'></a>
  </h2>
  <br><br>
 <a href='help.pathfinder.html#output'><b>Output format&nbsp;&nbsp;</a></b>
  <br><br>
  <table><tr>
  <td>Table: <input type='radio' name='outputchoice' value='pathsTable' checked='checked'></td>
  </tr><tr>
  <td>Graph: <input type='radio' name='outputchoice' value='outputgraph'></td>
  	<table>
  	<tr><td>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td>
  	    <td>Graph format:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  	         <select name='out_format'>
  			 <option selected value = 'flat'>tab-delimited format</option>
  			 <option value = 'GML'>GML format</option>
             </select>
  	    </td>
  	</tr>
  	<tr><td>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td>
  	    <td>Graph output type:
  	    <select name='outputType'>
        <option selected value = 'pathsGraphs'>each path as a separated component of a single graph</option>
        <option value = 'pathsUnion'>paths unified into one graph</option>
        <option value = 'pathsMarked'>input graph with paths highlighted</option>
         </select>
  	</td>
  	</tr>
  	</table>
  </tr>
  </table>
  <br><br>
  <table>
  <tr><td>
  <B><a href = 'help.pathfinder.html#email'>Email (optional)</a></B>
  </td><td>
  <input type='text' name='email' value='$default_email' size='30'></input>
  </td></tr>
  </table>
  <br>
  <hr>
  <table class='formbutton'>
  <TD><input type='submit' name='.submit' value='GO' /></TD>
  <TD><B><A HREF='pathfinder_form.php?demo=0'>RESET</A></B></TD>
  <TD><B><A HREF='pathfinder_form.php?demo=1'>DEMO1</A></B></TD>
  <TD><B><A HREF='pathfinder_form.php?demo2=1'>DEMO2</A></B></TD>
  <TD><B><A HREF='pathfinder_form.php?advanced=1'>ADVANCED</A></B></TD>
  </form>
  <TD><B><A HREF='help.pathfinder.html'>MANUAL</A></B></TD>
  <TD><B><A target = '_blank' HREF='".checkNeatTutorial("tutorials/neat_tutorial/Path_finding.html")."'>TUTORIAL</A></B></TD>
  <TD><B><A HREF='mailto:kfaust@ulb.ac.be'>MAIL</A></B></TD>
   <TD><B><A HREF='neat_webservices.php'>WSDL</A></B></TD>
  </TR></TABLE></ul></ul>");
?>
</body>
</html>