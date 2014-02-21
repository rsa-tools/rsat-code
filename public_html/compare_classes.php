<html>
<head>
   <title>Network Analysis Tools - convert-graph</title>
   <link rel="stylesheet" type="text/css" href = "main.css" media="screen">
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

  ## log file update
  UpdateLogFile("neat","","");  
  title('compare-classes - results');

  # File to store the commands
  $cmd_file = getTempFileName('commands_compare-classes', '.txt');
  $cmd_handle = fopen($cmd_file, 'a');

  ## Error status
  $error = 0;
  echo ($Qclass_file);
  echo ($Rclass_file);
  ## Get parameters
  if ($_FILES['Qclass_file']['name'] != "") {
    $Qclass_file = uploadFile('Qclass_file');
  } else if ($_REQUEST['pipe_Q_class_file'] != "") {
    $Qclass_file = $_REQUEST['pipe_Q_class_file'];
  }
  if ($_FILES['Rclass_file']['name'] != "") {
    $Rclass_file = uploadFile('Rclass_file');
  }
//   print_r($_FILES);
//   echo ($Qclass_file);
//   echo ($Rclass_file);
  
  $now = date("Ymd_His");
  $classesQ = $_REQUEST['classesQ'];
  $classesR = $_REQUEST['classesR'];
  $score_col = $_REQUEST['score_col'];
  $sort = "";

  ## self comparaison
  $self_compa = 0;
  $distinct = 0;
  $triangle = 0;
  if ($_REQUEST['self_compa'] == 'on') {
    $self_compa = 1;
  }
  if ($_REQUEST['distinct'] == 'on') {
    $distinct = 1;
  }
  if ($_REQUEST['triangle'] == 'on') {
    $triangle = 1;
  }

  ## RETURN FIELDS
  $return =  array("rank");

  ## Identify the column containing the sig score.
  ## This is required for the arc weight of the class-class graph.
  $sig_column = 6;
  if ($_REQUEST['occ'] == 'on') {
    array_push($return, "occ");
    $sig_column += 7;
  }
  if ($_REQUEST['freq'] == 'on') {
    array_push($return, "freq");
    $sig_column += 9;
  }
  if ($_REQUEST['proba'] == 'on') {
    array_push($return, "proba");
    if ($_REQUEST['sort'] == 'on') {
      $sort = "sig"; 
    }
  } else {
    $sig_column = 0; ## If proba is not computed, the graph will have no weight column
  }
  if ($_REQUEST['jac_sim'] == 'on') {
    array_push($return, "jac_sim");
  }  
  if ($_REQUEST['sor_sim'] == 'on') {
    array_push($return, "sor_sim");
  }  
  if ($_REQUEST['entropy'] == 'on') {
    array_push($return, "entropy");
  }  
  if ($_REQUEST['members'] == 'on') {
    array_push($return, "members");
  }  
  if ($_REQUEST['dotprod'] == 'on') {
    array_push($return, "dotprod");
  }
  $return = join(",", $return);

  # THRESHOLDS FIELDS
  $l_thr_val = array();
  $l_thr_field = array();
  $u_thr_val = array();
  $u_thr_field = array();  

  if ($_REQUEST['lth_q'] != "") {
    $thr = $_REQUEST['lth_q'];
    if (preg_match("/\d/", $thr)) {
      array_push($l_thr_field, "Q");
      array_push($l_thr_val, $thr);
    } else if ($thr != 'none') {
      warning("$thr is not a valid query class size lower threshold");
    }
  }
  
  if ($_REQUEST['lth_r'] != "") {
    $thr = $_REQUEST['lth_r'];
    if (preg_match("/\d/", $thr)) {
      array_push($l_thr_field, "R");
      array_push($l_thr_val, $thr);
    } else if ($thr != 'none'){
      warning("$thr is not a valid reference lower threshold");
    }
  }
  if ($_REQUEST['lth_qr'] != "") {
    $thr = $_REQUEST['lth_qr'];
    if (preg_match("/\d/", $thr)) {
      array_push($l_thr_field, "QR");
      array_push($l_thr_val, $thr);
    } else if ($thr != 'none'){
      warning("$thr is not a valid reference-query intersection lower threshold");
    }
  }
  if ($_REQUEST['lth_sig'] != "") {
    $thr = $_REQUEST['lth_sig'];
    if (preg_match("/\d/", $thr)) {
      array_push($l_thr_field, "sig");
      array_push($l_thr_val, $thr);
      ##      $sort = "sig"; ## JvH: I move this, because values can be sorted or not irrespective of the fact that a threshold was selected
    } else if ($thr != 'none'){
      warning("$thr is not a valid significance lower threshold");
    }
  }
  if ($_REQUEST['lth_pval'] != "") {
    $thr = $_REQUEST['lth_pval'];
    if (preg_match("/\d/", $thr)) {
      array_push($l_thr_field, "P_val");
      array_push($l_thr_val, $thr);
    } else if ($thr != 'none'){
      warning("$thr is not a valid P-value lower threshold");
    }
  }    
  if ($_REQUEST['lth_eval'] != "") {
    $thr = $_REQUEST['lth_eval'];
    if (preg_match("/\d/", $thr)) {
      array_push($l_thr_field, "E_val");
      array_push($l_thr_val, $thr);
    } else if ($thr != 'none'){
      warning("$thr is not a valid E-value lower threshold");
    }
  }    

  if ($_REQUEST['lth_jac'] != "") {
    $thr = $_REQUEST['lth_jac'];
    if (preg_match("/\d/", $thr)) {
      array_push($l_thr_field, "jac_sim");
      array_push($l_thr_val, $thr);
    } else if ($thr != 'none'){
      warning("$thr is not a valid Jaccard similarity lower threshold");
    }
  }      

  if ($_REQUEST['lth_sor'] != "") {
    $thr = $_REQUEST['lth_sor'];
    if (preg_match("/\d/", $thr)) {
      array_push($l_thr_field, "sor_sim");
      array_push($l_thr_val, $thr);
    } else if ($thr != 'none'){
      warning("$thr is not a valid Sorensen similarity lower threshold");
    }
  }      

  if ($_REQUEST['lth_mi'] != "") {
    $thr = $_REQUEST['lth_mi'];
    if (preg_match("/\d/", $thr)) {
      array_push($l_thr_field, "I(Q,R)");
      array_push($l_thr_val, $thr);
    } else if ($thr != 'none'){
      warning("$thr is not a valid mutual information lower threshold");
    }
  }
   
 if ($_REQUEST['uth_q'] != "") {
    $thr = $_REQUEST['uth_q'];
    if (preg_match("/\d/", $thr)) {
      array_push($u_thr_field, "Q");
      array_push($u_thr_val, $thr);
    } else if ($thr != 'none'){
      warning("$thr is not a valid query size upper threshold");
    }
  }
  if ($_REQUEST['uth_r'] != "") {
    $thr = $_REQUEST['uth_r'];
    if (preg_match("/\d/", $thr)) {
      array_push($u_thr_field, "R");
      array_push($u_thr_val, $thr);
    } else if ($thr != 'none'){
      warning("$thr is not a valid reference size upper threshold");
    }
  }
  if ($_REQUEST['uth_qr'] != "") {
    $thr = $_REQUEST['uth_qr'];
    if (preg_match("/\d/", $thr)) {
      array_push($u_thr_field, "QR");
      array_push($u_thr_val, $thr);
    } else if ($thr != 'none'){
      warning("$thr is not a valid reference-query intersection size upper threshold");
    }
  }
  if ($_REQUEST['uth_sig'] != "") {
    $thr = $_REQUEST['uth_sig'];
    if (preg_match("/\d/", $thr)) {
      array_push($u_thr_field, "sig");
      array_push($u_thr_val, $thr);
    } else if ($thr != 'none'){
      warning("$thr is not a valid significance upper threshold");
    }
  }
  if ($_REQUEST['uth_pval'] != "") {
    $thr = $_REQUEST['uth_pval'];
    if (preg_match("/\d/", $thr)) {
      array_push($u_thr_field, "P_val");
      array_push($u_thr_val, $thr);
    } else if ($thr != 'none'){
      warning("$thr is not a valid P-value upper threshold");
    }
  }    
  if ($_REQUEST['uth_eval'] != "") {
    $thr = $_REQUEST['uth_eval'];
    if (preg_match("/\d/", $thr)) {
      array_push($l_thr_field, "E_val");
      array_push($l_thr_val, $thr);
    } else if ($thr != 'none') {
      warning("$thr is not a valid E-value upper threshold");
    }
  }    

  if ($_REQUEST['uth_jac'] != "") {
    $thr = $_REQUEST['uth_jac'];
    if (preg_match("/\d/", $thr)) {
      array_push($u_thr_field, "jac_sim");
      array_push($u_thr_val, $thr);
    } else if ($thr != 'none'){
      warning("$thr is not a valid Jaccard similarity index threshold");
    }
  }      

  if ($_REQUEST['uth_sor'] != "") {
    $thr = $_REQUEST['uth_sor'];
    if (preg_match("/\d/", $thr)) {
      array_push($u_thr_field, "sor_sim");
      array_push($u_thr_val, $thr);
    } else if ($thr != 'none'){
      warning("$thr is not a valid Sorensen similarity index threshold");
    }
  }      

  if ($_REQUEST['uth_mi'] != "") {
    $thr = $_REQUEST['uth_mi'];
    if (preg_match("/\d/", $thr)) {
      array_push($u_thr_field, "I(Q,R)");
      array_push($u_thr_val, $thr);
    } else if ($thr != 'none'){
      warning("$thr is not a valid mutual information index threshold");
    }
  }  
   
  $l_thr_val = join(",", $l_thr_val);
  $u_thr_val = join(",", $u_thr_val);
  $l_thr_field = join(":", $l_thr_field);
  $u_thr_field = join(":", $u_thr_field);

  # output format
  $output_format = $_REQUEST['out_format'];
  $matrix = 0;
  
  if ($output_format == "matrix") {
    $matrix = $_REQUEST['matrix_field'];
    $l_thr_val = "";
    $u_thr_val = "";
    $l_thr_field = "";
    $u_thr_field = "";
    $sort = 0;
    $return = 0;
  } 

  # Place the query  file in variables
  if ($Qclass_file != "" && $classesQ == "") {
    $classesQ = storeFile($Qclass_file);
    $classesQ = trim_text($classesQ);
  }    
  # Place the reference file in variables
  if ($Rclass_file != "" && $classesR == "") {
    $classesR = storeFile($Rclass_file);
    $classesR = trim_text($classesR);
  }    
  
  ## classes R -> empty & self compa != 1 -> error
  if ($classesR == '' && $self_compa == 0) {
    $error = 1;
    error("You must submit a set of reference classes");
  }
  
  ## classes Q -> empty : error
  if ($classesQ == '') {
    $error = 1;
    error("You must submit a set of query classes");
  }  

  ## If the score column is given for 2 different files : error.
  if ($_REQUEST['dotprod'] == 'on' && $score_col == "") {
    $error = 1;
    error("You must submit a score column to compute the dot product ");
  }
  ## If -i option and reference classes are given -> error
  if ($self_compa && $classesR != "") {
    $error = 1;
    error("You must not specify reference classes when query classes are compared to themselves");
  }
  
  ## If -i option and query classes -> place the classes in a new variable
  if ($self_compa) {
    $input_classes = $classesQ;
    $classesQ = "";
  }  
  
   if (!$error) { 
     $cc_parameters = array( 
       "request" => array (
         "ref_classes"=>$classesR,
         "query_classes"=>$classesQ,
         "input_classes"=>$input_classes,
         "return_fields"=>$return,
         "score_column"=>$score_col,
         "upper_threshold_field"=>$u_thr_field,
         "lower_threshold_field"=>$l_thr_field,
         "lower_threshold_value"=>$l_thr_val,
         "upper_threshold_value"=>$u_thr_val,
         "sort"=>$sort,
         "distinct"=>$distinct,
         "triangle"=>$triangle,
         "matrix"=>$matrix
       )

     );       
         
    # Info message
    info("Results will appear below");
    echo"<hr>\n";
    hourglass("on");
    # Open the SOAP client
    $soap_client = new SoapClient(
                       $neat_wsdl,
                           array(
                                 'trace' => 1,
                                 'soap_version' => SOAP_1_1,
                                 'style' => SOAP_DOCUMENT,
                                 'encoding' => SOAP_LITERAL
                                 )
                           );


    # Execute the command and catch the errors
    try {
      $cc_echoed = $soap_client->compare_classes($cc_parameters);
      $soap_error = 0;
    } catch (Exception $soap_exception) {
      echo ("<pre>");
      echo "Error : \n\t",  $soap_exception->getMessage(), "\n";
      echo ("</pre>");
      $soap_error = 1;
      exit(1);
    }  
    $cc_response =  $cc_echoed->response;
    $cc_command = $cc_response->command;
    $cc_server = $cc_response->server;
    $cc_client = $cc_response->client;
    $cc_resultURL = rsat_path_to_url($cc_server);
    store_command($cc_command, "Class conversion", $cmd_handle);
    $URL['Class comparison (tab)'] = $cc_resultURL;

    ## Text-to-html
    $cc_server = rtrim ($cc_server);
    $cc_file = storeFile($cc_server, "html");
    $tth_parameters = array( 
        "request" => array(
        "inputfile"=>$cc_file,
        "chunk"=>1000,
      )
    );
    $tth_echoed = $soap_client->text_to_html($tth_parameters);
    $tth_response =  $tth_echoed->response;
    $tth_command = $tth_response->command;
    $tth_server = $tth_response->server;
    $tth_client = $tth_response->client;
    store_command($tth_command, "text-to-html", $cmd_handle);
    $URL['Class comparison (html)'] = rsat_path_to_url($tth_server);
    
    
    hourglass("off");

    ## Close command handle
    fclose($cmd_handle);
    $URL['Server commands'] = rsat_path_to_url($cmd_file);

    ## DISPLAY THE RESULT
    print_url_table($URL);

    ## PRESENT BUTTONS TO SEND RESULT TO OTHER TOOLS
    if ($output_format == 'matrix') {
      echo "
      <table class = 'nextstep'>
        <tr>
          <th colspan = 3>
            Next step
          </th>
        </tr>
        <tr>
          <td>
            <form method='post' action='contingency_stats_form.php'>
              <input type='hidden' name='pipe' value='1'>
              <input type='hidden' name='matrix_file' value='$cc_server'>
              <INPUT type='submit' value='Contingency-table statistics'>
            </form>
          </td>
          <td>
            <form method='post' action='draw_heatmap_form.php'>
              <input type='hidden' name='pipe' value='1'>
              <input type='hidden' name='table_file' value='$cc_server'>
              <INPUT type='submit' value='Display the matrix as heatmap'>
            </form>
          </td>
        </tr>
      </table>";

    } else {

      $server = $cc_server; ## Use the tab-delimited file as input for the next steps

      ## Open a table for the "next step" buttons
      echo("<table class = 'nextstep'>
    <tr>
      <th colspan = 3>
        Next step
      </th>
    </tr>
    <tr>");

      ## Send to display-graph
	echo ("
      <td>
        <form method='post' action='display_graph_form.php'>
          <input type='hidden' name='pipe' value='1'>
          <input type='hidden' name='graph_file' value='$server'>
          <input type='hidden' name='graph_format' value='tab'>
          <input type='hidden' name='scol' value='1'>
          <input type='hidden' name='tcol' value='2'>
          <input type='submit' value='Display the graph'>
        </form>
      </td>");

      ## Send to compare-graphs
      echo ("<td>
        <form method='post' action='compare_graphs_form.php'>
          <input type='hidden' name='pipe' value='1'>
          <input type='hidden' name='graph_file' value='$server'>
          <input type='hidden' name='graph_format' value='tab'>
          <input type='hidden' name='scol' value='1'>
          <input type='hidden' name='tcol' value='2'>
          <input type='submit' value='Compare this graph to another one'>
        </form>
      </td>");

      ## Send to convert-graph
      echo ("<td>
        <form method='post' action='convert_graph_form.php'>
          <input type='hidden' name='pipe' value='1'>
          <input type='hidden' name='graph_file' value='$server'>
          <input type='hidden' name='graph_format' value='tab'>
          <input type='hidden' name='scol' value='1'>
          <input type='hidden' name='tcol' value='2'>
          <input type='submit' value='Convert table to graph'>
        </form>
      </td>");

      ## New row
      echo("</tr><tr>");

      ## Send to graph-get-clusters
      echo ("<td>
        <form method='post' action='graph_get_clusters_form.php'>
          <input type='hidden' name='pipe' value='1'>
          <input type='hidden' name='graph_file' value='$server'>
          <input type='hidden' name='graph_format' value='tab'>
          <input type='hidden' name='scol' value='1'>
          <input type='hidden' name='tcol' value='2'>
          <input type='submit' value='Map clusters or extract a subnetwork'>
        </form>
      </td>");

      ## Send to graph-topology
      echo ("<td>
        <form method='post' action='graph_topology_form.php'>
          <input type='hidden' name='pipe' value='1'>
          <input type='hidden' name='graph_file' value='$server'>
          <input type='hidden' name='graph_format' value='tab'>
          <input type='hidden' name='scol' value='1'>
          <input type='hidden' name='tcol' value='2'>
          <input type='submit' value='Node topology statistics'>
        </form>
      </td>");

      ## Send to graph-neighbours
      echo ("<td>
        <form method='post' action='graph_neighbours_form.php'>
          <input type='hidden' name='pipe' value='1'>
          <input type='hidden' name='graph_file' value='$server'>
          <input type='hidden' name='graph_format' value='tab'>
          <input type='hidden' name='scol' value='1'>
          <input type='hidden' name='tcol' value='2'>
          <input type='submit' value='Neighbourhood analysis'>
        </form>
      </td> ");

      ## New table row
      echo("</tr><tr>");

      ## Send to MCL
      echo ("<td>
        <form method='post' action='mcl_form.php'>
          <input type='hidden' name='pipe' value='1'>
          <input type='hidden' name='graph_file' value='$server'>
          <input type='hidden' name='graph_format' value='tab'>
          <input type='hidden' name='scol' value='1'>
          <input type='hidden' name='tcol' value='2'>
          <input type='submit' value='MCL Graph clustering'>
        </form>
      </td>");

      ## Send to RNSC
      echo ("<td>
        <form method='post' action='rnsc_form.php'>
          <input type='hidden' name='pipe' value='1'>
          <input type='hidden' name='graph_file' value='$server'>
          <input type='hidden' name='graph_format' value='tab'>
          <input type='hidden' name='scol' value='1'>
          <input type='hidden' name='tcol' value='2'>
          <input type='submit' value='RNSC Graph clustering'>
        </form>
      </td>");


      ## Send to alter-graph
      echo ("<td>
        <form method='post' action='alter_graph_form.php'>
          <input type='hidden' name='pipe' value='1'>
          <input type='hidden' name='graph_file' value='$server'>
          <input type='hidden' name='graph_format' value='tab'>
          <input type='hidden' name='scol' value='1'>
          <input type='hidden' name='tcol' value='2'>
          <input type='submit' value='Graph alteration'>
        </form>
      </td>");

      ## New table row
      echo ("</tr><tr>");

      ## Send to path finder
      echo ("<td>
        <form method='post' action='pathfinder_form.php'>
          <input type='hidden' name='pipe' value='1'>
          <input type='hidden' name='graph_file' value='$server'>
          <input type='hidden' name='graph_format' value='tab'>
          <input type='hidden' name='scol' value='1'>
          <input type='hidden' name='tcol' value='2'>
          <input type='submit' value='Path Finding'>
        </form>
      </td>");

      ## Send to visant
      echo("<td>
          <form method='post' action='visant.php'>
          <input type='hidden' name='pipe' value='1'>
          <input type='hidden' name='visant_graph_file' value='$server'>
          <input type='hidden' name='graph_format' value='tab'>
          <input type='hidden' name='scol' value='1'>
          <input type='hidden' name='tcol' value='2'>
          <input type='submit' value='Load in VisANT'>
          </form>
        </td>");

      ## Close the table
      echo ("</tr></table>");
    }
         
  }
?>
