<html>
<head>
   <title>GrA-tools - convert-graph</title>
   <link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
</head>
<body class="results">
<?php 
  require ('functions.php');
  title('compare-graphs - results');
  # Error status
  $error = 0;
  # Get parameters
  if ($_FILES['Qclass_file']['name'] != "") {
    $Qclass_file = uploadFile('Qclass_file');
  } else if ($_REQUEST['pipe_Q_class_file'] != "") {
    $Qclass_file = $_REQUEST['pipe_Q_class_file'];
  }
  if ($_FILES['Rclass_file']['name'] != "") {
    $Rclass_file = uploadFile('Rclass_file');
  }
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
  if ($_REQUEST['occ'] == 'on') {
    array_push($return, "occ");
  }
  if ($_REQUEST['freq'] == 'on') {
    array_push($return, "freq");
  }
  if ($_REQUEST['proba'] == 'on') {
    array_push($return, "proba");
  }  
  if ($_REQUEST['jac'] == 'on') {
    array_push($return, "jac");
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
    } else {
      warning("$thr is not a valid query lower threshold");
    }
  }
  if ($_REQUEST['lth_r'] != "") {
    $thr = $_REQUEST['lth_r'];
    if (preg_match("/\d/", $thr)) {
      array_push($l_thr_field, "R");
      array_push($l_thr_val, $thr);
    } else {
      warning("$thr is not a valid reference lower threshold");
    }
  }
  if ($_REQUEST['lth_qr'] != "") {
    $thr = $_REQUEST['lth_qr'];
    if (preg_match("/\d/", $thr)) {
      array_push($l_thr_field, "QR");
      array_push($l_thr_val, $thr);
    } else {
      warning("$thr is not a valid reference-query intersection lower threshold");
    }
  }
  if ($_REQUEST['lth_sig'] != "") {
    $thr = $_REQUEST['lth_sig'];
    if (preg_match("/\d/", $thr)) {
      array_push($l_thr_field, "sig");
      array_push($l_thr_val, $thr);
      $sort = "sig";
    } else {
      warning("$thr is not a valid significance lower threshold");
    }
  }
  if ($_REQUEST['lth_pval'] != "") {
    $thr = $_REQUEST['lth_pval'];
    if (preg_match("/\d/", $thr)) {
      array_push($l_thr_field, "P_val");
      array_push($l_thr_val, $thr);
    } else {
      warning("$thr is not a valid P-value lower threshold");
    }
  }    
  if ($_REQUEST['lth_eval'] != "") {
    $thr = $_REQUEST['lth_eval'];
    if (preg_match("/\d/", $thr)) {
      array_push($l_thr_field, "E_val");
      array_push($l_thr_val, $thr);
    } else {
      warning("$thr is not a valid E-value lower threshold");
    }
  }    
  if ($_REQUEST['lth_jac'] != "") {
    $thr = $_REQUEST['lth_jac'];
    if (preg_match("/\d/", $thr)) {
      array_push($l_thr_field, "jac_sim");
      array_push($l_thr_val, $thr);
    } else {
      warning("$thr is not a valid Jaccard similarity lower threshold");
    }
  }      
  if ($_REQUEST['lth_mi'] != "") {
    $thr = $_REQUEST['lth_mi'];
    if (preg_match("/\d/", $thr)) {
      array_push($l_thr_field, "I(Q,R)");
      array_push($l_thr_val, $thr);
    } else {
      warning("$thr is not a valid mutual information lower threshold");
    }
  }
   
 if ($_REQUEST['uth_q'] != "") {
    $thr = $_REQUEST['uth_q'];
    if (preg_match("/\d/", $thr)) {
      array_push($u_thr_field, "Q");
      array_push($u_thr_val, $thr);
    } else {
      warning("$thr is not a valid query size upper threshold");
    }
  }
  if ($_REQUEST['uth_r'] != "") {
    $thr = $_REQUEST['uth_r'];
    if (preg_match("/\d/", $thr)) {
      array_push($u_thr_field, "R");
      array_push($u_thr_val, $thr);
    } else {
      warning("$thr is not a valid reference size upper threshold");
    }
  }
  if ($_REQUEST['uth_qr'] != "") {
    $thr = $_REQUEST['uth_qr'];
    if (preg_match("/\d/", $thr)) {
      array_push($u_thr_field, "QR");
      array_push($u_thr_val, $thr);
    } else {
      warning("$thr is not a valid reference-query intersection size upper threshold");
    }
  }
  if ($_REQUEST['uth_sig'] != "") {
    $thr = $_REQUEST['uth_sig'];
    if (preg_match("/\d/", $thr)) {
      array_push($u_thr_field, "sig");
      array_push($u_thr_val, $thr);
    } else {
      warning("$thr is not a valid significance upper threshold");
    }
  }
  if ($_REQUEST['uth_pval'] != "") {
    $thr = $_REQUEST['uth_pval'];
    if (preg_match("/\d/", $thr)) {
      array_push($u_thr_field, "P_val");
      array_push($u_thr_val, $thr);
    } else {
      warning("$thr is not a valid P-value upper threshold");
    }
  }    
  if ($_REQUEST['uth_eval'] != "") {
    $thr = $_REQUEST['uth_eval'];
    if (preg_match("/\d/", $thr)) {
      array_push($l_thr_field, "E_val");
      array_push($l_thr_val, $thr);
    } else {
      warning("$thr is not a valid E-value upper threshold");
    }
  }    
  if ($_REQUEST['uth_jac'] != "") {
    $thr = $_REQUEST['uth_jac'];
    if (preg_match("/\d/", $thr)) {
      array_push($u_thr_field, "jac_sim");
      array_push($u_thr_val, $thr);
    } else {
      warning("$thr is not a valid Jaccard similarity index threshold");
    }
  }      
  if ($_REQUEST['uth_mi'] != "") {
    $thr = $_REQUEST['uth_mi'];
    if (preg_match("/\d/", $thr)) {
      array_push($u_thr_field, "I(Q,R)");
      array_push($u_thr_val, $thr);
    } else {
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
    $l_thr_val = 0;
    $u_thr_val = 0;
    $l_thr_field = 0;
    $u_thr_field = 0;
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
    error("You must not specify reference classes when query classes are compare to themselves");
  }
  
  ## If -i option and query classes -> place the classes in a new variable
  if ($self_compa) {
    $input_classes = $classesQ;
    $classesQ = "";
  }  
  
   if (!$error) { 
     $parameters = array( 
       "request" => array (
         "ref_classes"=>$classesQ,
         "query_classes"=>$classesR,
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
  
    # Open the SOAP client
    $client = new SoapClient(
                       "$WWW_RSA"."/web_services/RSATWS.wsdl",
                           array(
                                 'trace' => 1,
                                 'soap_version' => SOAP_1_1,
                                 'style' => SOAP_DOCUMENT,
                                 'encoding' => SOAP_LITERAL
                                 )
                           );
    # Execute the command
    echo "<pre>";
    $echoed = $client->compare_classes($parameters);
    echo "</pre>"; 
    $response =  $echoed->response;
    $command = $response->command;
    $server = $response->server;
    $client = $response->client;
    $server = rtrim ($server);
    $temp_file = explode('/',$server);
    $temp_file = end($temp_file);
    $resultURL = $WWW_RSA."/tmp/".$temp_file;
    # The comment file has the same name as the
    # result file with ".comments" at the end of the string.
    # Comments
    echo "<pre>";
    echo "</pre>";
    # Display the results
    echo "The results is available at the following URL ";
    echo "<a href = '$resultURL'>$resultURL</a>"; 
    echo "<hr>\n";
 
  }
?>