<html>
<head>
<title>RSAT - Pathway-Extractor</title>
<link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">

</head>
   <body class="results" onload="javascript:document.getElementById('hourglass').hidden=true"> 
   <!--div id=hourglass><img src="../images/animated_hourglass.gif" alt="Please Wait!!"></img></div>-->
<?php
// Load RSAT configuration
   require('functions.php');

// print_r($properties);
UpdateLogFile("rsat","","");

// Import variables with prefix fs_ from form
import_request_variables('P','fs_');

//print_r($_POST);
// print_r($_FILES);

// Initialize variables
$basedir=$properties['RSAT']."/public_html/data/metabolic_networks";
#$outputdir = $properties['RSAT']."/public_html/tmp/";
#$outputdir =  $properties['RSAT']."/public_html/".getTempFileName('pathway-extractor', '.txt');
$outputdir =  getTempFileName('pathway-extractor');
$URL['Output directory'] = rsat_path_to_url($outputdir);
$outdir_created = mkdir($outputdir, 0775, TRUE);
if (!$outdir_created) {
  error("Could not create output directory");
 }

#$outputurl = $properties['rsat_ws_tmp'];
$outputurl = $rsat_path_to_url[$outputdir];

$CLASSPATH=$properties['RSAT']."/java/lib/NeAT_javatools.jar";
putenv("CLASSPATH=".$CLASSPATH);
//info("CLASSPATH=".str_replace($properties['RSAT'], '$RSAT', $CLASSPATH));

$cmda = $properties['RSAT'].'/perl-scripts/pathway-extractor';
$cmdb = $properties['RSAT'].'/perl-scripts/process-pathwayextractor-output';
$groupdescriptor = uniqid('PathwayExtractor');
$argumenta = " -v 2";
//$argumenta .= " -d ";

$argumentb = $argumenta;
if ($fs_directedgraph =="directedgraph") {
  $argumenta .= " -d ";
  $argumentb .= " -d ";
}else{
  $argumenta .= " -J -Q ";
}

$argumenta .= " -prefix ". $groupdescriptor;
#$exp_bed_file = "/^[\w\-\+\s,\.\#; \/]+$/";

////////////////////////////////////////////////////////////////
//Print <h3>
echo "<H3><a href='".$properties['rsat_www']."'>RSAT</a> - Pathway-Extractor - results</H3>";

////////////////////////////////////////////////////////////////
// Check arguments
$errors = 0;

// Check that gnn has been specified
if ($fs_gnn == "none" or $fs_gnn == "" ) {
  error( "You forgot to select the genome source for the mapping.");
  $errors=1;
} else {
	$GERdir = "$basedir/GER_files/$fs_gnn";
   $argumenta .= " -gnn ".$GERdir."/".$fs_gnn."-gene_ec.tab";
   $argumentb .= " -gnn ".$GERdir."/".$fs_gnn."-gene_ec.tab";
 }
// Check that gnn has been specified
if ($fs_network == "none" or $fs_network == "" ) {
  error( "You forgot to select network.");
  $errors=1;
} else {
  $neworkfilepattern= $fs_network;
  $neworkfilepattern = eregi_replace("_v.*","",$fs_network); 
  $neworkfilepattern = "$basedir/networks/$neworkfilepattern/$fs_network";
  $networknodenames = $neworkfilepattern."-metab-node_names.tab";
  $networkfile = $neworkfilepattern."-metab-network.tab";
  $argumenta .= " -nnn ".$networknodenames ." -g " . $networkfile;
  $argumentb .= " -nnn ".$networknodenames;
}
if ($fs_algorithm){
 	$argumenta .= " -a $fs_algorithm";
}
if ($fs_preproc){
 	$argumenta .= " -P";
}
if ($fs_postproc){
 	$argumenta .= " -C";
}
if ($fs_iter){
 	$argumenta .= " -I ".$fs_iter;
}
if ($fs_percentage){
 	$argumenta .= " -x ".$fs_percentage;
}
if ($fs_kwalksweight){
 	$argumenta .= " -k true";
}
$argumenta .= " -T pathsUnion -u";

//Check syntax of email address (ensure the texte netered in email box was an email address)
if($fs_output =="email") {
  if (!preg_match("#^[\wàáâãäåæçèéêëìíîïðñòóôõöøùúûüý._\-]+@([a-z]+.)+[a-z]{2,4}$#", $fs_user_email)) {
     error( "Email not valid");
     $errors=1;
  }
 }

////////////////////////////////////////////////////////////////
// Check that the BED has been specified once and only once

// Bed data pasted in text area
if (count($_POST["seeds"]) == 0) {
  error( "No seeds!");
  $errors=1;
 }else {
 	
  //$mystring = preg_replace( "#\r\n|\r|\t#", "\n", $_POST["seeds"] );
 	$mystring =""; 
  foreach( $_POST["seeds"] as $key=>$val){
//  		info($key."=>".$val);
  		$mystring.=str_replace(',', "\n", $val)."\n";
  }
  
//  info("SEEDS:".$mystring);
  #$array_line = explode("\n",$fs_seeds);
  $seed_file = $outputdir."/".$groupdescriptor."_seeds.tab";
  $file = fopen($seed_file, "w");
    fwrite($file, $mystring);
    fclose($file);
    $argumenta .= " -i $seed_file";
 }
 

$argumenta .= " -o $outputdir";

///////////////////////////////////////////
// Run commands
if ($errors == 0) { 
  // Add arguments to the command
  $argumentb .= " -o ".$outputdir;
  $cmdb .= $argumentb;		

  // Announce job starting
  $msg = "Starting job.";
  if ($fs_output =="email")  {
    $msg .= " After job completion, email will be sent to ".$fs_user_email;
  }
  info($msg);
  
  flush(); 

  // Report the command 
  $cmda .= $argumenta;
  info("Command : ".str_replace($properties['RSAT'], '$RSAT', $cmda));
  // info("Command : ".$cmda);

  //  Run the command
  exec($cmda, $output);

  foreach ($output as $line) {
    if ($line[0] == ";") {
      $info .= $line."<br/>\n";
    } else  if (eregi('^OUTPUTFILE =',$line)){
        $cmdb .= ' -i '. substr( $line , 12);
	info("Command : ".str_replace($properties['RSAT'], '$RSAT', $cmdb));
	// info("Command : ".$cmdb);
    	exec($cmdb, $outputb);
    } else {
    	echo $line;
    }
  }

  echo "<hr>";
  //  $result_files = new DirectoryIterator("glob:".$outputdir."/".$groupdescriptor."*");

  $result_files = glob($outputdir."/".$groupdescriptor."*");
  sort($result_files);
  // print_r($result_files);


  /*  $returned = "<table class='resultlink'>\n";
   $returned .= "<tr><th colspan='2'>Result file(s)</th></tr>\n";*/
  
  foreach($result_files as $f) {
    $returned .= "<tr><td>";
    if (preg_match("#png$#",$f)){
      // echo "<img width=800  title=\"inferedpathway\"src=\"".$outputurl."/".$f."\"><br/> <hr> <br/>\n";
      echo '<img width=800  title="inferedpathway" src="'.rsat_path_to_url($f).'"><br/> <hr> <br/>'."\n";
      $file_type = "Extracted pathway image file";
      //      $returned .= $file_type.": ";
    }elseif (preg_match("#pred_pathways.txt$#",$f)){
      $file_type = "Extracted pathway graph file";
      //$returned .= $file_type.": ";
    }elseif (preg_match("#pred_pathways_annot.dot$#",$f)){
      $file_type = "Extracted annotated pathway graph (dot format)";
      //$returned .= $file_type."./: ";
    }elseif (preg_match("#pred_pathways_annot.txt$#",$f)){
      $file_type = "Extracted pathway annotated graph";
    }elseif (preg_match("#pred_pathways_seeds_converted.txt$#",$f)){
      $file_type = "Mapped seeds";
    }elseif (preg_match("#seeds.tab$#",$f)){
      $file_type = "Seeds file: ";
    }else {
      info("Unknown file:". $f );
    }
    //$returned .= $file_type."</td>";
    $URL[$file_type] = rsat_path_to_url($f);
    //$returned .=  "</td><td><a href=\"".$outputurl."/".$f."\" >".$f."</a></td></tr>\n";
  }

  //$returned .= "</table>\n";
  //echo	$returned;
  //  echo "<hr>";
  // Display the command
  /*
  $cmd_report = str_replace($properties['RSAT'], '$RSAT', $cmd);
  info("Command : ".$cmd_report);
  echo "<hr>";
  */
  //display log file

  /*$info = "";
  $warning = "";
  $logs = file_get_contents(rsat_path_to_url($error_file));
  $logs = explode("\n", $logs);

  foreach ($logs as $line) {
    if ($line[0] == ";") {
      $info .= $line."<br/>\n";
    } else  {
       $warning .= $line."<br/>\n";
    } 
  }

  if ($info != "") {
    info($info);
    echo "<hr>";
  }    
  
  if ($warning != "<br/>\n") {
    warning($warning);
    echo "<hr>";
  }  

  */

  // Display the result
  print_url_table($URL);
    

  ################################################################
  // Send email with notification of task completion 
  if ($fs_output =="email")  {
    //    echo "Email outpout not available for now";
    // Parammeters for sending the mail
    $to = $fs_user_email;
    $subject = "[RSAT] Pathway-extractor result ".$now."_".$suffix;
    
    // Store the URL table in a variable
    $html_mail = 0; // Boolean variable indicating whether HTML format is supported in email
    $headers = ""; // Header (specifying Mime types)
    if ($html_mail) {
      $msg = $returned;
      
      $headers .= 'Mime-Version: 1.0'."\r\n";
      $headers .= 'Content-type: text/html; charset=utf-8'."\r\n";
      $headers .= "\r\n";
    } else {
      $msg = "fetch-sequences result\n\n";
      $msg .= "Result files:\n";
      foreach ($URL as $key => $value) {
	  $msg .= "\t".$key."\t".$value."\n";
      }
    }

    
    // Sending mail
    $smtp = $properties["smtp"];

    // Check that the SMTP was specificed in the property file of the server
    if ($smtp == "") {
      error("SMTP server is not specified in the RSAT config file. Please contact system administrator.");
    } else {
      info("Sending mail via SMTP ".$smtp);
      ini_set ( "SMTP", $smtp); 
      $mail_sent = mail($to, $subject, $msg, $headers);
      if ($mail_sent) {
	info("Job done, email sent to ".$to);
      } else {
	error("Notification mail could not be sent.\n\n.".$msg);
      }
    }  
  }  

  //Display pipe
  /*
  echo ('   <table class = "nextstep">
    <tr><th>next step</th></tr>
    <tr><td align=center>
	    <form method="post" action="peak-motifs_form.cgi">
	    <input type="hidden" name="title" value="'.str_replace('.bed', '', end(explode('/',$bed_file))).'">
	    <input type="hidden" name="sequence_url1" value="'.rsat_path_to_url($output_file).'">
	    <input type="submit" value="peak-motif">
	    </form></td>
	  </tr>
	  </table>');  
	  */
}	

?>
 
  </body>
</html>
