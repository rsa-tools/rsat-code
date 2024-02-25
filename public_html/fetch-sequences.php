<html>
<head>
<title>RSAT - fetch-sequences</title>
<link rel="stylesheet" type="text/css" href = "css/main_grat.css" media="screen">
   </head>
   <body class="results"> 

<?php
// Load RSAT configuration
   require('functions.php');
    
    include('menu.php');
    //printMenu();
// print_r($properties);
UpdateLogFile("rsat","","");

// print_r($_POST);
// print_r($_FILES);

// Initialize variables
$cmd = $properties['RSAT'].'/perl-scripts/fetch-sequences';
$argument = " -v 1";
$exp_bed_file = "/^[\w\-\+\s,\.\#; \/]+$/";

//Fill buffer
echo str_repeat(" ", 1024), "\n";

////////////////////////////////////////////////////////////////
//Print <h3>
echo "<H3><a href='".$properties['rsat_www']."'>RSAT</a> - fetch-sequences - results</H3>";

////////////////////////////////////////////////////////////////
// Check arguments
$errors = false;

// Check that genome has been specified
$genome = $_REQUEST['genome'];

if ($genome == "none" or $genome == "" ) {
  error( "You forgot to specify the genome.");
  $errors = true;
 } else {
  $argument .= " -genome $genome";
 }


//Check syntax of email address (ensure the texte netered in email box was an email address)
$output = $_REQUEST['output'];
$user_email = $_REQUEST['user_email'];
if($output =="email") {
  if (!preg_match("#^[^@\.]+(\.[^@]+)*@[a-zA-Z0-9-]+(\.[a-zA-Z0-9-]+)*(\.[a-zA-Z]{2,4})$#", $user_email)) {
     error( "Email not valid");
     $errors=true;
  }
 }

////////////////////////////////////////////////////////////////
// Check that the BED has been specified once and only once
$bed_specifications = 0;

// Bed data pasted in text area
$bed = $_REQUEST['bed'];
if ($bed != "") {
  $bed_specifications++;
 }
// Local bed file on client machine
if ($_FILES["bedfile"]['name'] != "") {
  $bed_specifications++;
 }
// Bed specified by URL of remote server
$sequence_url = $_REQUEST['sequence_url'];
if ( $sequence_url != "") {
  $bed_specifications++;
 }

// Report error if no BED has been specified
if ($bed_specifications == 0) {
  error("You forgot to specify the genomic coordinates (BED file)");
  $errors=true;
 } else if ($bed_specifications > 1) {
  error("The BED file can be specified only in 1 way (paste, upload or URL of external server)");
  $errors = true;
 }

// Header format
$header = $_REQUEST['header'];
if ($header == "galaxy") {
  $argument .= " -header_format galaxy";
 }

// Downstream extension
$downstr_ext = $_REQUEST['downstr_ext'];
if (($downstr_ext !="") &&
    ($downstr_ext != '0')) {
  if (!is_numeric($downstr_ext)) {
    error ($downstr_extr." Invalid value for downstream extension: should be an Integrer.");
    $errors = true;
  } else {
     $argument .= " -downstr_ext $downstr_ext";
  }
}

// Upstream extension
$upstr_ext = $_REQUEST['upstr_ext'];
if (($upstr_ext !="") &&
    ($upstr_ext != '0')) {
  if (!is_numeric($upstr_ext)) {
    error($upstr_ext." is not valid value for upstream extension : should be an Integrer.");
    $errors = true;
  }  else {
    $argument .= " -upstr_ext $upstr_ext";
  }
}

// Reference coordinates
$reference = $_REQUEST['reference'];
if ($reference != "segment") {
  $argument .= " -reference $reference";
} 

//////////////////////////////////////////////////
// Write bed file in temporary directory
if (!$errors) {

  // Should be replaced by getTempFileName()
  $suffix = "_".date("Ymd_His")."_".randchar(3).".bed";
	  
  // Upload file from client machine


  if ($_FILES["bedfile"]['name'] != "") {
    $bed_file_name = basename($_FILES['bedfile']['name']);

    // We need to keep the original extension in order to support
    // uncompression of .gz files
    $extension = end(explode( ".", $bed_file_name));
    if ($extension == "bed") {
      $bed_file_name = str_replace(".bed", "", $bed_file_name);
    }
    
    // Move uploaded bed file in tmp
    $bed_file = getTempFileName($bed_file_name, ".".$extension);

    if (move_uploaded_file($_FILES['bedfile']['tmp_name'], $bed_file)) {
      $argument .= " -i $bed_file";
    } else {
      error('File upload failed');
      $errors = true;
    }
  }
  
  // Bed data provided in text area
  if ($bed != "") {
    $array_line = explode("\n",$bed);
    // $bed_file = $properties['rsat_tmp']."/"."userbed".$suffix;
    $bed_file = getTempFileName("userbed", ".bed");
    $file = fopen ($bed_file, "w");
    $no_bed_line = true;
    $warnings = "";
    
    // Check input lines
    foreach($array_line as $line) {
      if (preg_match("/^[\w\-\+\s,\.\#; \/]+$/",$line)) {
	$no_bed_line = false;
	fwrite($file, $line."\n");
      } else if ($line != "") {
	$warnings .= htmlspecialchars($line)."<br/>\n";
      }
    }
    fclose($file);
    
    if ($warnings != "") {
      warning("Wrongly formatted line for bed format:<br/>\n".$warnings);
    }
    
    if ($no_bed_line) {
      error("Problem with pasted coordinates: not a single line is in bed format.");
      $errors = true;
    } else {
      $argument .= " -i $bed_file";
    }
  }

  // Bed data provided as an URL
  if ($sequence_url != "") {
    $url = explode("/",$sequence_url);
    $bed_file_name = end($url);

    // We need to keep the original extension in order to support
    // uncompression of .gz files
    $extension = end(explode( ".", $bed_file_name));
    if ($extension == "bed") {
      $bed_file_name = str_replace(".bed", "", $bed_file_name);
    }

    if ($url[0]=="http:" or $url[0]=='ftp:') {
      $argument .= " -u $sequence_url";
      
      //Add randum value to $bedfile for the outputfile
      //      $bed_file = $properties['rsat_tmp']."/".$bed_file;
      $bed_file = getTempFileName($bed_file_name, ".".$extension);
      
    } else {
      error($sequence_url."<p>Invalid URL (should start with http:// or ftp://)");
      $errors = true;
    }
  }
}

///////////////////////////////////////////
// Run fetch-sequences
if (!$errors) { 

  // Specify output file
  $output_file = $bed_file.".fasta";
  $output_file = str_replace(".bed.fasta",".fasta",$output_file);
  $output_file = str_replace(".gz.fasta",".fasta",$output_file);

  // Specify log file
  //  $error_file = str_replace(".bed","_log.txt",$bed_file);
  $error_file = $bed_file."_log.txt";
  $error_file = str_replace(".bed_log.txt","_log.txt",$error_file);
  $error_file = str_replace(".gz_log.txt","_log.txt",$error_file);
  
  $argument .= " -o $output_file";	  	
  if ($sequence_url == "") {
      $URL['Genomic coordinates (bed)'] = rsat_path_to_url($bed_file);
    } else {
      $URL['Genomic coordinates (bed)'] = $sequence_url;
    }
  $URL['Fetched sequences (fasta)'] = rsat_path_to_url($output_file);
  $URL['Log file (txt)'] = rsat_path_to_url($error_file);

  // Add arguments to the command
  $cmd .= $argument;	

  /* TEMPORARILY INACTIVATE BECAUSE IT DOES NOT WORK 
  // Announce job starting
  $msg = "Starting job.";
  if ($output =="email")  {
    $msg .= " After job completion, email will be sent to ".$user_email;
  }

  // Printing starting job
  info($msg);
  echo "<hr>";
  flush(); 
  */

  ////////////////////////////////////////////////////////////////
  // Send email with notification of starting task
  if ($output =="email")  {
    // Parammeters for sending the mail
    $to = $user_email;
    $subject = "[RSAT] fetch-sequences start ".$now."_".$suffix;
    
    // Store the URL table in a variable
    $html_mail = 0; // Boolean variable indicating whether HTML format is supported in email
    $headers = ""; // Header (specifying Mime types)
    if ($html_mail) {
      $msg = "<table class='resultlink'>\n";
      $msg .= "<tr><th colspan='2'>Result file(s)</th></tr>\n";
      foreach ($URL as $key => $value) {
	$msg .= "<tr><td>".$key."</td><td><a href = '".$value."'>".$value."</a></td></tr>\n";
      }
      $msg .= "</table>\n";
      
      $headers .= 'Mime-Version: 1.0'."\r\n";
      $headers .= 'Content-type: text/html; charset=utf-8'."\r\n";
      $headers .= "\r\n";
    } else {
      $msg = "fetch-sequences starting\n";
      $msg .= "Result files :\n";
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
      $mail_sent = mail($to, $subject." ; Job submitted", $msg, $headers);
      if ($mail_sent) {
	info("Job started, submission mail sent to ".$to);
      } else {
	error("Notification mail could not be sent.\n\n.".$msg);
      }
    }
  }

  // Display the command
  $cmd_report = str_replace($properties['RSAT'], '$RSAT', $cmd);
  info("Command : ".$cmd_report);
  echo "<hr>";

  // Run the command
  exec("$cmd >/dev/null", $error);
       
  /*
       
  //display log file
  
  $info = "";
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
//  echo $cmd;

  ////////////////////////////////////////////////////////////////
  // Send email with notification of task completion 
  if ($output =="email")  {
    //    echo "Email outpout not available for now";
    // Parammeters for sending the mail
    $to = $user_email;
    $subject = "[RSAT] fetch-sequences result ".$now."_".$suffix;
    
    // Store the URL table in a variable
    $html_mail = 0; // Boolean variable indicating whether HTML format is supported in email
    $headers = ""; // Header (specifying Mime types)
    if ($html_mail) {
      $msg = "<table class='resultlink'>\n";
      $msg .= "<tr><th colspan='2'>Result file(s)</th></tr>\n";
      foreach ($URL as $key => $value) {
	$msg .= "<tr><td>".$key."</td><td><a href = '".$value."'>".$value."</a></td></tr>\n"; 
      }
      $msg .= "</table>\n";
      
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
      $mail_sent = mail($to, $subject." ; Job completed", $msg, $headers);
      if ($mail_sent) {
	info("Job done, completion mail sent to ".$to);
      } else {
	error("Notification mail could not be sent.\n\n.".$msg);
      }
    }  
  }  

  //Display pipe
  echo ('
  <table class = "nextstep">
    <tr><th colspan=4>next step</th></tr>
    <tr>
      <td align=center>
	    <form method="post" action="peak-motifs_form.cgi">
	    <input type="hidden" name="title" value="'.str_replace('.bed', '', end(explode('/',$bed_file))).'">
	    <input type="hidden" name="sequence_url1" value="'.rsat_path_to_url($output_file).'">
	    <input type="submit" value="peak-motifs">
	    </form>
      </td>
         
      <td align=center>
	    <form method="post" action="random-genome-fragments_form.cgi">
	    <input type="hidden" name="title" value="'.str_replace('.bed', '', end(explode('/',$bed_file))).'">
	    <input type="hidden" name="sequence_url1" value="'.rsat_path_to_url($output_file).'">
	    <input type="submit" value="random-genome-fragments">
            <br>(<font color="red">Note: select organism manually</font>)
	    </form>
      </td>
    </tr>

  </table>');  
}	

?>
 
  </body>
</html>
