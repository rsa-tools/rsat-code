<html>
<head>
<title>RSAT - fetch-sequences</title>
<link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
   </head>
   <body class="results"> 
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
if ($fs_genome == "none" or $fs_genome == "" ) {
  error( "You forgot to specify the genome.");
  $errors = true;
 } else {
  $argument .= " -genome $fs_genome";
 }

//Check syntax of email address (ensure the texte netered in email box was an email address)
if($fs_output =="email") {
  if (!preg_match("#^[^@\.]+(\.[^@]+)*@[a-zA-Z0-9-]+(\.[a-zA-Z0-9-]+)*(\.[a-zA-Z]{2,4})$#", $fs_user_email)) {
     error( "Email not valid");
     $errors=true;
  }
 }

////////////////////////////////////////////////////////////////
// Check that the BED has been specified once and only once
$bed_specifications = 0;

// Bed data pasted in text area
if ($fs_bed != "") {
  $bed_specifications++;
 }
// Local bed file on client machine
if ($_FILES["bedfile"]['name'] != "") {
  $bed_specifications++;
 }
// Bed specified by URL of remote server
if ( $fs_sequence_url != "") {
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
if ($fs_header == "galaxy") {
  $argument .= " -header_format galaxy";
 }

// Downstream extension
if (($fs_downstr_ext !="") &&
    ($fs_downstr_ext != '0')) {
  if (!is_numeric($fs_downstr_ext)) {
    error ($fs_downstr_extr." Invalid value for downstream extension: should be an Integrer.");
    $errors = true;
  } else {
     $argument .= " -downstr_ext $fs_downstr_ext";
  }
}

// Upstream extension
if (($fs_upstr_ext !="") &&
    ($fs_upstr_ext != '0')) {
  if (!is_numeric($fs_upstr_ext)) {
    error($fs_upstr_ext." is not valid value for upstream extension : should be an Integrer.");
    $errors = true;
  }  else {
		$argument .= " -upstr_ext $fs_upstr_ext";
  }
}

// Reference coordinates
if ($fs_reference != "segment") {
  $argument .= " -reference $fs_reference";
} 

//////////////////////////////////////////////////
// Write bed file in temporary directory
if (!$errors) {
	$suffix = "_".date("Ymd_His")."_".randchar(3).".bed";
	  
  // Upload file from client machine
  if ($_FILES["bedfile"]['name'] != "") {
    $bed_file_name = basename($_FILES['bedfile']['name']);
    
    // Move uploaded bed file in tmp
    $bed_file = $properties['rsat_tmp']."/".$bed_file_name;
    $extension = end( explode( ".", $bed_file));
    
    if ($extension == "bed") {
    	$bed_file = str_replace(".bed",$suffix,$bed_file);
    } else {
    	$bed_file = $bed_file.$suffix.".".$extension; ## We hold the original extension in order to support uncompression of .gz files
    }
    
    if(move_uploaded_file($_FILES['bedfile']['tmp_name'], $bed_file)) {
			$argument .= " -i $bed_file";
    } else {
			error('File upload failed');
			$errors = true;
    }
  }
  
  // Bed data provided in text area
  if ($fs_bed != "") {
  	$array_line = explode("\n",$fs_bed);
  	$bed_file = $properties['rsat_tmp']."/"."userbed".$suffix;
  	$file = fopen ($bed_file, "w");
  	$no_bed_line = true;
  	$warnings = "";
  
  	foreach($array_line as $line) {
  		if (preg_match("/^[\w\-\+\s,\.\#; \/]+$/",$line)) {
  			$no_bed_line = false;
  			fwrite($file, $line."\n");
  		} else {
  			$warnings .= htmlspecialchars($line)."<br/>\n";
  		}
  	}
  	fclose($file);
  
  	if ($warnings != "") {
  		warning("Not bed line :<br/>\n".$warnings);
  	}
  
  	if ($no_bed_line) {
  		error("All your line are not bed format");
  		$errors = true;
  	} else {
  		$argument .= " -i $bed_file";
  	}
  }

  // Check URL
  if ($fs_sequence_url != "") {
    $url = explode("/",$fs_sequence_url);
    $bed_file = end($url);

    if ($url[0]=="http:" or $url[0]=='ftp:') {
      $argument .= " -u $fs_sequence_url";
	  		
      //Add randum value to $bedfile for the outputfile
      $bed_file = $properties['rsat_tmp']."/".$bed_file;
      $extension = end( explode( ".", $bed_file));
    
      if ($extension == "bed") {
	$bed_file = str_replace(".bed",$suffix,$bed_file);
      } else {
	//	$bed_file = $bed_file.$suffix;
    	$bed_file = $bed_file.$suffix.".".$extension; ## We hold the original extension in order to support uncompression of .gz files
      }
	    
    } else {
      error($fs_sequence_url." is not a valid URL (should start with http: or ftp:.");
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
  $URL['Genomic coordinates (bed)'] = rsat_path_to_url($bed_file);
  $URL['Fetched sequences (fasta)'] = rsat_path_to_url($output_file);
  $URL['Log file (txt)'] = rsat_path_to_url($error_file);

  // Add arguments to the command
  $cmd .= $argument;	

  /* TEMPORARILY INACTIVATE BECAUSE IT DOES NOT WORK 
  // Announce job starting
  $msg = "Starting job.";
  if ($fs_output =="email")  {
    $msg .= " After job completion, email will be sent to ".$fs_user_email;
  }

  // Printing starting job
  info($msg);
  echo "<hr>";
  flush(); 
  */

  ################################################################
  // Send email with notification of starting task
  if ($fs_output =="email")  {
  	// Parammeters for sending the mail
  	$to = $fs_user_email;
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
  		$msg .= "Result files will :\n";
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
  			info("Job start, email sent to ".$to);
  		} else {
  			error("Notification mail could not be sent.\n\n.".$msg);
  		}
  	}
  }
  
  // Run the command
  exec("$cmd >/dev/null", $error);

  /*
  // Display the command
  $cmd_report = str_replace($properties['RSAT'], '$RSAT', $cmd);
  info("Command : ".$cmd_report);
  echo "<hr>";
  
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

  ################################################################
  // Send email with notification of task completion 
  if ($fs_output =="email")  {
    //    echo "Email outpout not available for now";
    // Parammeters for sending the mail
    $to = $fs_user_email;
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
      $mail_sent = mail($to, $subject, $msg, $headers);
      if ($mail_sent) {
	info("Job done, email sent to ".$to);
      } else {
	error("Notification mail could not be sent.\n\n.".$msg);
      }
    }  
  }  

  //Display pipe
  echo ('   <table class = "nextstep">
    <tr><th>next step</th></tr>
    <tr><td align=center>
	    <form method="post" action="peak-motifs_form.cgi">
	    <input type="hidden" name="title" value="'.str_replace('.bed', '', end(explode('/',$bed_file))).'">
	    <input type="hidden" name="sequence_url1" value="'.rsat_path_to_url($output_file).'">
	    <input type="submit" value="peak-motifs">
	    </form></td>
	  </tr>
	  </table>');  
}	

?>
 
  </body>
</html>
