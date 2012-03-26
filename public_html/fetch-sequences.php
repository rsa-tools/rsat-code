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
$argument = " -v 2";
$confict_extension = False;
$exp_bed_file = "/^[\w\-\+\s,\.\#; \/]+$/";
?>

<html>
<head>
<title>RSAT - fetch-sequences</title>
<link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
   </head>
   <body class="results">
   
   <?php
////////////////////////////////////////////////////////////////
//Print <h3>
echo "<H3><a href='".$properties['rsat_www']."'>RSAT</a> - fetch-sequences - results</H3>";

////////////////////////////////////////////////////////////////
// Check arguments
$errors = 0;

// Check that genome has been specified
if ($fs_genome == "none" or $fs_genome == "" ) {
  error( "You forgot to specify the genome.");
  $errors=1;
 } else {
  $argument .= " -genome $fs_genome";
 }

//Check if email is right
if($fs_output =="email") {
  if (!preg_match("#^[\wàáâãäåæçèéêëìíîïðñòóôõöøùúûüý._\-]+@[a-z]+.[a-z]{2,4}$#", $fs_user_email)) {
     error( "Email not valid");
     $errors=1;
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
  $errors++;
 } else if ($bed_specifications > 1) {
  error("The BED file can be specified only in 1 way (paste, upload or URL of external server)");
  $errors++;
 }


// Header format
if ($fs_header == "galaxy") {
  $argument .= " -header galaxy";
 }


// Downstream extension
if (($fs_downstr_ext !="") &&
    ($fs_downstr_ext != '0')) {
  if (!is_numeric($fs_downstr_ext)) {
    error ($fs_downstr_extr." Invalid value for downstream extension: should be an Integrer.");
    $errors++;
  } else {
    if (!$confict_extension) {
      $argument .= " -downstr_ext $fs_downstr_ext";
    }
  }
 }

// Upstream extension
if (($fs_upstr_ext !="") &&
    ($fs_upstr_ext != '0')) {
  if (!is_numeric($fs_upstr_ext)) {
    error($fs_upstr_ext." is not valid value for upstream extension : should be an Integrer.");
    $errors++;
  }  else {
    if (!$confict_extension) {
      $argument .= " -upstr_ext $fs_upstr_ext";
    }
  }
 }

// Reference coordinates
if ($fs_reference != "segment") {
  $argument .= " -reference $fs_reference";
 }
 

/////////////////////////////////
// Write bed file in temporary directory
if (!$errors) {
  
  $now = date("Ymd_His");
  $suffix = randchar(3);
  
  // Upload file from client machine
  if ($_FILES["bedfile"]['name'] != "") {
    $bed_file_name = basename($_FILES['bedfile']['name']);
    $extension = end(explode(".", $bed_file_name));
    
    // Move uploaded bed file in tmp
    if ($extension =="bed") {
      $bed_file = $properties['rsat_tmp']."/".$bed_file_name;
      $bed_file = str_replace(".bed","_".$now."_".$suffix.".".$extension,$bed_file);
      if(move_uploaded_file($_FILES['bedfile']['tmp_name'], $bed_file)) {
	$argument .= " -i $bed_file";
      } else {
	error('File upload failed');
	$error = 1;
      }
    }	else {
      error("Wrong file extension $extention");
      $error = 1;
    }
  }
  
	 
  // Bed data provided in text area
  if ($fs_bed != "") {
    $array_line = explode("\n",$fs_bed);
    $bed_file = $properties['rsat_tmp']."/"."userbed_".$now."_".$suffix.".bed";
    $file = fopen ($bed_file, "w");
    
    foreach($array_line as $line) {
      if (preg_match($exp_bed_file,$line)) {
	fwrite($file, $line);
      } else {
	warning ("Not bed line ".htmlspecialchars($line));
      }
    }
    fclose($file);
    $argument .= " -i $bed_file";
  }	  

  // Check URL
  if ($fs_sequence_url != "") {
    $url = explode("/",$fs_sequence_url);
    $bed_file = end($url);
    $bed = explode(".",$bed_file);
    $extension = end($bed);

    if (($url[0]=="http:" or $url[0]=='ftp:') and $extension =="bed") {
      $argument .= " -u $fs_sequence_url";
  		
      //Add randum value to $bedfile for the outputfile
      $bed_file = $properties['rsat_tmp']."/".$bed_file;
      $bed_file = str_replace(".bed","_".$now."_".$suffix.".bed",$bed_file);
  		
    } else {
      error($fs_sequence_url." is not a valid URL (should start with http: or ftp:.");
      $errors++;
    }
  }
 }

///////////////////////////////////////////
// Run fetch-sequences
if ($errors == 0) { 

  // Specify output file
  $output_file = str_replace(".bed",".fasta",$bed_file);
  $error_file = str_replace(".bed","_log.txt",$bed_file);
  $argument .= " -o $output_file";	  	
  $URL['Genomic coordinates (bed)'] = rsat_path_to_url($bed_file);
  $URL['Fetched sequences (fasta)'] = rsat_path_to_url($output_file);
  $URL['Log file (txt)'] = rsat_path_to_url($error_file);

  // Add arguments to the command
  $cmd .= $argument;	

  // Run the command
  exec($cmd, $error);
  // print_r($error);
  // echo ("Done");

  if ($fs_output =="display")  {

    // Display the command
    $cmd_report = str_replace($properties['RSAT'], '$RSAT', $cmd);
    info("Command : ".$cmd_report);
    echo "<hr>";

    // Display the result
    print_url_table($URL);
    
  } else {
    echo "Email outpout not available for now";
/*
    //parammetter of mail
    $to = $fs_user_email;
    $subject = "fetch-sequences results";

    //put result in variable
    $msg = "<table class='resultlink'>\n";
    $msg .= "<tr><th colspan='2'>Result file(s)</th></tr>\n";
    foreach ($URL as $key => $value) {
      $msg .= "<tr><td>".$key."</td><td><a href = '".$value."'>".$value."</a></td></tr>\n"; 
    }
    $msg .= "</table>\n";
    
    $headers = 'Mime-Version: 1.0'."\r\n";
    $headers .= 'Content-type: text/html; charset=utf-8'."\r\n";
    $headers .= "\r\n";

    //sending mail
    mail($to, $subject, $msg, $headers);*/
  }  
}	

?>

</body>
</html>
