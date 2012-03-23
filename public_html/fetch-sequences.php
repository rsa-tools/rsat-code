<html>
  <head>
   <title>RSAT - fetch-sequences</title>
   <link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
  </head>
  <body class="results">
    <H3><a href='http://rsat.ulb.ac.be/'>RSAT</a> - fetch-sequences - results</H3>
  
  
<?php

## Load RSAT configuration
require ('functions.php');
#print_r($properties);
UpdateLogFile("rsat","","");

## Import variables with prefix fs_ from form
import_request_variables('P','fs_');

#print_r($_POST);
#print_r($_FILES);

## Initialize variables
$cmd = $properties['RSAT'].'/perl-scripts/fetch-sequences';
$argument = " -v 1";
$confict_extension = False;
$exp_bed_file = "/^[\w\-\+\s,\.\#; \/]+$/";

################################################################
## Check arguments
$errors = 0;

## Check that genome has been specified
if ($fs_genome == "none" or $fs_genome == "" ) {
  error( "You forgot to specify the genome", "<br/>\n");
  $errors = 1;
} else {
  $argument .= " -genome $fs_genome";
}


## Check that the BED has been specified once and only once
$bed_specifications = 0;
if ($fs_bed != "") {
  $bed_specifications++;
}
if ($_FILES["bedfile"]['name'] != "") {
  $bed_specifications++;
}

if ( $fs_sequence_url != "") {
  $bed_specifications++;
}

## Report error if no BED has been specified
if ($bed_specifications == 0) {
  error("You forgot to specify the genomic coordinates (BED file)");
  $errors = 1;
} else if ($bed_specifications > 1) {
  error("The BED file can be specified only in 1 way (paste, upload or URL of external server)");
  $errors = 1;
}


##Check optional arguments
if ($fs_header == "galaxy") {
  $argument .= " -header galaxy";
 }

if ($fs_downstr_ext != '0') {
  if (!is_numeric($fs_downstr_ext)) {
    error ("Downstr. extend $fs_downstr_ext is not a integrer");
    $errors = 1;
  } else {
    if (!$confict_extension) {
      $argument .= " -downstr_ext $fs_downstr_ext";
    }
  }
}

if ($fs_upstr_ext != '0') {
  if (!is_numeric($fs_upstr_ext)) {
     error ("Upstr. extend $fs_upstr_ext is not a integrer");
     $errors = 1;
  }  else {
    if (!$confict_extension) {
      $argument .= " -upstr_ext $fs_upstr_ext";
    }
  }
}

if ($fs_reference != "segment") {
  $argument .= " -reference $fs_reference";
}
 
#################################
##Write bed file
if (!$errors) {
  
	$now = date("Ymd_His");
	$suffix = randchar(3);
	
	##Move upload bed file in tmp
	if ($_FILES["bedfile"]['name'] != "") {
	  $bed_file_name = basename($_FILES['bedfile']['name']);
	  $extension = end(explode(".", $bed_file_name));

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
 
  ##Write bed file in tmp
	if ($fs_bed != "") {
	 	$array_ligne = explode("\n",$fs_bed);
	 	$bed_file = $properties['rsat_tmp']."/"."userbed_".$now."_".$suffix.".bed";
	 	$file = fopen ($bed_file, "w");

	 	foreach($array_ligne as $ligne) {
	 	  if (preg_match($exp_bed_file,$ligne)) {
	 			fwrite($file, $ligne);
	 	  } else {
	      warning ("Not bed line", $ligne);
	 	  }
	  }
	  fclose($file);
	  $argument .= " -i $bed_file";
  }	  

  ##Check url
  if ($fs_sequence_url != "") {
  	$url = explode("/",$fs_sequence_url);
  	$bed_file = end($url);
  	$bed = explode(".",$bed_file);
  	$extension = end($bed);

  	if (($url[0]=="http:" or $url[0]=='ftp:') and $extension =="bed") {
  		$argument .= " -u $fs_sequence_url";
  		
  		##Add randum value to $bedfile for the outputfile
  		$bed_file = $properties['rsat_tmp']."/".$bed_file;
  		$bed_file = str_replace(".bed","_".$now."_".$suffix.".bed",$bed_file);
  		
  	} else {
  		error("$fs_sequence_url is not a URL of bed file.");
  		$error = 1;
  	}
  }
  
  ##Add otput argument
	if (!$errors) {  
		$output_file = str_replace(".bed",".fasta",$bed_file);
		$argument .= " -o $output_file";	  
	
		$URL['Genomic coordinates (bed)'] = rsat_path_to_url($bed_file);
		$URL['Fetched sequences (fasta)'] = rsat_path_to_url($output_file);
	}	
}

if (!$errors) {
	####Faire fichier
	$cmd .= $argument;	
	exec($cmd, $error);

	##DISPLAY command
  $cmd = str_replace($properties['RSAT'], '$RSAT', $cmd);
	info("Command : ".$cmd);
	echo "<hr>";
	## DISPLAY THE RESULT
	print_url_table($URL);  
}	


?>
  </body>
</html>