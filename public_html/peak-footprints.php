<?php
// Load RSAT configuration
require ('functions.php');
UpdateLogFile("rsat","","");

//print_r($properties);
//print_r($_POST);
//print_r($_FILES);


// Import variables with prefix pf_ from form
import_request_variables('P','pf_');

//Fill buffer for flush()
echo str_repeat(" ", 1024), "\n";

//Print header
$prog_name = "peak-footprints";
$result = true;
require ('RSAT_header.php');

// Initialize variables
$cmd = 'python '.$properties['RSAT'].'/contrib/peak-footprints/peak-footprints';
$argument = " --v 2";
$argument .= " --pipeline ".$properties['RSAT']."/contrib/peak-footprints/default_pipeline.xml";
$argument .= " --db_root_path ".$properties['RSAT'].'/data/motif_databases';
$errors = false;

////////////////////////////////////////////////////////////////
// Check that the BED has been specified once and only once
$bed_specifications = 0;

// Bed data pasted in text area
if ($pf_bed != "") {
	$bed_specifications++;
}

// Local bed file on client machine
if ($_FILES["bed_file"]['name'] != "") {
	$bed_specifications++;
}

// Local bed file on client machine
if ($pf_bed_url != "") {
	$bed_specifications++;
}

// Report error if no BED has been specified
if ($bed_specifications == 0) {
	error("You forgot to specify the genomic coordinates (BED file)");
	$errors = true;
} else if ($bed_specifications > 1) {
	error("The BED file can be specified only in 1 way (paste, upload or URL of external server)");
	$errors = true;
}

// Check that the custom_motifs has been specified 0 or once.
$custom_motifs_specifications = 0;

// Bed data pasted in text area
if ($pf_custom_motifs != "") {
	$custom_motifs_specifications++;
}
// Local bed file on client machine
if ($_FILES["custom_motifs_file"]['name'] != "") {
	$custom_motifs_specifications++;
}

// Report error if no BED has been specified
if ($custom_motifs_specifications > 1) {
	error("The Custom_Motifs file can be specified only in 1 way (paste, upload)");
	$errors = true;
}

//////////////////////////////////////////////////////////////
//Function for checking same type of argument
function numeric($name, $value) {
	if (!is_numeric($value)) {
		error("$name '$value' is not a number");
		return false;
	}	else {
		if (0>=$value or $value>=1) {
			error("$name must be between 0 and 1");
			return false;
		} else {
			return true;
		}
	}
}

function integrer($name, $value) {
	if (!ereg("^[0-9]+$", $value)) {
		error("$name '$value' is not a integrer");
		return false;
	}	else {
		if (1>$value) {
			error("$name must be >1");
			return false;
		} else {
			return true;
		}
	}
}

//Check argument
if ($pf_genome == "0")  {
	error("You forgot to choose a genome");
	$errors = true;
} else {
	$argument .= " --ref_species $pf_genome";
}

if ($_POST['species_ali_keep'] == "") {
	error("You forgot to specify aligned species");
	$errors = true;
} else {
	$al_species = ' --align_species "';
	foreach($pf_species_ali_keep as $key => $value) {
		if ($value!=end($pf_species_ali_keep)) {
			$al_species .= "$value ";
		} else {
			$al_species .= $value.'"';
		}
	}
	$argument .= $al_species;
}

if ($pf_r_motif == "") {
	error("You forgot to specify Reference Motif");
	$errors = true;
} else {
	if (preg_match("/^[\w\.]+$/",$pf_r_motif)) {
		$argument .= " --ref_motif $pf_r_motif";
	}	
}

if ($pf_nb_peaks != '') {
	if (integrer("Maximun numbers of peaks", $pf_nb_peaks)) {
		$argument .= " --peak_number $pf_nb_peaks";
	}	else {
		$errors = true;
	}
}

if (numeric("Column-wise conservation threshold", $pf_cons_thres)) {
	if ($pf_cons_thres != 0.7) {
		$argument .= " --conservation_threshold $pf_cons_thres";
	}
} else {
	$errors = true;
}

if (numeric("Block-wise conservation threshold", $pf_window_cons_thres)) {
	if ($pf_window_cons_thres != 0.7) {
		$argument .= " --window_conservation_threshold $pf_window_cons_thres";
	}
} else {
	$errors = true;
}

if (integrer("Initial sliding window size", $pf_window_size)) {
	if ($pf_window_size != 5) {
		$argument .= " --window_size $pf_window_size";
	}
}	else {
		$errors = true;
}

if (integrer("Maximum number of motif reported", $pf_motif_number)) {
	if ($pf_motif_number != 50) {
		$argument .= " --max_motif_number $pf_motif_number";
	}
} else {
	$errors = true;
}

if (integrer("Maximum number of motif reported by family ", $pf_motif_number_family)) {
	if ($pf_motif_number_family != 4) {
		$argument .= " --window_size $pf_motif_number_family";
	}
}	else {
	$errors = true;
}

if($pf_output =="email") {
	if (!preg_match("#^[^@]+@([a-z]+\.)+[a-z]{2,4}$#", $pf_user_email)) {
		error( "Email not valid");
		$errors=true;
	}
}

//////////////////////////////////////////////////////////////////////////////
//Writing bed/custom_motifs file in tmp
if (!$errors) {
	$suffix = "_".date("Ymd_His")."_".randchar(3);
	$user_name = "users_pipeline$suffix";
	$argument .= " --pipeline_name ".$user_name;
	
	// Bed data provided in text area
	if ($pf_bed != "") {
		$array_line = explode("\n",$pf_bed);
		$bed_file = $properties['rsat_tmp']."/"."userbed".$suffix.".bed";
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
			$argument .= " --input_peaks $bed_file";
		}
	}
	
	// Upload bed file from client machine
	if ($_FILES["bed_file"]['name'] != "") {
		$bed_file_name = basename($_FILES['bed_file']['name']);

		// Move uploaded bed file in tmp
		$bed_file = $properties['rsat_tmp']."/".$bed_file_name;
		$extension = end( explode( ".", $bed_file));
		
		if ($extension == "bed") {
			$bed_file = str_replace(".bed",$suffix.".bed",$bed_file);
		} else {
			$bed_file = $bed_file.$suffix.".bed";
		}

		if(move_uploaded_file($_FILES['bed_file']['tmp_name'], $bed_file)) {
			$argument .= " --input_peaks $bed_file";
		} else {
			error('File upload failed');
			$errors = true;
		}
	}

	// Upload bed file from url
	if ($pf_bed_url != "") {
		$url = explode("/",$pf_bed_url);
		$bed_file = end($url);

		if ($url[0]=="http:" or $url[0]=='ftp:') {
			$web_bed = file_get_contents($pf_bed_url);
  		
			//Add randum value to $bedfile for the outputfile
			$bed_file = $properties['rsat_tmp']."/".$bed_file;
			$extension = end( explode( ".", $bed_file));
		
			if ($extension == "bed") {
				$bed_file = str_replace(".bed",$suffix.".bed",$bed_file);
			} else {
				$bed_file = $bed_file.$suffix.".bed";
			}

		$file = fopen ($bed_file, "w");
		fwrite($file, $web_bed);
		fclose($file);
		$argument .= " --input_peaks $bed_file";

		} else {
			error($fs_sequence_url." is not a valid URL (should start with http: or ftp:.");
			$errors = true;
		}
	}
	




	//Custom_Motifs data provided in text area
	if ($pf_custom_motifs != "") {
		$array_line = explode("\n",$pf_custom_motifs);
		$custom_motifs_file = $properties['rsat_tmp']."/"."userbed".$suffix.".transfac";
		$file = fopen ($custom_motifs_file, "w");
		$no_custom_motifs_line = true;
	
		foreach($array_line as $line) {
			if (preg_match("/^[\w\-\+\s,\.\#; \/]+$/",$line)) {
				$no_custom_motifs_line = false;
				fwrite($file, $line."\n");
			} else {
				warning (htmlspecialchars($line)." is not a bed line.\n");
			}
		}
		fclose($file);
	
		if ($no_custom_motifs_line) {
			error("All your line are not custom_motifs format");
			$errors = true;
		} else {
			$argument .= " --custo_db_file_path $custom_motifs_file";  
		}
	}	
	
	// Upload custom_motifs file from client machine
	if ($_FILES["custom_motifs_file"]['name'] != "") {
		$custom_motifs_file_name = basename($_FILES['custom_motifs_file']['name']);
		$extension = end( explode( ".", $custom_motifs_file_name));
		
		// Move uploaded custom_motifs file in tmp
		$custom_motifs_file = $properties['rsat_tmp']."/".$custom_motifs_file_name;
		
		if ($extension == "tf") {
			$custom_motifs_file = str_replace(".tf",$suffix.".tf",$custom_motifs_file);
		} else {
			$custom_motifs_file .= $suffix.".tf";
		}
		
		
		
		
		if(move_uploaded_file($_FILES['custom_motifs_file']['tmp_name'], $custom_motifs_file)) {
			$argument .= " --custo_db_file_path $custom_motifs_file"; 
		} else {
			error('File upload failed');
			$errors = true;
		}
	}
}

/////////////////////////////////////////////////////////////////////////////
if (!$errors) {	

	$maf_path = $properties['RSAT']."/".$pf_multiz_path."/";
	$argument .= " --maf_path $maf_path"; 

	$outpout_path = $properties['rsat_tmp']."/peak-footprints_output/output/$user_name";
	
	//Prepare URL result table
	$URL['Genomic coordinates (bed)'] = rsat_path_to_url($bed_file);
	if ($custom_motifs_specifications == 1) {
		$URL['Custom motifs (transfac format)'] = rsat_path_to_url($custom_motifs_file);
	}
	$URL['Progress can bo follow here'] = rsat_path_to_url("$outpout_path/progression.xml");
	$URL['When progress is complete, result page appears here'] = rsat_path_to_url($outpout_path."/".$user_name."_9_FinalOutputProcessor/".$user_name."_MotifClassification.xml");
	
	// Add arguments to the command
	$cmd .= $argument;

	// Announce job starting
	$msg = "Starting job.";
	if ($fs_output =="email")  {
		$msg .= " After job completion, email will be sent to ".$fs_user_email;
	}
	info($msg);
	echo "<hr>";
	
	// Display the command
	$cmd_report = str_replace($properties['RSAT'], '$RSAT', $cmd);
	info("Command : ".$cmd_report);
	echo "<hr>";

	flush();
	
	// Display the result
	print_url_table($URL);


	///////////////////////////////////////////////////////////////
	// Send email with notification of starting task
	if ($fs_output =="email")  {
		// Parammeters for sending the mail
		$to = $fs_user_email;
		$subject = "[RSAT] $prog_name begin ".$suffix;
	
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
			$msg = "$prog_name result\n\n";
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
				info("Job start, email sent to ".$to);
			} else {
				error("Notification mail could not be sent.\n\n.".$msg);
			}
		}
	}
				
	// Run the command
	exec($cmd, $error);	

	///////////////////////////////////////////////////////////////
	// Send email with notification of task completion
	if ($fs_output =="email")  {
		$subject = "[RSAT] $prog_name result ".$suffix;

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
}

?>
	</body>
</html>
