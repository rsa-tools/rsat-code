<?php
// Load RSAT configuration
require ('functions.php');
UpdateLogFile("rsat","","");

//print_r($properties);
//print_r($_POST);
//print_r($_FILES);

// Import variables with prefix fs_ from form
import_request_variables('P','pf_');

//Fill buffer for flush()
echo str_repeat(" ", 1024), "\n";

//Print header
$prog_name = "peak-footprints";
$result = true;
require ('RSAT_header.php');

// Initialize variables
$cmd = $properties['RSAT'].'/python-scripts/peak-footprints';
$argument = " -v 2";
$errors = false;

////////////////////////////////////////////////////////////////
// Check that the BED has been specified once and only once
$bed_specifications = 0;

// Bed data pasted in text area
if ($pf_bed != "") {
	$bed_specifications++;
}
// Local bed file on client machine
if ($_FILES["bedfile"]['name'] != "") {
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

// Check that the BED has been specified 0 or once.
$transfac_specifications = 0;

// Bed data pasted in text area
if ($pf_transfac != "") {
	$transfac_specifications++;
}
// Local bed file on client machine
if ($_FILES["transfac_file"]['name'] != "") {
	$transfac_specifications++;
}

// Report error if no BED has been specified
if ($transfac_specifications > 1) {
	error("The Transfac file can be specified only in 1 way (paste, upload)");
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

function int($name, $value) {
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
	$error++;
}

if (int("Initial sliding window size", $pf_window_size)) {
	if ($pf_window_size != 5) {
		$argument .= " --window_size $pf_window_size";
	}
}	else {
		$errors = true;
}

if (int("Maximum number of motif reported", $pf_motif_number)) {
	if ($pf_motif_number != 50) {
		$argument .= " --max_motif_number $pf_motif_number";
	}
}	else {
	$errors = true;
}

if (int("Maximum number of motif reported by family ", $pf_motif_number_family)) {
	if ($pf_motif_number_family != 4) {
		$argument .= " --window_size $pf_motif_number_family";
	}
}	else {
	$errors = true;
}

//////////////////////////////////////////////////////////////////////////////
//Writing bed/transfac file in tmp
if (!$errors) {
	$suffix = "_".date("Ymd_His")."_".randchar(3);
	
	// Bed data provided in text area
	if ($pf_bed != "") {
		$array_line = explode("\n",$pf_bed);
		$bed_file = $properties['rsat_tmp']."/"."userbed".$suffix.".bed";
		$file = fopen ($bed_file, "w");
		$no_bed_line = true;
		
		foreach($array_line as $line) {
			if (preg_match("/^[\w\-\+\s,\.\#; \/]+$/",$line)) {
				$no_bed_line = false;
				fwrite($file, $line);
			} else {
						warning (htmlspecialchars($line)." is not a bed line.\n");
			}
		}
		fclose($file);
		
		if ($no_bed_line) {
			error("All your line are not bed format");
			$errors = true;
		} else {
			$argument .= " --input_peaks $bed_file";
		}
	}
	
	// Upload bed file from client machine
	if ($_FILES["bedfile"]['name'] != "") {
		$bed_file_name = basename($_FILES['bedfile']['name']);

		// Move uploaded bed file in tmp
		if (end(explode(".", $bed_file_name)) =="bed") {
			$bed_file = $properties['rsat_tmp']."/".$bed_file_name;
			$bed_file = str_replace(".bed",$suffix.".bed",$bed_file);
			if(move_uploaded_file($_FILES['bedfile']['tmp_name'], $bed_file)) {
				$argument .= " --input_peaks $bed_file";
			} else {
				error('File upload failed');
				$errors = true;
			}
		}	else {
			error("Wrong file extension");
			$errors = true;
		}
	}

	//Transfac data provided in text area
	if ($pf_transfac != "") {
		$array_line = explode("\n",$pf_transfac);
		$transfac_file = $properties['rsat_tmp']."/"."userbed".$suffix.".transfac";
		$file = fopen ($transfac_file, "w");
		$no_transfac_line = true;
	
		foreach($array_line as $line) {
			if (preg_match("/^[\w\-\+\s,\.\#; \/]+$/",$line)) {
				$no_transfac_line = false;
				fwrite($file, $line);
			} else {
				warning (htmlspecialchars($line)." is not a bed line.\n");
			}
		}
		fclose($file);
	
		if ($no_transfac_line) {
			error("All your line are not transfac format");
			$errors = true;
		} else {
			$argument .= " --??? $bed_file";  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		}
	}	
	
	// Upload transfac file from client machine
	if ($_FILES["transfacfile"]['name'] != "") {
		$transfac_file_name = basename($_FILES['transfacfile']['name']);
	
		// Move uploaded transfac file in tmp
		if (end(explode(".", $transfac_file_name)) == "transfac") {
			$transfac_file = $properties['rsat_tmp']."/".$transfac_file_name;
			$transfac_file = str_replace(".transfac",$suffix.".transfac",$transfac_file);
			if(move_uploaded_file($_FILES['transfacfile']['tmp_name'], $transfac_file)) {
				$argument .= " --??? $transfac_file";   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			} else {
				error('File upload failed');
				$errors = true;
			}
		}	else {
			error("Wrong file extension"); 
			$errors = true;
		}
	}
}

/////////////////////////////////////////////////////////////////////////////
if (!$errors) {
	
	// Specify output file
	$output_file = str_replace(".bed",".fasta",$bed_file);
	$argument .= " -o $output_file";
	$URL['Genomic coordinates (bed)'] = rsat_path_to_url($bed_file);
	$URL['Motif (transfac)'] = rsat_path_to_url($transfac_file);
	$URL['Peak-footprints (fasta)'] = rsat_path_to_url($output_file);
	
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
/* $cmd_report = str_replace($properties['RSAT'], '$RSAT', $cmd);
	info("Command : ".$cmd_report);
	echo "<hr>";
*/
	flush();
	
	// Display the result
	print_url_table($URL);
	
	///////////////////////////////////////////////////////////////
	// Send email with notification of task completion
	if ($fs_output =="email")  {
		// Parammeters for sending the mail
		$to = $fs_user_email;
		$subject = "[RSAT] $prog_name result ".$now."_".$suffix;
	
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