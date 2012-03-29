<?php
// Load RSAT configuration
require ('functions.php');
UpdateLogFile("rsat","","");

//print_r($properties);
print_r($_POST);
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
}

if ($_POST['species_ali_keep'] == "") {
	error("You forgot to specify aligned species");
	$errors = true;
}

if ($pf_r_motif == "") {
	error("You forgot to specify Reference Motif");
	$errors = true;
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

echo $argument;

//////////////////////////////////////////////////////////////////////////////
//Writing bed/transfac file in tmp
if (!$errors) {
	$suffix = "_".date("Ymd_His")."_".randchar(3);
	
	// Bed data provided in text area
	if ($pf_bed != "") {
	}
}



/////////////////////////////////////////////////////////////////////////////
if (!$errors) {
	// Announce job starting
	$msg = "Starting job.";
	if ($fs_output =="email")  {
		$msg .= " After job completion, email will be sent to ".$fs_user_email;
	}
	
	info($msg);
	echo "<hr>";
	flush();
}

?>
	</body>
</html>