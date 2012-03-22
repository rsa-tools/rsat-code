<?php
## Load RSAT configuration
require ('functions.php');
#print_r($properties);
UpdateLogFile("rsat","","");

#Import variable from formulaire with prefixe fs_
import_request_variables('P','fs_');

#print_r($_POST);
#print_r($_FILES);

#define variable
$cmd = $properties['RSAT'].'/perl-scripts/fetch-sequences';
$argument = "";
$confict_extension = False;

##Check arguments
if ($fs_genome == "none") {
  echo "You have forgetten to indicate genome", "<br/>\n";
} else {
	$argument .= " -genome $fs_genome";
}

if ($fs_bed == "" and !isset($_FILES["bedfile"]) and $_FILES["bedfile"]["name"]=="") {
  echo "You have forgetten to indicate bed file", "<br/>\n";
}

##Write/move bed file in tmp
if (isset($_FILES["bedfile"])) {
	$file_name = basename($_FILES['bedfile']['name']);
	$file = $properties['rsat_tmp']."/".$file_name;
	if(move_uploaded_file($_FILES['bedfile']['tmp_name'], $file)) {
		$argument .= " -i $file";
		$output_file = str_replace(".bed",".fasta",$file);
		$argument .= " -o $output_file";
	} else {
		echo 'Upload Fail';
	}
}

##Check optional arguments
if ($fs_header == "galaxy") {
	$argument .= " -header galaxy";
}
/*
if (($fs_downstr_ext !="" or $fs_upstr_ext !="") and $fs_extend !="") {
	$confict_extension = True;
	echo "You have indicate one side extention and both side extention";
}
*/
if ($fs_downstr_ext !="") {
	if (!is_numeric($fs_downstr_ext)) {
		echo "Downstr. extend $fs_downstr_ext is not a integrer", "<br/>\n";
	} else {
		if (!$confict_extension) {
			$argument .= " -downstr_ext $fs_downstr_ext";
		}
	}
}

if ($fs_upstr_ext !="") {
	if (!is_numeric($fs_upstr_ext)) {
		echo "Upstr. extend $fs_upstr_ext is not a integrer", "<br/>\n";
	}  else {
		if (!$confict_extension) {
			$argument .= " -upstr_ext $fs_upstr_ext";
		}
	}
}

/*
if ($fs_extend !="") {
	if (!is_numeric($fs_extend)) {
		echo "Both side extend $fs_extend is not a integrer", "<br/>\n";
	} else {
		if (!$confict_extension) {
			$argument .= " -extend $fs_extend";
		}
	}
}
*/

if ($fs_reference != "segment") {
	$argument .= " -reference $fs_reference";
}

####Faire fichier

echo $argument, "<br/>\n";
echo $cmd.$argument, "<br/>\n";
$cmd = $properties['RSAT'].'/perl-scripts/supported-organisms-ucsc';
exec($cmd.$argument, $yy);
print_r($yy);
echo "fini";

?>