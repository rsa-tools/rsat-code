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
$argument = "";
$confict_extension = False;
$exp_bed_file = "/^[\w-\+\r,\n\t\.\#; \/]+$/";

################################################################
## Check arguments

$errors = 0;

## Check that genome has been specified
if ($fs_genome == "none") {
  error( "You forgot to specify the genome", "<br/>\n");
  $errors++;
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
  error("You forgot to specify the genomic coordinates (BED file)", "<br/>\n");
  $errors++;
 } else if ($bed_specifications > 1) {
  error("The BED file can be specified only in 1 way (paste, upload or URL of external server)", "<br/>");
  $errors++;
 }


## Write/move bed file in tmp
if ($_FILES["bedfile"]['name'] != "") {
  $file_name = basename($_FILES['bedfile']['name']);
  $file = $properties['rsat_tmp']."/".$file_name;
  if(move_uploaded_file($_FILES['bedfile']['tmp_name'], $file)) {
    $argument .= " -i $file";
    $output_file = str_replace(".bed",".fasta",$file);
    $argument .= " -o $output_file";
  } else {
    echo 'File upload failed';
  }
 }

if ($fs_bed != "") {
  $array_ligne = explode("\n",$fs_bed);
#print_r($array_ligne);
  foreach($array_ligne as $ligne) {
    if (preg_match($exp_bed_file,$ligne)) {
#echo "ok ", $ligne, "<br/>";
    } else 
      echo "Not bed line", $ligne, "<br/>\n";
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
$cmd .= $argument;
echo $cmd, "<br/>\n";
exec($cmd, $error);
print_r($error);
echo "fini";

?>