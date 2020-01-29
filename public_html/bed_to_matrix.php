<?php
    session_start();
    $_SESSION['form_uri_bed_to_matrix'] = $_SESSION['form_uri_bed_to_matrix'];
?>
<!DOCTYPE html>
<!-- <head>
        <title>RSAT - <?php print basename(__FILE__); ?></title>
</head>
<body>
    <h3 align='center'>           
        <?php
        print "<a href='" . $_SERVER['REQUEST_URI'] . "/rsat/'>RSAT</a> ("
                . basename(__FILE__) . ")";
        ?>
    </h3> -->
    
 <html>
<head>
<title>RSAT - Features matrix</title>
<link rel="stylesheet" type="text/css" href = "css/main_grat.css" media="screen">
   </head>
   <body class="results">     
    
 
    <?php
   
    require('functions.php');
    $RSAT_config_props=load_props(__DIR__."/RSAT_config.props");
    $RSAT = $RSAT_config_props["RSAT"];
    $rand_dir = randchar(6);
    $rsat_tmp = $RSAT_config_props["rsat_tmp"];
    
    $workingdir = $rsat_tmp."/SVM/svmrsat" . $rand_dir;
    $_SESSION['workingdir_bed_to_matrix'] = $workingdir;
    $_SESSION['previous_bed_to_matrix'] = basename(__FILE__);
    $_SESSION['expected_workingdir_file_nb_bed_to_matrix'] = 11;
    shell_exec("/bin/mkdir -m 777 $workingdir");

    echo "<H3><a href='".$properties['rsat_www']."'>RSAT</a> - Features matrix - results</H3>";


////////////////////////////////////////////////////////////////////////////////////////////////////////////

$stat = $_REQUEST['stat'];
$background_seqs = $_REQUEST['background_seqs'];
$genome = $_REQUEST['genome'];

if ($genome == "none" or $genome == "" ) {
  error( "You forgot to specify the genome.");
  $errors = true;
 } 
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////MATRICES////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Check that the custom_motifs has been specified 0 or once.
$custom_motifs_specifications = 0;

// Matrix data pasted in text area
if ($_REQUEST['matrix'] != "") {
	$custom_motifs_specifications++;
}
// Local matrix file on client machine
if ($_FILES["matrix_file"]['name'] != "") {
	$custom_motifs_specifications++;
}
// Matrix specified by URL of remote server
$matrix_file_url = $_REQUEST['matrix_file_url'];
if ( $matrix_file_url != "") {
  $custom_motifs_specifications++;
 }

// Report error if no Matrix has been specified
if ($custom_motifs_specifications > 1) {
	error("The Custom_Motifs file can be specified only in 1 way (paste, upload)");
	$errors = true;
}
//$matrix_format = $_REQUEST['matrix_format'];

//Custom_Motifs data provided in text area
  if ($_REQUEST['matrix'] != "") {
    $array_line = explode("\n",$_REQUEST['matrix']);
    $tmp_tf_matrix = $workingdir . "/matrix.tf";
    $file = fopen ($tmp_tf_matrix, "w");
    $no_custom_motifs_line = true;
	
    foreach($array_line as $line) {
      if (preg_match("/^[\w\-\+\s,\.\#; \/]+$/",$line)) {
	$no_custom_motifs_line = false;
	fwrite($file, $line."\n");
      } 
    }
    fclose($file);
    if ($no_custom_motifs_line) {
      error("All your line are not custom_motifs format");
      $errors = true;
    } 
  }	
	
  // Upload custom_motifs file from client machine
  if ($_FILES["matrix_file"]['name'] != "") {
    $custom_motifs_file_name = basename($_FILES['matrix_file']['name']);		
    // Move uploaded custom_motifs file in tmp
    $tmp_tf_matrix = $workingdir."/".$custom_motifs_file_name;
    if(move_uploaded_file($_FILES['matrix_file']['tmp_name'], $tmp_tf_matrix)) {
    } else {
      error('File upload failed');
      $errors = true;
    }
  }
  
  // Upload positive matrix file from url
  if ($matrix_file_url != "") {
    $matrix_url = explode("/",$matrix_file_url);
    $matrixfn = end($matrix_url);
    
    if ($matrix_url[0]=="http:" or $matrix_url[0]=='ftp:') {
      $web_matrix = file_get_contents($matrix_file_url);
      $tmp_tf_matrix = $workingdir . "/" . $matrixfn;
      $file = fopen ($tmp_tf_matrix, "w");
      fwrite($file, $web_matrix);
      fclose($file);

    } else {
      error(htmlentities($matrix_file_url)." is not a valid URL (should start with http: or ftp:.");
      $errors = true;
    }
  }
    
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////POSITIVE BED FILE///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Check that the positive BED has been specified once and only once
$bed_specifications = 0;

// Bed data pasted in text area
$textarea_bed = $_REQUEST['bed'];
if ($textarea_bed != "") {
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

 
// Bed data pasted in text area
 if ($textarea_bed != "") {
    $array_line = explode("\r",$textarea_bed);
    // $bed_file = $properties['rsat_tmp']."/"."userbed".$suffix.".bed";
    $tmp_bed = $workingdir . "/Positive_userbed.bed";

    $file = fopen ($tmp_bed, "w");
    $no_bed_line = true;
    $warnings = "";

    foreach($array_line as $line) {
      if (preg_match("/^[\w\-\+\s,\.\#; \/]+$/",$line)) {
	$no_bed_line = false;
	fwrite($file, $line);
      } else {

	//$warnings .= htmlspecialchars($line)."<br/>\n";
      }
    }
    fclose($file);

    if ($warnings != "") {
      print ("Not bed line :<br/>\n".$warnings);
    }
    
    if ($no_bed_line) {
      print ("All your line are not bed format in positive bed input.");
      $errors = true;
      
    } 
  }

  
  // Upload positive bed file from url
  if ($sequence_url != "") {
    $url = explode("/",$sequence_url);
    $bedfn = end($url);
   
    if ($url[0]=="http:" or $url[0]=='ftp:') {
      $web_bed = file_get_contents($sequence_url);
      $tmp_bed = $workingdir . "/Positive_userbed.bed";

      $file = fopen ($tmp_bed, "w");
      fwrite($file, $web_bed);
      fclose($file);

    } else {
      error(htmlentities($sequence_url)." is not a valid URL (should start with http: or ftp:).");
      $errors = true;
    }
  }
  
  // Upload positive bed file from client machine
  if ($_FILES["bedfile"]['name'] != "") {
    $bedfn = basename($_FILES['bedfile']['name']);
    
    // Upload bed file to the directory for this query
    // $bed_file = $properties['rsat_tmp']."/".$bed_file_name;
    $tmp_bed = $workingdir . "/Positive_userbed.bed";
    $extension = end( explode( ".", $tmp_bed));
    
    if ($extension == "bed") {
      $tmp_bed = str_replace(".bed",$suffix.".bed",$tmp_bed);
    } else {
      $tmp_bed = $tmp_bed.$suffix.".bed";
    }
    
    if(move_uploaded_file($_FILES['bedfile']['tmp_name'], $tmp_bed)) {
      //$argument .= " --input_peaks $bed_file";
    } else {
      error('Positive bed file upload failed');
      $errors = true;
    }
  }
 
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////NEGATIVE BED FILE///////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////// 
 
 // Check that the positive BED has been specified once and only once

$Negbed_specifications = 0;

// Bed data pasted in text area
$textarea_Negbed = $_REQUEST['Negbed'];
if ($textarea_Negbed != "") {
  $Negbed_specifications++;
 }
// Local bed file on client machine
if ($_FILES["Negbedfile"]['name'] != "") {
  $Negbed_specifications++;
 }

// Bed specified by URL of remote server
$Negequence_url = $_REQUEST['Negsequence_url'];
if ( $Negequence_url != "") {
  $Negbed_specifications++;
 }
 
$shuffling = $_REQUEST['shuffling'];
 if ( $shuffling != "") {
  $Negbed_specifications++;
 }
// error if no BED has been specified
if ($Negbed_specifications == 0) {
	$shuffling = $_REQUEST['shuffling'];
	if($_REQUEST['shuffling'] = ''){
		$shuffling = 2;
  		warning("No negative bed was specified ".$shuffling." times random sequence set will be generated");
  		}
	
 } else if ($Negbed_specifications > 1) {
  error("The negative BED file can be specified only in 1 way (paste, upload or URL of external server)");
  $errors = true;
 }

// Negative Bed data provided in text area
	
  if ($textarea_Negbed != "") {
    $array_line = explode("\n",$textarea_Negbed);
    $Negbed_file = $workingdir . "/Negative_userbed.bed";
    $Negfile = fopen ($Negbed_file, "w");
    $no_bed_line = true;
    $warnings = "";
    
    foreach($array_line as $line) {
      if (preg_match("/^[\w\-\+\s,\.\#; \/]+$/",$line)) {
	$no_bed_line = false;
	fwrite($Negfile, $line."\n");
      } else {
	$warnings .= htmlspecialchars($line)."<br/>\n";
      }
    }
    fclose($Negfile);
    
    if ($warnings != "") {
      warning("Not bed line :<br/>\n".$warnings);
    }
    
    if ($no_bed_line) {
      error("All your line are not bed format in negative bed file.");
      $errors = true;
    } 
  }
  
  // Upload Negative bed file from url
  if ($Negequence_url != "") {
    $Negurl = explode("/",$Negequence_url);
    $Negbedfn = end($Negurl);
    
    if ($Negurl[0]=="http:" or $Negurl[0]=='ftp:') {
      $Negweb_bed = file_get_contents($Negequence_url);
      
      // $bed_file = $properties['rsat_tmp']."/".$bed_file;
    	$Negbed_file = $workingdir . "/Negative_userbed.bed";
      $Negfile = fopen ($Negbed_file, "w");
      fwrite($Negfile, $Negweb_bed);
      fclose($Negfile);

    } else {
      error(htmlentities($Negequence_url)." is not a valid URL (should start with http: or ftp:).");
      $errors = true;
    }
  }
    
  // Upload positive bed file from client machine
  if ($_FILES["Negbedfile"]['name'] != "") {
    $Negbedfn = basename($_FILES['Negbedfile']['name']);
    
    // Upload bed file to the directory for this query
    // $bed_file = $properties['rsat_tmp']."/".$bed_file_name;
    $Negbed_file = $workingdir . "/Negative_userbed.bed";
    $extension = end( explode( ".", $Negbed_file));
    
    if ($extension == "bed") {
      $Negbed_file = str_replace(".bed",$suffix.".bed",$Negbed_file);
    } else {
      $Negbed_file = $Negbed_file.$suffix.".bed";
    }
    
    if(move_uploaded_file($_FILES['Negbedfile']['tmp_name'], $Negbed_file)) {
      //$argument .= " --input_peaks $bed_file";
    } else {
      error('Positive bed file upload failed');
      $errors = true;
    }
  }
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////WIG FILE(S) UPLOAD//////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Upload wig file from url
  if ($wig_url != "") {
    $wigurl = explode("/",$wig_url);
    $wigfn = end($wigurl);
    
    if ($wigurl[0]=="http:" or $wigurl[0]=='ftp:') {
      $wigweb_wig = file_get_contents($wig_url);
      
      // $bed_file = $properties['rsat_tmp']."/".$bed_file;
    	$wig_file = $workingdir . "/" . $wigfn;
      $wigfile = fopen ($wig_file, "w");
      fwrite($wigfile, $wigweb_wig);
      fclose($wigfile);

    } else {
      error(htmlentities($wig_url)." is not a valid URL (should start with http: or ftp:).");
      $errors = true;
    }
  }

// upload multiple wig files
$fileUploaded = false;

//You will have to check for each member in array if it was a successful file upload
foreach($_FILES["wigs"]["error"] as $key => $value){
    if($value == 0){
        $fileUploaded = true;
        break;
    }
}

//if(isset($_FILES['wigs'])){
if($fileUploaded){
	foreach($_FILES['wigs']['tmp_name'] as $key => $tmp_name)
	{
   	$file_name = $_FILES['wigs']['name'][$key];
   	$file_tmp =$_FILES['wigs']['tmp_name'][$key];
	   $wig_file = $workingdir . "/" . $file_name ;
    	if(move_uploaded_file($_FILES['wigs']['tmp_name'][$key], $wig_file)) {
      	$tmp_wig .= $wig_file .",";
    	} 
    	else {
      	error('wig file upload failed');
      	$errors = true;
    	}
  }
  $tmp_wig = rtrim($tmp_wig, ",");
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////background freq file BED FILE///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

$BGfreq_specifications = 0;

if ($_FILES["backgroundFreqs"]['name'] != "") {
  $BGfreq_specifications++;
 }

$backgroundFreqs_url = $_REQUEST['backgroundFreqs_url'];
echo $backgroundFreqs_url;
if ( $backgroundFreqs_url != "") {
  $BGfreq_specifications++;
 }

if ($BGfreq_specifications == 0) {
  error("You forgot to specify the background frequencies.");
  $errors=true;
 } else if ($BGfreq_specifications > 1) {
  error("The background frequencies file can be specified only in 1 way (paste, upload or URL of external server)");
  $errors = true;
 }
 
if ($backgroundFreqs_url != "") {
    $BGFreqsurl = explode("/",$backgroundFreqs_url);
    $BGFreqsfn = end($BGFreqsurl);
    $BGFreqweb = file_get_contents($backgroundFreqs_url);
    $tmp_bgFreq = $workingdir . "/background_frequencies.freq";
    $BGFreqfile = fopen ($tmp_bgFreq, "w");
    fwrite($BGFreqfile, $BGFreqweb);
    fclose($BGFreqfile);

    } else {
      $errors = true;
    }
  
if ($_FILES["backgroundFreqs"]['name'] != "") {
 
  $bgFreq = basename($_FILES['backgroundFreqs']['name']);
  $tmp_bgFreq = $workingdir . "/background_frequencies.freq";
  $extension = end( explode( ".", $tmp_bgFreq));
    
    if ($extension == "freq") {
      $tmp_bgFreq = str_replace(".freq",$suffix.".freq",$tmp_bgFreq);
    } else {
      $tmp_bgFreq = $tmp_bgFreq.$suffix.".freq";
    }
    
    if(move_uploaded_file($_FILES['backgroundFreqs']['tmp_name'], $tmp_bgFreq)) {
      //$argument .= " --input_peaks $bed_file";
    } else {
      error('background frequencies file upload failed');
      $errors = true;
    }
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////background seq file BED FILE///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Check that the positive BED has been specified once and only once
$BGseq_specifications = 0;

// Local bed file on client machine
if ($_FILES["backgroundFreqs"]['name'] != "") {
  $BGseq_specifications++;
 }

// Bed specified by URL of remote server
$backgroundSeqs_url = $_REQUEST['background_seqs_url'];

if ( $backgroundFreqs_url != "") {
  $BGseq_specifications++;
 }
 // Report error if no BED has been specified
if ($BGseq_specifications == 0) {
  $tmp_bgSeq = "";
 } else if ($BGseq_specifications > 1) {
  error("The background sequence file can be specified only in 1 way (paste, upload or URL of external server)");
  $errors = true;
 }
 
if ($backgroundSeqs_url != "") {
    $BGSequrl = explode("/",$backgroundSeqs_url);
    $BGSeqfn = end($BGSequrl);
    $BGSeqweb = file_get_contents($backgroundSeqs_url);
    $tmp_bgSeq = $workingdir . "/background_sequences.bed";
    $BGSeqfile = fopen ($tmp_bgSeq, "w");
    fwrite($BGSeqfile, $BGSeqweb);
    fclose($BGSeqfile);
}
if ($_FILES["background_seqs"]['name'] != "") {
  
  $bgSeq = basename($_FILES['background_seqs']['name']);
  $tmp_bgSeq = $workingdir . "/background_sequences.bed";
  $extension = end( explode( ".", $tmp_bgSeq));
    
    if ($extension == "bed") {
      $tmp_bgSeq = str_replace(".bed",$suffix.".bed",$tmp_bgSeq);
    } else {
      $tmp_bgSeq = $tmp_bgSeq.$suffix.".bed";
    }
    
    if(move_uploaded_file($_FILES['background_seqs']['tmp_name'], $tmp_bgSeq)) {
      //$argument .= " --input_peaks $bed_file";
    } else {
      error('background sequences file upload failed');
      $errors = true;
    }
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////COMMAND BUILDING////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
///chrom info


    
    
    $cmd = "export RSAT=".$RSAT;
    $cmd .= " ; /bin/bash " . $RSAT . "/R-scripts/R-scripts_SVM/bed_to_matrix.sh ";
    $cmd .= " --positive_bed " . $tmp_bed;
    if (isset($Negbed_file)){
    	$cmd .= " --negative_bed " . $Negbed_file;
    	}
    $cmd .= " -m " .  $tmp_tf_matrix;
    if ($fileUploaded){
    	$cmd .= " --ngs ".$tmp_wig;
    	}
    if ($stat){
    	$cmd .= " --stat ".$stat;
    	}
    $cmd .= " --background_freqs ".$tmp_bgFreq;
    $cmd .= " --genome " . $genome;
    if ($shuffling != 0){
    	$cmd .= " --shuffling " . $shuffling;
	if ($tmp_bgSeq != ''){
    		$cmd .= " --background_seqs ". $tmp_bgSeq;
		}
	}

    $cmd .= " --outdir " . $workingdir;
    $cmd .= " > /dev/null &";

   print "<p>Command: " . $cmd . "</p><br><br>";
   print "Thanks for your submission. Your task has been submitted to the RSAT server.<br><br>";
   print "Results will be available <a href='bed_to_matrix_results.php' >here</a>";

    exec($cmd);
    ?>
</body>
</html>
