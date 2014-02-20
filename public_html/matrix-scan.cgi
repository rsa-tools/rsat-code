#!/usr/bin/perl
############################################################
#
# $Id: matrix-scan.cgi,v 1.46 2013/11/03 19:33:31 jvanheld Exp $
#
# Time-stamp: <2003-06-16 00:59:07 jvanheld>
#
############################################################
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}

use CGI;
use CGI::Carp qw/fatalsToBrowser/;
#### redirect error log to a file
BEGIN {
    $ERR_LOG = "/dev/null";
#    $ERR_LOG = &RSAT::util::get_pub_temp()."/RSA_ERROR_LOG.txt";
    use CGI::Carp qw(carpout);
    open (LOG, ">> $ERR_LOG")
	|| die "Unable to redirect log\n";
    carpout(*LOG);
}
require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

#$ENV{rsat_echo}=1;

### print the result page
&RSA_header("matrix-scan result", "results");
&ListParameters() if ($ENV{rsat_echo} >=2);

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

$command = $SCRIPTS."/matrix-scan -v 1 ";
$prefix = "matrix-scan";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); $tmp_file_name = &ShortFileName($tmp_file_path);
@result_files = ();

################################################################
## quick  mode
my $quick_mode = 0;
if ($query->param("quick")) {
	$parameters .= " -quick";
	$quick_mode = 1;
}

################################################################
## sequence file
if ($quick_mode) {
  ($sequence_file, $sequence_format) = &MultiGetSequenceFile(1, $tmp_file_path.".fasta", 1);
} else {
  ($sequence_file,$sequence_format) = &GetSequenceFile();
}


#### matrix-scan parameters
&ReadMatrixScanParameters();
$parameters .= " -i $sequence_file -seq_format $sequence_format";
push @result_files, ("Input sequence",$sequence_file);

################################################################
### vertically print the matrix
if ($query->param('vertically_print')) {
  $parameters .= " -p";
}


################################################################
## Treatment of N characters

if ($query->param('n_score')) {
 my $n_treatment = $query->param('n_score');
  $parameters .= " -n ".$n_treatment;
} else {
	$parameters .= " -n score ";
}

$command .= " $parameters";

## Define output file
$result_file =  $tmp_file_path.".ft";
push @result_files, ("Scan result (FT)",$result_file);

#$command .= " -o ".$result_file;

################################################################
#### echo the command (for debugging)
&ReportWebCommand($command);

################################################################
### execute the command ###
if ($query->param('output') eq "display") {

  unless ($query->param('table')) {
    &PipingWarning() unless ($query->param("analysis_type") eq "analysis_occ");
  }



  ### Print the result on Web page
  open RESULT, "$command |";
  print "<PRE>";
  &PrintHtmlTable(RESULT, $result_file, true, 1000);
  print "</PRE>";
  print "<HR SIZE = 3>";

#  system("$command");


  ################################################################
  ## Convert features to GFF3 format
  $gff3_file =  $tmp_file_path.".gff3";
  $command = "${SCRIPTS}/convert-features -from ft -to gff3 ";
  $command .= " -i ".$result_file;
  $command .= " -o ".$gff3_file;
  &RSAT::message::Info("Converting features to GFF3 format", $command) if ($ENV{rsat_echo} >=2);
  system($command);
  push @result_files, ("Features (GFF3)", $gff3_file);

  ################################################################
  ## Convert features for loading as genome tracks
  if ($origin eq "genomic") {
    $bed_file =  $tmp_file_path.".bed";
    $command = "${SCRIPTS}/convert-features -from ft -to bed ";
    $command .= " -i ".$result_file;
    $command .= " -o ".$bed_file;
    &RSAT::message::Info("Converting features to BED format", $command) if ($ENV{rsat_echo} >=2);;
    system($command);
    push @result_files, ("Features (BED)", $bed_file);
  }


  ################################################################
  ## Table with links to the raw result files in various formats
  &PrintURLTable(@result_files);

  if ($origin eq "genomic") {

    ################################################################
    ## Build the URL to view the features on a Genome Browser
    my $genomic_format = `grep 'Genomic coordinate format' $result_file`;
    chomp($genomic_format);

    my $browser = `grep -i 'browser url' $result_file`;
    chomp($browser);
    $browser =~ s/.*browser url\s+(\S+)/$1/i;

    if ($genomic_format =~ /Genomic coordinate format\s+(\S+)/) {
      $genomic_format = lc($1);
      &RSAT::message::Warning("Genomic format: $genomic_format") if ($ENV{rsat_echo} >=2);;
    }
    if ($genomic_format eq "ucsc") {
#      $browser_url = "<a target='_blank' href='http://genome.ucsc.edu/cgi-bin/hgCustom'>UCSC genome browser</a>",
      $browser_url = "<a target='_blank' href='";
      $browser_url .= $browser;
      $browser_url .= "&hgt.customText=".$result_URL.".bed";
      $browser_url .= "'><img border=0 height='20' src='images/UCSC_icon.jpg' alt='UCSC'></a>";
    } elsif ($genomic_format =~ /ensembl/i) {
#      &RSAT::message::Warning($browser);
      $browser_url = "<a target='_blank' href='";
      $browser_url .= $browser;
      $browser_url .= ";contigviewbottom=url:".$result_URL.".bed";
      $browser_url .= "=normal'><img border=0 height='20' src='images/e-ensembl_icon.png' alt='EnsEMBL'> genome browser</a>";
    }

    print ("<h3>Upload features in genome viewer ($genomic_format)</h3>",
	   $browser_url);
  }
  print "</table></center>";

  ## Collect genes for piping the results to gene-info
  $genes = `grep -v '^;' $result_file | grep -v '^#' | cut -f 1 | sort -u `;
  &PipingForm() unless ($query->param("analysis_type") eq "analysis_occ");

  print "<HR SIZE = 3>";

} else {
  &EmailTheResult($command , $query->param('user_email'), $result_file);
}


print $query->end_html;

exit(0);


sub PipingForm {
  ### prepare data for piping
print <<End_of_form;
<CENTER>
<TABLE class="nextstep">
<TR>
  <TD colspan=2>
    <H3>Next step</H3>
  </TD>
  </TR>
  <TR>
  <TD valign=top>
    <FORM METHOD="POST" ACTION="feature-map_form.cgi">
    <INPUT type="hidden" NAME="feature_file" VALUE="$result_file">
    <INPUT type="hidden" NAME="format" VALUE="feature-map">
    <INPUT type="hidden" NAME="handle" VALUE="none">
    <INPUT type="hidden" NAME="fill_form" VALUE="on">
    <INPUT type="submit" value="feature map">
    </FORM>
  </TD>
  <TD>
    <FORM METHOD="POST" ACTION="gene-info_form.cgi">
    <INPUT type="hidden" NAME="queries" VALUE="$genes">
    <INPUT type="submit" value="gene information"> Specify the source organism of the scanned sequences<br/>
End_of_form
	&OrganismPopUp;
	print '
    
    </FORM>
  </TD>
</TR>
</TABLE>
</CENTER>';
#End_of_form
}


################################################################
#
# read matrix scanning parameters
#
sub ReadMatrixScanParameters {

  ################################################################
  ## matrix specification
  unless ($query->param('matrix') =~ /\S/) { ### empty matrix
    &RSAT::error::FatalError("You did not enter the matrix");
  }
  $matrix_file = $tmp_file_path.".matrix";

  $matrix_format = lc($query->param('matrix_format'));
  $parameters .= " -matrix_format ".$matrix_format;

  open MAT, "> ".$matrix_file;
  print MAT $query->param('matrix');
  close MAT;
  &DelayedRemoval($matrix_file);
  $parameters .= " -m $matrix_file";
  push @result_files, ("Matrix file (".$matrix_format.")",$matrix_file);

  ## Use consensus as matrix name
  if ($query->param("consensus_as_name") eq "on") {
    $parameters .= " -consensus_name ";
  }


  ################################################################
  ## Pseudo-counts
  if (&IsReal($query->param('pseudo_counts'))) {
    $parameters .= " -pseudo ".$query->param('pseudo_counts');
  } else {
    &FatalError("Pseudo-count should be a real number");
  }
  if ($query->param('pseudo_distribution') eq "equi_pseudo") {
    $parameters .= " -equi_pseudo ";
  }
  #   ################################################################
  #   ## pseudo-counts and weights are mutually exclusive
  #   if (&IsReal($query->param('pseudo_counts'))) {
  #     $parameters .= " -pseudo ".$query->param('pseudo_counts');
  #   }
  #   if ($query->param('pseudo_distribution') eq "equi_pseudo") {
  #     $parameters .= " -equi_pseudo ";
  #   }

  ################################################################
  ## decimals
  if (&IsReal($query->param('decimals'))) {
    $parameters .= " -decimals ".$query->param('decimals');
  }

  ################################################################
  ## strands
  my $str = "-1str";
  if ($query->param('strands') =~ /both/i) {
    $str = "-2str";
    $parameters .= " -2str";
  } else {
    $str = "-1str";
    $parameters .= " -1str";
  }

  ################################################################
  #### origin
  $origin = $query->param('origin');
  $parameters .= " -origin ".$origin;

  ################################################################
  #### Offset
  my $offset = $query->param('offset');
  if ((&IsInteger($offset)) && ($offset != 0)) {
    $parameters .= " -offset ".$offset;
  }


  ################################################################
  ## Markov order
  my $markov_order = $query->param('markov_order');
  &RSAT::error::FatalError("Markov model should be a Natural number") unless &IsNatural($markov_order);

  ################################################################
  ## Background model method
  my $bg_method = $query->param('bg_method');
  if ($bg_method eq "bginput") {
    $parameters .= " -bginput";
    $parameters .= " -markov ".$markov_order;

  } elsif ($bg_method eq "window") {
    my $window_size = $query->param('window_size');
    &RSAT::message::FatalError(join("\t",$window_size, "Invalid value for the window size. Should be a Natural number." )) unless (&IsNatural($window_size));

    $parameters .= " -window ".$window_size;
    $parameters .= " -markov ".$markov_order;

  } elsif ($bg_method eq "bgfile") {
    ## Select pre-computed background file in RSAT genome directory
    my $organism_name = $query->param("organism");
    my $noov = "ovlp";
    my $background_model = $query->param("background");
    my $oligo_length = $markov_order + 1;
    $bg_file = &ExpectedFreqFile($organism_name,
				 $oligo_length, $background_model,
				 noov=>$noov, str=>"-1str");
    $parameters .= " -bgfile ".$bg_file;

  } elsif ($bg_method =~ /upload/i) {
    ## Upload user-specified background file
    my $bgfile = $tmp_file_path."_bgfile.txt";
    push @result_files, ("Background model", $bgfile);
    my $upload_bgfile = $query->param('upload_bgfile');
    if ($upload_bgfile) {
      if ($upload_bgfile =~ /\.gz$/) {
	$bgfile .= ".gz";
      }
      my $type = $query->uploadInfo($upload_bgfile)->{'Content-Type'};
      open BGFILE, ">$bgfile" ||
	&cgiError("Cannot store background file in temp dir.");
      while (<$upload_bgfile>) {
	print BGFILE;
      }
      close BGFILE;
      $parameters .= " -bgfile $bgfile";
      $parameters .= " -bg_format ".$query->param('bg_format');
    } else {
      &FatalError ("If you want to upload a background model file, you should specify the location of this file on your hard drive with the Browse button");
    }
		
  } else {
    &RSAT::error::FatalError($bg_method," is not a valid method for background specification");
  }

  ################################################################
  ## bg_pseudo
  if (&IsReal($query->param('bg_pseudo'))) {
    $parameters .= " -bg_pseudo ".$query->param('bg_pseudo');
  }

  ################################################################
  ## Return fields

  ################################################################
  ## Return sites
  if ($query->param('analysis_type') eq "analysis_sites") {

    if ($query->param("return_site_limits") eq "on") {
      $parameters .= " -return limits";
    }

    if ($quick_mode){
      if ($query->param("return_field")){
	$parameters .= " -return ".$query->param("return_field");
      }
      my $th = $query->param("thresh_value");
      &RSAT::error::FatalError($th." is not a valid value for the threshold. Should be a number. ") unless (&IsReal($th));
      if ($query->param("return_field") eq "sites"){
	$parameters .= " -lth score ".$th;
      }
      if ($query->param("return_field") eq "pval"){
	$parameters .= " -uth pval ".$th;
      }
    } else {

      my @return_fields = qw(sites pval rank normw weight_limits bg_residues);
      foreach my $field (@return_fields) {
	if ($query->param("return_".$field) eq "on") {
	  $parameters .= " -return ".$field;
	}
      }
    

    ## thresholds
    my @threshold_fields = qw(score pval sig rank proba_M proba_B normw);
    foreach my $field (@threshold_fields) {
      if ($query->param("lth_".$field) ne "none") {
	my $lth = $query->param("lth_".$field);
	&RSAT::error::FatalError($lth." is not a valid value for the lower $field threshold. Should be a number. ") unless (&IsReal($lth));
	$parameters .= " -lth $field $lth ";
      }

      if ($query->param("uth_".$field) ne "none") {
	my $uth = $query->param("uth_".$field);
	&RSAT::error::FatalError($uth." is not a valid value for the upper $field threshold. Should be a number. ") unless (&IsReal($uth));
	$parameters .= " -uth $field $uth ";
      }
    }
  }
  }

  ################################################################
  ## Return distribution
  if ($query->param('analysis_type') eq "analysis_occ") {
    my @return_fields = qw(distrib occ_proba);
    foreach my $field (@return_fields) {
      if ($query->param("return_".$field) eq "on") {
	$parameters .= " -return ".$field;
      }
    }
    	
    if ($query->param('sort_distrib') eq "occ_sig") {
      $parameters .= " -sort_distrib ";
    }
    	
    ## thresholds
    my @threshold_fields = qw(inv_cum exp_occ occ_pval occ_eval occ_sig occ_sig_rank);
    foreach my $field (@threshold_fields) {
      if ($query->param("lth_".$field) ne "none") {
	my $lth = $query->param("lth_".$field);
	&RSAT::error::FatalError($lth." is not a valid value for the lower $field threshold. Should be a number. ") unless (&IsReal($lth));
	$parameters .= " -lth $field $lth ";
      }

      if ($query->param("uth_".$field) ne "none") {
	my $uth = $query->param("uth_".$field);
	&RSAT::error::FatalError($uth." is not a valid value for the upper $field threshold. Should be a number. ") unless (&IsReal($uth));
	$parameters .= " -uth $field $uth ";
      }
    }
    my @threshold_fields = qw(occ_score);
    foreach my $field (@threshold_fields) {
      $ms_field = $field;
      $ms_field =~ s/occ_//;
      if ($query->param("lth_".$field) ne "none") {
	my $lth = $query->param("lth_".$field);
	&RSAT::error::FatalError($lth." is not a valid value for the lower $ms_field threshold. Should be a number. ") unless (&IsReal($lth));
	$parameters .= " -lth $ms_field $lth ";
      }

      if ($query->param("uth_".$field) ne "none") {
	my $uth = $query->param("uth_".$field);
	&RSAT::error::FatalError($uth." is not a valid value for the upper $ms_field threshold. Should be a number. ") unless (&IsReal($uth));
	$parameters .= " -uth $ms_field $uth ";
      }
    }
  }

  ################################################################
  ## Return crer
  if ($query->param('analysis_type') eq "analysis_crer") {
  	$parameters .= " -n skip ";
    my @return_fields = qw(crer normw);
    foreach my $field (@return_fields) {
      if ($query->param("return_".$field) eq "on") {
	$parameters .= " -return ".$field;
      }
    }
    if ($query->param("return_crer_limits") eq "on") {
      $parameters .= " -return limits";
    }
    
    if ($query->param("return_crer_sites") eq "on") {
      $parameters .= " -return sites";
    }
    
    if ($query->param("crer_ids") eq "on") {
      $parameters .= " -crer_ids";
    }
    	
    ## thresholds
    my @threshold_fields = qw(crer_size crer_sites crer_pval crer_sig);
    foreach my $field (@threshold_fields) {
      if ($query->param("lth_".$field) ne "none") {
	my $lth = $query->param("lth_".$field);
	&RSAT::error::FatalError($lth." is not a valid value for the lower $field threshold. Should be a number. ") unless (&IsReal($lth));
	$parameters .= " -lth $field $lth ";
      }

      if ($query->param("uth_".$field) ne "none") {
	my $uth = $query->param("uth_".$field);
	&RSAT::error::FatalError($uth." is not a valid value for the upper $field threshold. Should be a number. ") unless (&IsReal($uth));
	$parameters .= " -uth $field $uth ";
      }
    }

    my @threshold_fields = qw(site_pval);
    foreach my $field (@threshold_fields) {
      $ms_field = $field;
      $ms_field =~ s/site_//;
      if ($query->param("lth_".$field) ne "none") {
	my $lth = $query->param("lth_".$field);
	&RSAT::error::FatalError($lth." is not a valid value for the lower $ms_field threshold. Should be a number. ") unless (&IsReal($lth));
	$parameters .= " -lth $ms_field $lth ";
      }

      if ($query->param("uth_".$field) ne "none") {
	my $uth = $query->param("uth_".$field);
	&RSAT::error::FatalError($uth." is not a valid value for the upper $ms_field threshold. Should be a number. ") unless (&IsReal($uth));
	$parameters .= " -uth $ms_field $uth ";
      }
    }
  }

  my @return_fields = qw(matrix freq_matrix weight_matrix bg_model);
  foreach my $field (@return_fields) {
    if ($query->param("return_".$field) eq "on") {
      $parameters .= " -return ".$field;
    }
  }
 
}
