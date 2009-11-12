#!/usr/bin/perl
############################################################
#
# $Id: matrix-scan.cgi,v 1.29 2009/11/12 09:33:18 jvanheld Exp $
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
#    $ERR_LOG = "$TMP/RSA_ERROR_LOG.txt";
    use CGI::Carp qw(carpout);
    open (LOG, ">> $ERR_LOG")
	|| die "Unable to redirect log\n";
    carpout(*LOG);
}
require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

$command = $SCRIPTS."/matrix-scan -v 1 ";
$tmp_file_name = sprintf "matrix-scan.%s", &AlphaDate();
$result_file =  $TMP."/".$tmp_file_name.".ft";

### Read the CGI query
$query = new CGI;

#$ENV{rsat_echo}=1;

### print the result page
&RSA_header("matrix-scan result", "results");
&ListParameters() if ($ENV{rsat_echo} >=2);

#### update log file ####
&UpdateLogFile();

################################################################
## sequence file
($sequence_file,$sequence_format) = &GetSequenceFile();

#### matrix-scan parameters
&ReadMatrixScanParameters();
$parameters .= " -i $sequence_file -seq_format $sequence_format";

################################################################
### vertically print the matrix
if ($query->param('vertically_print')) {
  $parameters .= " -p";
}


################################################################
## Treatment of N characters
$parameters .= ' -n score';

$command .= " $parameters";
#$command .= " -o ".$result_file;

################################################################
#### echo the command (for debugging)
print "<pre>$command</pre>" if ($ENV{rsat_echo} >= 1);

################################################################
### execute the command ###
if ($query->param('output') eq "display") {

  unless ($query->param('table')) {
    &PipingWarning() unless ($query->param("analysis_type") eq "analysis_occ");
  }



  ### Print the result on Web page
  open RESULT, "$command |";
  print "<PRE>";
  &PrintHtmlTable(RESULT, $result_file, true);
  print "</PRE>";
  print "<HR SIZE = 3>";

#  system("$command");


  ################################################################
  ## Convert features to GFF3 format
  $gff3_file =  $TMP."/".$tmp_file_name.".gff3";
  $command = "${SCRIPTS}/convert-features -from ft -to gff3 ";
  $command .= " -i ".$result_file;
  $command .= " -o ".$gff3_file;
  &RSAT::message::Info("Converting features to GFF3 format", $command) if ($ENV{rsat_echo} >=2);
  system($command);

  ################################################################
  ## Convert features for loading as genome tracks
  if ($origin eq "genomic") {
    $bed_file =  $TMP."/".$tmp_file_name.".bed";
    $command = "${SCRIPTS}/convert-features -from ft -to bed ";
    $command .= " -i ".$result_file;
    $command .= " -o ".$bed_file;
    &RSAT::message::Info("Converting features to BED format", $command) if ($ENV{rsat_echo} >=2);;
    system($command);
  }


  ################################################################
  ## Table with links to the raw result files in various formats
  $result_URL = $ENV{rsat_www}."/tmp/".$tmp_file_name;
  print "<center><TABLE class=\"nextstep\">\n";
  print "<tr><td colspan='3'><h3>Raw result files</h3> </td></tr>";
#  print ("<tr>",
#	 "<th>Format</th>",
#	 "<th>URL</th>",
#	 "<th>Usage</th>",
#	 "</tr>");
  print ("<tr>",
	 "<td>FT</td>",
	 "<td>","<a href='".$result_URL.".ft'>".$result_URL.".ft</a>","</td>",
	 "<td>RSAT feature-map</td>",
	 "</tr>");
  print ("<tr>",
	 "<td>GFF3</td>",
	 "<td>","<a href='".$result_URL.".gff3'>".$result_URL.".gff3</a>","</td>",
	 "<td>General feature format</td>",
	 "</tr>");
  if ($origin eq "genomic") {
    ################################################################
    ## Identify the source of the features
    my $genomic_format = `grep 'Genomic coordinate format' $result_file`;
    chomp($genomic_format);
    if ($genomic_format =~ /Genomic coordinate format\s+(\S+)/) {
      $genomic_format = lc($1);
      &RSAT::message::Warning("Genomic format: $genomic_format") if ($ENV{rsat_echo} >=2);;
    }
    if ($genomic_format eq "ucsc") {
      $browser_url = "<a target='_blank' href='http://genome.ucsc.edu/cgi-bin/hgCustom'>UCSC genome browser</a>",
    } elsif ($genomic_format =~ /ensembl/i) {
      my $browser = `grep -i 'browser url' $result_file`;
      chomp($browser);
      $browser =~ s/.*browser url\s+(\S+)/$1/i;
      &RSAT::message::Warning($browser);
      $browser_url = "<a target='_blank' href='";
      $browser_url .= $browser;
      $browser_url .= ";contigviewbottom=url:".$result_URL.".bed";
      $browser_url .= "=normal'>EnsEMBL genome browser</a>";
    }

    print ("<tr>",
	   "<td>BED</td>",
	   "<td>","<a href='".$result_URL.".bed'>".$result_URL.".bed</a>","</td>",
	   "<td>".$browser_url."</td>",
	   "</tr>");
  }
  print "</table></center>";

  ## Collect genes for piping the results to gene-info
  $genes = `grep -v '^;' $result_file | grep -v '^#' | cut -f 1 | sort -u `;
  &PipingForm() unless ($query->param("analysis_type") eq "analysis_occ");

  print "<HR SIZE = 3>";

#} elsif ($query->param('output') =~ /server/i) {
  #&ServerOutput("$command ", $query->param('user_email'));
} else {
  &EmailTheResult("$command ", $query->param('user_email'));
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
# read patser parameters
#
sub ReadMatrixScanParameters {

  ################################################################
  ## matrix specification
  unless ($query->param('matrix') =~ /\S/) { ### empty matrix
    &RSAT::error::FatalError("You did not enter the matrix");
  }
  $matrix_file = "$TMP/$tmp_file_name.matrix";

  $matrix_format = lc($query->param('matrix_format'));
  $parameters .= " -matrix_format ".$matrix_format;

  open MAT, "> ".$matrix_file;
  print MAT $query->param('matrix');
  close MAT;
  &DelayedRemoval($matrix_file);
  $parameters .= " -m $matrix_file";

  ## Use consensus as matrix name
  if ($query->param("consensus_as_name") eq "on") {
    $parameters .= " -consensus_name ";
  }

  ################################################################
  ## pseudo-counts and weights are mutually exclusive
  if (&IsReal($query->param('pseudo_counts'))) {
    $parameters .= " -pseudo ".$query->param('pseudo_counts');
  }
 
  if ($query->param('pseudo_distribution') eq "equi_pseudo") {
    $parameters .= " -equi_pseudo ";
  }
 
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
    my $bgfile = "${TMP}/${tmp_file_name}_bgfile.txt";
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
    my @return_fields = qw(sites pval rank normw weight_limits bg_residues);
    foreach my $field (@return_fields) {
      if ($query->param("return_".$field) eq "on") {
	$parameters .= " -return ".$field;
      }
    }
    
     if ($query->param("return_site_limits") eq "on") {
      $parameters .= " -return limits";
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
  	$parameters .= " -n score ";
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
