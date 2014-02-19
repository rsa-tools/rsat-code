#!/usr/bin/perl
############################################################
#
# info-gibbs
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

### print the result page
&RSA_header("info-gibbs result", "results");
&ListParameters() if ($ENV{rsat_echo} >=2);


## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

@result_files = ();
$command = "$ENV{RSAT}/contrib/info-gibbs/info-gibbs";
#$command = "$ENV{RSAT}/python-scripts/info-gibbs-python";
#$convert_matrix_command = "$SCRIPTS/convert-matrix -from gibbs -return counts";
$convert_seq_command = "$SCRIPTS/convert-seq";
$prefix = "info-gibbs";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); ($tmp_file_dir, $tmp_file_name) = &SplitFileName($tmp_file_path);
#$tmp_file_name = sprintf "info-gibbs.%s", &AlphaDate();

################################################################
#
# read parameters
#
$parameters = '';

## Scan single or both strands
$strand="-1str";

#### strand
if (lc($query->param("two_strands")) eq "on") {
  $strand="-2str";
  $parameters .= ' --strand=+-';
} else {
  $parameters .= ' --strand=+';
}


## sequence file
($sequence_file,$sequence_format) = &GetSequenceFile("fasta", no_format=>1, add_rc=>0);
$parameters .= " -i ".$sequence_file;
push @result_files, ("Input sequence",$sequence_file);

### matrix length
if (&IsNatural($query->param('length'))) {
    $parameters .= ' -w ' . $query->param('length');
}

### expected number of matches
# if (&IsNatural($query->param('expected'))) {
#     $parameters .= ' --words='.$query->param('expected');
# }

if (&IsNatural($query->param('expected'))) {
    $parameters .= ' --mean_sps='.$query->param('expected');
}

if (&IsNatural($query->param('iter'))) {
    $parameters .= ' --iter='.$query->param('iter');
}
if (&IsNatural($query->param('nrun'))) {
    $parameters .= ' --nrun='.$query->param('nrun');
}

if (&IsNatural($query->param('motifs'))) {
    $parameters .= ' --motifs='.$query->param('motifs');
}

if ($query->param('freq_estimate') =~ /background/i) {
  if (&IsNatural($query->param('bg_order'))) {
    $oligo_length = $query->param('bg_order') + 1;
  }
  #    $oligo_length = 3 + 1; # MM3
  %supported_background = (
			   "upstream"=>1,
			   "upstream-noorf"=>1,
#			   "intergenic"=>1
			  );

  ### check organism
  unless ($organism = $query->param('organism')) {
    &cgiError("You should specify an organism to use intergenic frequency calibration");
  }
  unless (%{$supported_organism{$organism}}) {
    &cgiError("Organism $org is not supported on this site");
  }
  my $background = $query->param("background");
  unless ($supported_background{$background}) {
    &cgiError("$background is not supported as background model");
  }
  ### Taxon-specific background model
  if ($query->param('bg_level') eq 'taxon') {
    unless ($taxon = $query->param('taxon')) {
      &cgiError("You should specify an taxon to use intergenic frequency calibration");
    }
    &CheckTaxon($taxon);
    $organism = $taxon;
  }

  ################################################################
  ## Convert the background model because info-gibbs requires bg in MotifSampler (inclusive) format
  #$exp_freq_file = "$ENV{RSAT}/public_html/data/genomes/$organism/oligo-frequencies/" . "$oligo_length" . "nt_" . "$background" . "_" . "$organism$overlap$strand.freq.gz";
  $exp_freq_file = &ExpectedFreqFile($organism, $oligo_length, $background, type=>$oligotype, noov=>$overlap, str=>$strand, taxon=>$taxon);
  $convert_bg_cmd = "$SCRIPTS/convert-background-model -from oligo-analysis -to MotifSampler -i $exp_freq_file -o ${TMP}/$tmp_file_name.bg";
  &ReportWebCommand($convert_bg_cmd);
  system "$convert_bg_cmd";

  $parameters .= "--bgfile=${TMP}/$tmp_file_name.bg ";
  #print $exp_freq_file;
  #$freq_option = " -bg $background -org $organism";
  #$freq_option = " --bgoligo=$exp_freq_file.gz";


  # } elsif ($query->param('freq_estimate') =~ /upload/i) {
  #     $exp_freq_file = "${TMP}/$tmp_file_name.expfreq";
  #     $upload_freq_file = $query->param('upload_freq_file');
  #     if ($upload_freq_file) {
  #       ## Support compressed .gz files
  #       if ($upload_freq_file =~ /\.gz$/) {
  #   $exp_freq_file .= ".gz";
  #       }
  #       $type = $query->uploadInfo($upload_freq_file)->{'Content-Type'};
  #       open FREQ, ">$exp_freq_file" ||
  #   &cgiError("Cannot store expected frequency file in temp dir.");
  #       while (<$upload_freq_file>) {
  #   print FREQ;
  #       }
  #       close FREQ;
  #       $freq_option = " --bgoligo=$exp_freq_file";
  #     } else {
  #       &FatalError ("If you want to upload an expected frequency file, you should specify the location of this file on your hard drive with the Browse button");
  #     }
}

## Result file
$result_file = $tmp_file_path.".tab";
push @result_files, ('info-gibbs result', $result_file);

### additional parameters
#$parameters .= ' --finalcycle';
&ReportWebCommand($command." ".$parameters);

if ($query->param('output') eq "display") {  
    &PipingWarning();

    ### execute the command ###
    #$matrix_file = &RSAT::util::get_pub_temp()."/$tmp_file_name.matrix";
    #print("$command $parameters\n");
    system "$command $parameters > $result_file";
    &DelayedRemoval($result_file);

    ### Print result on the web page
    print '<h4>info-gibbs result</h4>';
    print "<pre>";
    my ($res) = &OpenInputFile($result_file);
    while (<$res>) {
      s|$ENV{RSAT}/||g;
      print $_;
    }
    print "</pre>";

    ## Display matrices with logos and links
    my ($out_matrix_file) = &display_matrices_web($result_file, "infogibbs");
    push @result_files, ('converted matrices', $out_matrix_file);

    &PrintURLTable(@result_files);
    &PipingForm();

    print "<hr size=\"3\">";

} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'), $result_file);
}

print $query->end_html;
exit(0);



### prepare data for piping
sub PipingForm {
  my $command = "$ENV{RSAT}/perl-scripts/convert-matrix -i $result_file -from tab -to tab -top 1 -return counts";
  my $matrix_content = `$command`;
  $matrix_content =~ s|//\n||gm;
  $matrix_content =~ s|;.*\n||gm;
#  print "<pre>".$command."</pre>";
#  print "<pre>".$matrix_content."</pre>";

  $title = $query->param('title');
  $title =~ s/\"/\'/g;
    print <<End_of_form;
<hr size="3">
<table class="Nextstep">
<tr>
<td colspan="3">
<h3>next step</h3>
</td>
</tr><tr>
<!--
<td valign="bottom" align="center">
<form method="post" action="patser_form.cgi">
<input type="hidden" name="title" value="$title">
<input type="hidden" name="matrix_file" value="$result_file">
<input type="hidden" name="matrix_format" value="tab">
<input type="hidden" name="sequence_file" value="$sequence_file">
<input type="hidden" name="sequence_format" value="$sequence_format">
<input type="submit" value="pattern matching (patser)">
</form>
</td>
-->
<td valign="bottom" align="center">
<b><font color=red>new</a></b>
<form method="POST" action="matrix-scan_form.cgi">
<input type="hidden" name="title" value="$title">
<input type="hidden" name="matrix_file" value="$result_file">
<input type="hidden" name="matrix_format" value="tab">
<input type="hidden" name="sequence_file" value="$sequence_file">
<input type="hidden" name="sequence_format" value="fasta">
<input type="submit" value="pattern matching (matrix-scan)">
</form>
</td>

<td valign=bottom align=center>
<form method="post" action="convert-matrix_form.cgi">
<input type="hidden" name="title" value="$title">
<input type="hidden" name="matrix_file" value="$result_file">
<input type="hidden" name="matrix_format" value="tab">
<input type="hidden" name="logo" value="on" checked="checked">
<input type="submit" value="convert-matrix">
</form>
</td>


<td valign=bottom align=center>
<form method="post" target='_blank' action="http://meme.nbcr.net/meme4_3_0/cgi-bin/tomtom.cgi">
<input type="hidden" name="query" value="$matrix_content">
<input type="hidden" name="DIST" value="sandelin">
<input type="submit" value="TOMTOM">
</form>
Compare a single matrix to a motif database.
</td>
</tr>


</tr>
</table>
End_of_form
}
