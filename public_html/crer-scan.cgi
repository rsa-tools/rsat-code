#!/usr/bin/perl
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

### Print the header
&RSA_header("crer-scan result", "results");
&ListParameters() if ($ENV{rsat_echo} >= 2);

## Check security issues
&CheckWebInput($query);

$command = "python3 ".$ENV{RSAT}."/python-scripts/crer_scan.py";
$prefix = "crer_scan";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); ($tmp_file_dir, $tmp_file_name) = &SplitFileName($tmp_file_path);
@result_files = ();

#### read parameters ####
$parameters = " -v 1 -s ";

################################################################
## Sites
if ($query->param('sites') =~ /\S/) {
  $query_file = $tmp_file_path."_query.txt";
  push @result_files, ("input sites",$query_file);
  open QUERY, ">".$query_file;
  print QUERY $query->param('sites');
  close QUERY;
  &DelayedRemoval($query_file);
  $parameters .= " -i ".$query_file;
} else {
    &cgiError("Please enter sites.\n");
}

################################################################
## site format
if ($query->param('in_format')) {
    my ($in_format) = split " ", $query->param('in_format'); ### take the first word
    $parameters .= " -in_format ".$in_format;
}

## Lower threshold on CRER size
if (&IsNatural($query->param('lth_crer_size'))) {
  $parameters .= " -lth_crer_size ".$query->param('lth_crer_size');
}
## Upper threshold on CRER size
if (&IsNatural($query->param('uth_crer_size'))) {
  $parameters .= " -uth_crer_size ".$query->param('uth_crer_size');
}

## Lower threshold on number of sites per CRER
if (&IsNatural($query->param('lth_crer_sites'))) {
  $parameters .= " -lth_crer_sites ".$query->param('lth_crer_sites');
}
## Upper threshold on number of sites per CRER
if (&IsNatural($query->param('uth_crer_sites'))) {
  $parameters .= " -uth_crer_sites ".$query->param('uth_crer_sites');
}

## Lower threshold on distance between sites
if (&IsNatural($query->param('lth_crer_sites_distance'))) {
  $parameters .= " -lth_crer_sites_distance ".$query->param('lth_crer_sites_distance');
}
## Upper threshold on distance between sites
if (&IsNatural($query->param('uth_crer_sites_distance'))) {
  $parameters .= " -uth_crer_sites_distance ".$query->param('uth_crer_sites_distance');
}

## Lower threshold on site p-value
if (&IsReal($query->param('lth_site_pval'))) {
  if ($query->param('lth_site_pval') < 0) {
    &RSAT::error::FatalError($query->param('lth_site_pval'), "is not a valid threshold on p-value (must be positive)");
  } elsif ($query->param('lth_site_pval') > 1) {
    &RSAT::error::FatalError($query->param('lth_site_pval'), "is not a valid threshold on p-value (must be <= 1)");
  }
  $parameters .= " -lth_site_pval ".$query->param('lth_site_pval');
}
## Upper threshold on site p-value
if (&IsReal($query->param('uth_site_pval'))) {
  if ($query->param('uth_site_pval') < 0) {
    &RSAT::error::FatalError($query->param('uth_site_pval'), "is not a valid threshold on p-value (must be positive)");
  } elsif ($query->param('uth_site_pval') > 1) {
    &RSAT::error::FatalError($query->param('uth_site_pval'), "is not a valid threshold on p-value (must be <= 1)");
  }
  $parameters .= " -uth_site_pval ".$query->param('uth_site_pval');
}

## Lower threshold on CRER e-value
if (&IsReal($query->param('lth_crer_eval'))) {
  if ($query->param('lth_crer_eval') < 0) {
    &RSAT::error::FatalError($query->param('lth_crer_eval'), "is not a valid threshold on e-value (must be positive)");
  }
  $parameters .= " -lth_crer_eval ".$query->param('lth_crer_eval');
}
## Upper threshold on CRER e-value
if (&IsReal($query->param('uth_crer_eval'))) {
  if ($query->param('uth_crer_eval') < 0) {
    &RSAT::error::FatalError($query->param('uth_crer_eval'), "is not a valid threshold on e-value (must be positive)");
  }
  $parameters .= " -uth_crer_eval ".$query->param('uth_crer_eval');
}

## Lower threshold on CRER significance
if (&IsReal($query->param('lth_crer_sig'))) {
  $parameters .= " -lth_crer_sig ".$query->param('lth_crer_sig');
}
## Upper threshold on CRER significance
if (&IsReal($query->param('uth_crer_sig'))) {
  $parameters .= " -uth_crer_sig ".$query->param('uth_crer_sig');
}


## Lower threshold on CRER score
if (&IsReal($query->param('lth_score'))) {
  $parameters .= " -lth_score ".$query->param('lth_score');
}
## Upper threshold on CRER score
if (&IsReal($query->param('uth_score'))) {
  $parameters .= " -uth_score ".$query->param('uth_score');
}

## Lower threshold on overlap between sites
if (&IsNatural($query->param('lth_overlap'))) {
  $parameters .= " -lth_overlap ".$query->param('lth_overlap');
}
## Upper threshold on overlap between sites
if (&IsNatural($query->param('uth_overlap'))) {
  $parameters .= " -uth_overlap ".$query->param('uth_overlap');
}

## Only return limits of sequences containing at least one CRER
$parameters .= " -return_limits_filtered";


## Output file
$result_file = $tmp_file_path.".tab";
#$parameters .= " -o ".$result_files;
push @result_files, ("CRERs",$result_file);

&ReportWebCommand($command." ".$parameters);


## Update log file
&UpdateLogFile();

################################################################
## Run the command
if ($query->param('output') eq "display") {
  ## Print a warning telling that the results will appear below
  &PipingWarning();
  
  ## Run the command, print the result on Web page, and store a copy in result file
  open RESULT, "$command $parameters |";
  print '<H2>Result</H2>';
  &PrintHtmlTable(RESULT, $result_file, 1);
  close(RESULT);
  print "<HR SIZE = 3>";

  &PrintURLTable(@result_files);

  ## Collect genes for piping the results to gene-info
  $genes = `grep -v '^;' $result_file | grep -v '^#' | cut -f 1 | sort -u `;
  &PipingForm();

  print "<HR SIZE = 3>";

} else {
  &EmailTheResult("$command $parameters", $query->param('user_email'), $result_file);
}
print $query->end_html();

exit(0);


##############f##################################################
##
## Pipe the result to other commands
##
sub PipingForm {
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

