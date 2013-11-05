#!/usr/bin/perl
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

### Read the CGI query
$query = new CGI;

### Print the header
&RSA_header("pattern-assembly result", "results");

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

&ListParameters if ($ENV{rsat_echo} >=2);

@result_files = ();

$pattern_assembly_command = "$SCRIPTS/pattern-assembly";
$prefix = "pattern-assembly";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); ($tmp_file_dir, $tmp_file_name) = &SplitFileName($tmp_file_path);
#$tmp_file_name = sprintf "pattern-assembly.%s", &AlphaDate;

#### read parameters ####
$parameters .= " -v 1";

################################################################
#### patterns
my $pattern_file = $tmp_file_path.".pat";
push @result_files, "Patterns", $pattern_file;
open PATTERNS, "> ".$pattern_file;
if ( $query->param('patterns') =~ /\S/) {
    print PATTERNS $query->param('patterns');
} elsif ($query->param('pattern_file') =~ /\S/) {
    ### upload file from the client
    $upload_data_file = $query->param('pattern_file');
    $type = $query->uploadInfo($upload_data_file)->{'Content-Type'};
    while (<$upload_data_file>) {
	print PATTERNS;
    }
} else {
    &cgiError("You should either enter query patterns in the box, or specify a pattern file\n");
}
close PATTERNS;
&DelayedRemoval($pattern_file);
$parameters .= " -i ".$pattern_file;

################################################################
####  max number of flanking residues
$maxfl = $query->param('maxfl');
if (&IsNatural($maxfl)) {
    if ($maxfl < 0) {
	&cgiError("Maximum flanking residues should be positive ($maxfl is invalid)");
    } else {
	$parameters .= " -maxfl $maxfl";
    }
} else {
    &cgiError("Maximum flanking residues should be a positive natural number ($maxfl is invalid)");
}

################################################################
####  max number of substitutions
$subst = $query->param('subst');
if (&IsNatural($subst)) {
    if ($subst < 0) {
	&cgiError("Maximum substitutions should be positive ($subst is invalid)");
    } else {
	$parameters .= " -subst $subst";
    }
} else {
    &cgiError("Maximum substitutions should be a positive natural number ($subst is invalid)");
}

################################################################
####  max number of patterns
$maxpat = $query->param('maxpat');
if (&IsNatural($maxpat)) {
    if ($maxpat < 0) {
	&cgiError("Maximum number of patterns should be positive ($maxpat is invalid)");
    } else {
	$parameters .= " -maxpat $maxpat";
    }
} else {
    &cgiError("Maximum number of patterns should be a positive natural number ($maxpat is invalid)");
}

################################################################
####  max number of patterns
$sc = $query->param('sc');
if (&IsNatural($sc)) {
    if ($sc < 2) {
	&cgiError("Score column should be >= 2 ($sc is invalid)");
    } else {
	$parameters .= " -sc $sc";
    }
}

################################################################
#### single or double strand assembly
if ($query->param('strand') =~ /insensitive/) {
  $parameters .= " -2str";
} else {
  $parameters .= " -1str";
}

## Output file
$result_file = $tmp_file_path.".asmb";
push @result_files, "Assembly", $result_file;

$pattern_assembly_command .= " ".$parameters;

################################################################
#### run the command
&ReportWebCommand($pattern_assembly_command);
if ($query->param('output') eq "display") {
  open RESULT, "$pattern_assembly_command |";
  print '<H2>Result</H2>';
  &PrintHtmlTable(RESULT, $result_file, true);
  close(RESULT);
  &PrintURLTable(@result_files);

  print "<HR SIZE = 3>";

} else { 
    &EmailTheResult($pattern_assembly_command, $query->param('user_email'), $result_file);
}
print $query->end_html();

exit(0);


