#!/usr/bin/perl

############################################ imports
### redirect error log to a file
if ($0 =~ /([^(\/)]+)$/) {
  push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
### redirect error log to a file
BEGIN {
    $ERR_LOG = "/dev/null";
#    $ERR_LOG = "/tmp/RSA_ERROR_LOG.txt";
    use CGI::Carp qw(carpout);
    open (LOG, ">> $ERR_LOG")
	|| die "Unable to redirect log\n";
    carpout(*LOG);
}
require "RSA.lib";
require "RSA.disco.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

############################################ configuration
$chip_seq_analysis_command = "$ENV{RSAT}/perl-scripts/chip-seq-analysis";
#$convert_seq_command = "$SCRIPTS/convert-seq";
#$purge_sequence_command = "$SCRIPTS/purge-sequence";
$tmp_file_name = sprintf "chip-seq-analysis.%s", &AlphaDate();

############################################ result header
### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("chip-seq_analysis result", "results");
&ListParameters() if ($ENV{rsat_echo} >=2);

### update log file
&UpdateLogFile();

############################################ command line paramters
### read parameters
$parameters = "";

### peak sequences file
($sequence_file, $sequence_format) = &GetSequenceFile(1);

### control sequences file
($control_sequence_file, $control_sequence_format) = &GetSequenceFile(2);

### peak only or peak+control analysis ?
if ($control_sequence_file eq '') {
    $analysis = 'peaks';
    $input_data = "-i $sequence_file";
} else {
    $analysis = 'peaks+control';
    $input_data = "-i $sequence_file -ctl $control_sequence_file";
}

### tasks
@tasks = ();
if ($query->param('oligo-analysis') =~ /yes/ ) {
    push(@tasks, "oligos");
    #$parameters .= '-task oligos ';
    $oligo_length = $query->param('oligo_length');
    &FatalError("$oligo_length Invalid oligonucleotide length") unless &IsNatural($oligo_length);
    $parameters .= "-l $oligo_length";
}
if ($query->param('dyad-analysis') =~ /yes/ ) {
    push(@tasks, "dyads");
    #$parameters .= '-task dyads ';
}
if ($query->param('position-analysis') =~ /yes/ ) {
    push(@tasks, "positions");
    #$parameters .= 'positions ';
}

### add -task
$parameters .= "-task " . join(",", @tasks);

### task specific parameters
if (&IsNatural($query->param('markov_order'))) {
  $param .= " --markov=".$query->param('markov_order');
}

### verbosity
$parameters .= " -v 5";

############################################ construct command
$command .=  "$chip_seq_analysis_command $input_data $parameters";

print "<pre>command: $command<P>\n</pre>" if ($ENV{rsat_echo} >=1);

&SaveCommand("$command", "$TMP/$tmp_file_name");

############################################ display or send result
if ($query->param('output') =~ /server/i) {
    &ServerOutput("$command", $query->param('user_email'), $tmp_file_name);
} else {
    &EmailTheResult("$command", $query->param('user_email'), $tmp_file_name);
}

############################################ result footer
print $query->end_html;

exit(0);

