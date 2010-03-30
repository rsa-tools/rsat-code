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
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
############################################ configuration
$command = "$ENV{RSAT}/perl-scripts/chip-seq-analysis";
$output_directory = sprintf "chip-seq-analysis.%s", &AlphaDate();
$output_prefix = "ChIP-seq_analysis_";
$output_path = "$TMP/$output_directory";
$output_path =~ s|\/\/|\/|g;
`mkdir -p $output_path`;
# $ENV{'PATH'} = $ENV{'PATH'} . ":$ENV{RSAT}/perl-scripts" . ":$ENV{RSAT}/python-scripts";
# $ENV{'PATH'} = $ENV{'PATH'} . ":$ENV{RSAT}/bin";

############################################ result page header
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
($sequence_file, $sequence_format) = &MultiGetSequenceFile(1, "$output_path/$output_prefix" . "peak_seq", 1);

### control sequences file
($control_sequence_file, $control_sequence_format) = &MultiGetSequenceFile(2, "$output_path/$output_prefix" . "ctl_seq", 0);

#print "<pre>sequence_file: [$control_sequence_file]\n</pre>" if ($ENV{rsat_echo} >=1);

### peak only or peak+control analysis ?
if ($control_sequence_file eq '') {
    $analysis = 'peaks';
    $parameters .= "-i $sequence_file";
} else {
    $analysis = 'peaks+control';
    $parameters .= "-i $sequence_file -ctl $control_sequence_file";
}

### tasks
@tasks = ("purge", "seqlen", "profiles", "synthesis");

if ($analysis eq "peaks") {
    if ($query->param('oligo-analysis') =~ /on/) {
        push(@tasks, "oligos");
        #$oligo_length = $query->param('oligo_length');
        #$oligo_length = $query->param('oligo_length');
        #&FatalError("$oligo_length Invalid oligonucleotide length") unless &IsNatural($oligo_length);
        #$parameters .= "-l $oligo_length";
    }
    if ($query->param('dyad-analysis') =~ /on/) {
        push(@tasks, "dyads");
    }
    if ($query->param('local-word-analysis') =~ /on/) {
        push(@tasks, "local_words");
    }
} else {
    push(@tasks, "oligos-diff");
}

push(@tasks, "positions");
push(@tasks, "word_compa");
push(@tasks, "motif_compa");

### add -task
$parameters .= " -task " . join(",", @tasks);

### task specific parameters
if (&IsNatural($query->param('markov_order'))) {
  $param .= " -markov ".$query->param('markov_order');
}

### output directory
$parameters .= " -outdir $output_path";

### output prefix
$parameters .= " -prefix $output_prefix";

### verbosity
$parameters .= " -v 1";

############################################ display or send result
$index_file = $output_directory."/".$output_prefix."synthesis.html";
#$index_file = $output_directory."/".$output_prefix."index.html";
my $mail_title = join (" ", "[RSAT]", "chip-seq-analysis", &AlphaDate());
&EmailTheResult("$command $parameters", $query->param('user_email'), $index_file, title=>$mail_title);
#&ServerOutput("$command $parameters", $query->param('user_email'), $tmp_file_name);
# $debug = "$command $parameters 2> $TMP/log.txt";
# print $debug;
# `$debug`;

############################################ result page footer
print $query->end_html;

exit(0);

