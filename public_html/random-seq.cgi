#!/usr/bin/perl

## CVS
## added the possibility to specify the expected frequency for each nucleotide separately

if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
#require "cgi-lib.pl";
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
$command = "$SCRIPTS/random-seq";
$tmp_file_name = sprintf "random-seq.%s", &AlphaDate();

$size_limit = 5e+6;

### Read the CGI query
$query = new CGI;

### print the header
&RSA_header("Random sequence result", "results");

#### update log file ####
&UpdateLogFile();


&ListParameters() if ($ENV{rsat_echo} >= 2);

################################################################
#### read parameters ####
$parameters = "";

#### sequence length
$length = $query->param('length');
if (&IsNatural($length)) {
    $parameters .= " -l $length ";
} else {
    &FatalError("Sequence length must be a natural number");
}

#### number of repetitions
$repet = $query->param('repet');
if (&IsNatural($repet)) {
    $parameters .= " -n $repet";
} else {
    &FatalError("Repetitions must be a natural number");
}

#### check lengths and repetitions
if ($repet * $length == 0) {
    &FatalError("Sequence length and repeitions must be non-null");
}
if ($repet*$length > $size_limit) {
    &FatalError("The web interface does not support queries of this size. Maximum size per query (length * repetitions) = ".$size_limit);
}

#### line width
if (&IsNatural($query->param('lw'))) {
    $parameters .= " -lw ".$query->param('lw');
}

#### output format
$out_format = $query->param('format');
&CheckOutputSeqFormat($out_format);
$parameters .= " -format $out_format ";

#### alphabet
if ($query->param('proba') eq "alphabet") {

    $freq{A} = $query->param('Afreq');
    $freq{C} = $query->param('Cfreq');
    $freq{T} = $query->param('Tfreq');
    $freq{G} = $query->param('Gfreq');

    ## Check the values 
    foreach my $letter (keys %freq) {
	unless (&IsReal($freq{$letter})) {
	    &FatalError("Invalid frequency value ".$freq{$letter}." for residue ".$letter);
	}
    }

    ## Print residue frequencies in a file
    $alphabet_file = $TMP."/".$tmp_file_name.".alphabet";
    open ALPHA, ">".$alphabet_file;
    foreach my $letter (keys %freq) {
	print ALPHA $letter, "\t", $freq{$letter}, "\n";
    }
    close ALPHA;
    $parameters .= " -expfreq ".$alphabet_file;
    &DelayedRemoval($alphabet_file);

## Pre-calibrated Markov models
} elsif (($query->param('proba') =~ /upstream/i) ||
	 ($query->param('proba') =~ /protein/i)) {
    ### check organism
    unless ($organism = $query->param('organism')) {
	&cgiError("You should specify an organism to use upstream frequency calibration");
    }
    unless (defined(%{$supported_organism{$organism}})) {
	&cgiError("Organism $organism is not supported on this site");
    }
    if ($query->param('proba') =~ /protein/i) {
      $oligopept_size = $query->param("oligopept_size");
      unless (&IsNatural($oligopept_size)) {
	&cgiError("Invalid oligopeptide length $oligopept_size");
      }
      $seq_type = "protein"; ## Used for the piping form
      $parameters .= " -bg protein -org $organism -ol $oligopept_size -type protein";
    } else {
      $oligo_size = $query->param("oligo_size");
      unless (&IsNatural($oligo_size)) {
	&cgiError("Invalid oligonucleotide length $oligo_size");
      }
      $seq_type = "dna"; ## Used for the piping form
      $parameters .= " -bg upstream-noorf -org $organism -ol $oligo_size -type dna";
    }
  }

print "<PRE>command: $command $parameters<P>\n</PRE>" if ($ENV{rsat_echo} >= 1);


### execute the command ###
if ($query->param('output') eq "display") {
    $sequence_file = "$TMP/$tmp_file_name.res";

    open RESULT, "$command $parameters |";

    ### open the mirror file ###
    if (open MIRROR, ">$sequence_file") {
	$mirror = 1;
	&DelayedRemoval($sequence_file);
    }

    ### print the result ### 
    &PipingWarning();
    print '<H4>Result</H4>';
    print "<PRE>";
    while (<RESULT>) {
	print "$_";
	print MIRROR $_ if ($mirror);
    }
    print "</PRE>";
    close RESULT;

    ### prepare data for piping
    &PipingFormForSequence();
    print "<HR SIZE = 3>";

    &DelayedRemoval($sequence_file);

} elsif ($query->param('output') =~ /server/i) {
    &ServerOutput("$command $parameters", $query->param('user_email'), $tmp_file_name);
} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'), $tmp_file_name);
}
print $query->end_html;

exit(0);

