#!/usr/bin/perl
# the first line of the script must tell us which language interpreter to use,
# in this case its perl

use strict;

# import the modules we need for this test; XML::Compile is included on the server
# by default.
use XML::Compile::SOAP11;
use XML::Compile::WSDL11;
use XML::Compile::Transport::SOAPHTTP;

eval
{
    # Retriving and processing the WSDL
    my $wsdl  = XML::LibXML->new->parse_file('http://rsat.bigre.ulb.ac.be/rsat/web_services/RSATWS.wsdl');
    my $proxy = XML::Compile::WSDL11->new($wsdl);

    # Generating a request message based on the WSDL
    my $getSequence = $proxy->compileClient('peak_motifs');

    # Defining a few parameters
    my $test = `cat peak_motifs_test1_demo_seq.fasta`;

my $verbosity = 1;
#my $format = 'fasta';  ## The format of input sequences
my $min_length = 6;  ## Length of patterns to be discovered
my $max_length = 7;  ## Length of patterns to be discovered
my $top_peaks = 5000;
my $max_seq_length = 500;
my $noov = 1;  ## Do not allow overlapping patterns
my $str = 2;  ## Search on both strands
my $image_format = 'png';
my $task = 'purge,seqlen,composition,merge_words,collect_motifs,timelog,synthesis,oligos,positions';
my $graph_title = 'Oct4_Chen2008';

my %args = (
	    'verbosity' => $verbosity,
            'test' => $test,
#            'format' => $format,
            'min_length' => $min_length,
            'max_length' => $max_length,
	    'top_peaks' => $top_peaks,
	    'max_seq_length' => $max_seq_length,
	    'graph_title' => $graph_title,
            'noov' => $noov,
            'str' => $str,
	    'image_format' => $image_format,
	    'task' => $task
           );

    # Calling the service and getting the response
    my $answer = $getSequence->( request => {%args});
#    print "Answer: ".$answer."\n";
    # If the response arrived, look for a specific pattern
    # If the pattern is present, return 0 because the test passed.
    # If the result is something else, return 2 to indicate a warning.
    # If no answer has arrived, return 1 to indicate the test failed.
    if ( defined $answer ) {
	print "Result : ".$answer->{output}->{response}->{client}."\n";
	print "Command : ".$answer->{output}->{response}->{command}."\n";
	print "Server : ".$answer->{output}->{response}->{server}."\n";
# 	if ($answer->{output}->{response}->{client} =~ 'tgccaa'){
# 	    print "Passed\n";
# 	    exit 0;
# 	} else {
# 	    print "Unexpected data\n";
# 	    exit 2;
# 	}
    } else {
	print "Failed\n";
	exit 1;
    }
};

if ($@)
{
    print "Caught an exception\n";
    print $@."\n";
    exit 1;
}
