#!/usr/bin/perl
# the first line of the script must tell us which language interpreter to use,
# in this case its perl

use strict;

# import the modules we need for this test; XML::Compile is included on the server
# by default.
use XML::Compile::SOAP11;
use XML::Compile::WSDL11;
use XML::Compile::Transport::SOAPHTTP;
use Data::Dumper;

eval
{
    # Retriving and processing the WSDL
    my $wsdl  = XML::LibXML->new->parse_file('http://rsat.ulb.ac.be/rsat/web_services/RSATWS.wsdl');
    my $proxy = XML::Compile::WSDL11->new($wsdl);
    
    # Generating a request message based on the WSDL
    my $getSequence = $proxy->compileClient('matrix_scan');
    
    #Defining a few parameters
    my $sequence = '>Hoxb1  ENSMUSG00000018973; upstream from -2000 to -1; size: 2000; location: chromosome:NCBIM36:11:1:121798632:1 96179917 96181916 D; upstream neighbour: ENSMUSP00000098092 (distance: 12200)
GTGCCTCTGAGGATGCCGCATGGAAAGAGGAGGGGCAACTAGAGAGCCAGACAAGGATGT
GGGTGGGGAGAGGGGGGAGGGGAGCCAGCAAGGTGAAGAGGCCGGGGAGATAGAGACAAC
TTTGGGCCCTTGAAGGGAGTGGAATGCATGCAAGCCTGGCTTTTTTTTTTTTTTTTTTTT
TTTAGGGTAAGGAAACCTGAGAAGTTTCAAGGGTTACCCCAGACCCGTCCCACCTGCTTT
GGCTTCCTCCGCAGAGTGCCACTGTTTACGGAGATCCCTCCCTGAACTCTTGCCCTCCTG
GACTTGCCCTAGCTCAGGCCCCAGGCCTGTGGCCAGGCAGACACCCTGACAAGTTACAAA
TGAGAGTGGGTGTTGGATTCTTGTCTTCAGAGTCTGGAGGAGGAGACATCAATGAGCTCT
ACTACCCAAGAGCATCTCTTCTAATTCCAAACTGCCTGCTGCATTTCAGAGAGTGAGCAA
GGCTCTTCCTGTTTTCCCTCCCTTGATCTTAACCCAAGAGAAAGGAAGAAAGAAAGAAAG
AGAGAGAGAGAGAGAGAGAGTGGGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAAAG
GAAGGAAGAGAGAGAGAGAGAGAGATCCATTTCTAAGCCTTAGTTTAGAAGCTCTACACC
AATTCTCCCATTCTCTGAGCAAACTAAACTCCCTTTATCTATTCTAGAACATGCATGGAA
ATTTTGGATAGTCTTTGAGACTTAACTAGCCTTGGATTTATGTATTTGAGGCTTGGAACT
TAAAAAAAGAAACTTTTAGGGAAACTAAAGATAATGGAGGGGTGTCAGGGAGGAGACCAG
AGAAGTGCAGAGGATCTTTACTAGTTTTGAACACTTTTGAAAAGAGAATAAAAATTCGCT
TCCTTCTCACCCATCTGTAAATATCACCTGGCTCTCTTCATTGCTGCTCTCTGGAGAGAT
ATAGGGAGTTATTTAAGTCTGAGTGAGTGAGCCAGCTACAGCTACGGATGGCTAGAGAAA
GGATAAATACATCTAGAGAGATGAGTTGATGGATAGGTAGGGAGCGAATATAAATATGGA
GAGAGATGGACGAACACAGAGATGAATGGAGACAGAGGGAATATAAAGAGAGGCTGAGGG
AGAGAAATGGCCTCGCGAATGAGGAGAGAGCGTTTATGCAGAGAGATGGAGAGAGGCGTC
CCTCATCTTCCCTCCACCTTTCCATTTCCCTTCTCCTGTTCCCTTGAAATGCCATCGTTT
TCCCTCCTCTTGATTTGTCATGAATTCTTCTGGATTTTAGACCCCTGGCTCACACCCTGA
GGTGTTTTTACAGTCTCCCTCTCCCTTCCTAGTCATCCTTTTGTCCCAAGATAACAGACC
TAACCAGGCTGTCCTAGAACACCATCCCAAGACAAACTGGGAGAAAGAAACATGGAATGG
GGGTTCTGGGATAAGTAAGGAGTCTGGCTCTCCAGTCTCTAGCCCTACAGCCTTGGGGTG
GGGGTAGGCTTCTTTGAGGGGAGAAGAGAGTGAGGGGCAGCTGGGTTGAATTTGGCCAAA
TCTAATAATCCAAGAACCTATTGAAGGCCTTGGGGGGTTGGGAGGGGAGTAAAAGTCTTG
AATATTCTTCCAACTTCTCCCCCCCCTCCCTCTGGTCCCTTCTTTCCAAAAAGTCTTTGA
AGAAAGATGTTTTTGACGCTTCCATGTCGCTCTCAGATGGATGGGCTCAGAGTGATTGAA
GTGTCTTTGTCATGCTAATGATTGGGGGGTGATGGATGGGCGCTGGGACTGCCAAACTCT
GGCCCGCTTAGCCCATTGGCCTGGGAGAGATCACATGTGCCCCCCCACCCCCACTCCCTA
GCCCCTTCCTAGGGGATCGCTGGCGGGGCCAAGCTGGCCCGGGCCATGGGCTCAAGCTTC
AGCTCTGTGACATACTGCCGAAAGGTTGTAGGGCAAGAGGGTGTCTCCCCCAAACGGCCC
GACCCTCCTTCGGCCTCTAC';

    my $matrix = 'a	2	0	13	0	2	2	13	1	0	1	0
c	2	1	0	0	0	0	0	0	0	0	1
g	8	11	0	13	3	0	0	7	0	9	8
t	1	1	0	0	8	11	0	5	13	3	4
//
a	0	2	2	1	2	0	0	13	0	1	13
c	1	0	2	0	0	1	0	0	0	0	0
g	11	3	8	7	0	8	0	0	13	9	0
t	1	8	1	5	11	4	13	0	0	3	0';

    my $background = 'upstream';
    my $organism = 'Mus_musculus_EnsEMBL';
    #my $tmp_background_infile="/home/rsat/rsat/public_html/data/genomes/Homo_sapiens_ensembl_74_GRCh37/oligo-frequencies/1nt_upstream-noorf_Homo_sapiens_ensembl_74_GRCh37-ovlp-1str.freq.gz";
    my $markov = 0;
    my $background_pseudo = 0.01;
    my @uth = ('pval 0.001');
    my $str = 2;
    my $origin = 'start';
    my $pseudo = 1;
    my $n_treatment = 'score';
#   my $quick = '1';

    my %args = (
#	'verbosity' => 1,
	'sequence' => $sequence, 
	'matrix' => $matrix,
	'background' => $background,
	'organism' => $organism,
#	'quick' => $quick,
	'markov' => $markov,
	#'tmp_background_infile' => $tmp_background_infile,
	'background_pseudo' => $background_pseudo,
	'uth' => \@uth,
	'str' => $str,
	'origin' => $origin,
	'pseudo' => $pseudo,
	'n_treatment' => $n_treatment
	);

    # Calling the service and getting the response
    my $answer = $getSequence->( request => {%args});

    # If the response arrived, look for a specific match
    # If the match is correct, return 0 because the test passed.
    # If the result is something else, return 2 to indicate a warning.
    # If no answer has arrived, return 1 to indicate the test failed.
    if ( defined $answer ) {
    warn ("Server command : ".$answer->{output}->{response}->{command}."\n");
#	print "Result : ".$answer->{output}->{response}->{client}."\n";
	if ($answer->{output}->{response}->{client} =~ 'TGGGTGTTGGA'){
	    print "Passed\n";
	    print "\nResult :\n\n", $answer->{output}->{response}->{client}, "\n";
	    exit 0;
	} else {
	    print "Unexpected data\n";
	    print "\nResult : ", Dumper($answer), "\n";
	    exit 2;
	}
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
