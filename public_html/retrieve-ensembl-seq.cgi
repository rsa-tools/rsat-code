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

use SOAP::WSDL; ## Requires version 2.0 or later of SOAP::WSDL
use lib '../ws_clients/perl_clients/RSATWS';
use MyInterfaces::RSATWebServices::RSATWSPortType;

use File::Basename;

$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

$tmp_file_name = sprintf "retrieve-ensembl-seq.%s", &AlphaDate();

### Read the CGI query
$query = new CGI;


$ENV{rsat_echo}=2;


### print the header
&RSA_header("retrieve-ensembl-seq result", 'results');


#### update log file ####
&UpdateLogFile();

&ListParameters() if ($ENV{rsat_echo} >= 2);

################################################################
## Single or multi-genome query

#### organism
my $organism = $query->param('organism');

#### Ortho
my $ortho = 0;
if ($query->param('single_multi_org') eq 'multi') {
    $ortho = 1;
}

### taxon
my $taxon = '';
if ($query->param('taxon_selection')) {
    $taxon = $query->param('taxon_selection');
}

### homology type
my $homology_type = '';
if ($query->param('homology_selection')) {
    $homology_type = $query->param('homology_selection');
}

#### queries ####
my @gene_selection = ();
if ($query->param('uploaded_file')) {
    $upload_file = $query->param('uploaded_file');
#    $gene_list_file = "${TMP}/${tmp_file_name}.genes";
#    if ($upload_file =~ /\.gz$/) {
#	$gene_list_file .= ".gz";
#    }
    $type = $query->uploadInfo($upload_file)->{'Content-Type'};
#    open SEQ, ">$gene_list_file" ||
#	&cgiError("Cannot store gene list file in temporary directory");
    while (<$upload_file>) {
#	print SEQ;
	chomp;
	push @gene_selection, $_;
    }
#    close SEQ;
    
} else {
    my $gene_selection = $query->param('gene_selection');
    $gene_selection =~ s/\r/\n/g;
    $gene_selection =~ s/\n\n/\n/g;
    if ($gene_selection =~ /\S/) {
	@gene_selection = split ("\n", $gene_selection);
    } else {
	&cgiError("You should enter at least one gene identifier in the query box..");
    }
}

### feature type
my $feattype = '';
my $first_intron = 0;
my $non_coding = 0;
if ($query->param('feattype')) {
    $feattype = $query->param('feattype');
    if ($feattype eq 'first intron') {
	$feattype = 'intron';
	$first_intron = 1;
    } elsif ($feattype eq 'non-coding exon') {
	$feattype = 'exon';
	$non_coding = 1;
    }
}

### sequence type
my $sequence_type = '';
if ($query->param('sequence_type')) {
    $sequence_type = $query->param('sequence_type');
}

### all transcripts
my $all_transcripts = 0;
if (lc($query->param('alltranscripts')) eq 'on') {
    $all_transcripts = 1;
}

### limits ###
my $from = '';
if (&IsInteger($query->param('from'))) {
    $from = $query->param('from');
}
my $to = '';
if (&IsInteger($query->param('to'))) {
    $to = $query->param('to');
}

### prevent orf/gene overlap ###
my $noorf = 0;
my $nogene = 0;
if ($query->param('prevent_overlap')) {
    if ($query->param('prevent_overlap') eq 'ORF') {
	$noorf = 1;
    } elsif ($query->param('prevent_overlap') eq 'gene') {
	$nogene = 1;
    }
}

### repeat masking
my $rm = 0;
if (lc($query->param('rm')) eq "on") {
    $rm = 1;
}

### mask coding
my $mask_coding = 0;
if (lc($query->param('maskcoding')) eq "on") {
    $mask_coding = 1;
}

#print  "<PRE><B>Command :</B> $command $parameters</PRE><P>" if ($ENV{rsat_echo} >= 1);


my %args = (
            'organism' => $organism,
            'query' => \@gene_selection,
            'noorf' => $noorf,
            'nogene' => $nogene,
            'from' => $from,
            'to' => $to,
            'feattype' => $feattype,
            'all_transcripts' => $all_transcripts,
            'first_intron' => $first_intron,
            'non_coding' => $non_coding,
            'type' => $sequence_type,
            'repeat' => $rm,
            'mask_coding' => $mask_coding,
            'ortho' => $ortho,
            'taxon' => $taxon,
            'homology_type' => $homology_type
    );

#### execute the command #####
if (($query->param('output') =~ /display/i) ||
    ($query->param('output') =~ /server/i)) {

    ## Service call
    my $soap=MyInterfaces::RSATWebServices::RSATWSPortType->new();

    ## Send the request to the server
    my $som = $soap->retrieve_ensembl_seq({'request' => \%args});

    ### print the result ### 
    &PipingWarning();

    print '<H4>Result</H4>';

    ## Get the result
    unless ($som) {
	print "A fault ".$som->get_faultcode()." occured: ". $som->get_faultstring()."\n";
    } else {
	my $results = $som->get_response();

	## Report the remote command
	my $command = $results -> get_command();
	print "Command used on the server: ".$command, "\n";
	## Report the result
	$server_file = $results -> get_server();
	$result = $results -> get_client();
    }

    print "<PRE>";
    print $result;
    print "</PRE>";

    if ($query->param('output') =~ /server/i) {
	my $server_file_name = basename($server_file);
	$result_URL = "$ENV{rsat_www}/tmp/$server_file_name";
	print ("The result is available at the following URL: ", "\n<br>",
	       "<a href=${result_URL}>${result_URL}</a>",
	       "<p>\n");
    }

    ### prepare data for piping
    &PipingFormForSequence();

    print "<HR SIZE = 3>";

}
#} elsif 
#    &ServerOutput("$command $parameters", $query->param('user_email'));
#} else {
#    &EmailTheResult("$command $parameters", $query->param('user_email'));
#}

print $query->end_html;

exit(0);
