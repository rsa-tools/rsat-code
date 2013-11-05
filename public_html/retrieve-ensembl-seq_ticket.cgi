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

### Read the CGI query
$query = new CGI;

#### update log file ####
#&UpdateLogFile();

&ListParameters() if ($ENV{rsat_echo} >= 2);


#$tmp_file_name = sprintf "retrieve-ensembl-seq.%s", &AlphaDate();
$prefix = "retrieve-ensembl-seq";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1,0); $tmp_file_name = &ShortFileName($tmp_file_path);

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
    if ($query->param('homology_selection') eq 'all') {
	$homology_type = '';
    } elsif ($query->param('homology_selection') eq 'orthologs') {
	$homology_type = 'ortholog';
    } elsif ($query->param('homology_selection') eq 'paralogs') {
	$homology_type = 'paralog';
    } else {
	$homology_type = $query->param('homology_selection');
    }
}

#### queries ####
my @gene_selection = ();
#$gene_list_file = "${TMP}/${tmp_file_name}.genes";
$gene_list_file = $tmp_file_path.".genes";
if ($query->param('uploaded_file')) {
    $upload_file = $query->param('uploaded_file');
#    if ($upload_file =~ /\.gz$/) {
#	$gene_list_file .= ".gz";
#    }
    $type = $query->uploadInfo($upload_file)->{'Content-Type'};
    open QUERY, ">$gene_list_file" || &cgiError("Cannot store gene list file in temporary directory");
    while (<$upload_file>) {
	chomp;
	print QUERY;
#	push @gene_selection, $_;
    }
    close QUERY;
    
} else {
    my $gene_selection = $query->param('gene_selection');
    $gene_selection =~ s/\r/\n/g;
    $gene_selection =~ s/\n\n/\n/g;
    if ($gene_selection =~ /\S/) {
	@gene_selection = split ("\n", $gene_selection);
	open QUERY, ">$gene_list_file" || &cgiError("Cannot store gene list file in temporary directory");;
	foreach my $row (@gene_selection) {
	    print QUERY $row, "\n";
	}
	close QUERY;
    } else {
	&cgiError("You should enter at least one gene identifier in the query box..");
    }
}

### feature type
my $feattype = '';
if ($query->param('feattype')) {
    $feattype = $query->param('feattype');
}

### sequence position
my $sequence_position = '';
if ($query->param('sequence_position')) {
    $sequence_position = $query->param('sequence_position');
}

### sequence type
my $sequence_type = '';
my $type = '';
my $first_intron = 0;
my $non_coding = 0;
my $utr = 'all';
if ($query->param('sequence_type')) {
    $sequence_type = $query->param('sequence_type');
    if ($sequence_type eq 'upstream/downstream') {
	$type = $sequence_position;
    } elsif ($sequence_type eq 'gene') {
	$feattype = 'gene';
	$type = 'feature';
    } elsif ($sequence_type eq 'intron') {
	$feattype = 'intron';
	$type = 'feature';
    } elsif ($sequence_type eq 'first intron') {
	$feattype = 'intron';
	$type = 'feature';
	$first_intron = 1;
    } elsif ($sequence_type eq 'exon') {
	$feattype = 'exon';
	$type = 'feature';
    } elsif ($sequence_type eq 'non-coding exon') {
	$feattype = 'exon';
	$type = 'feature';
	$non_coding = 1;
    } elsif ($sequence_type eq 'UTRs') {
	$feattype = 'utr';
	$type = 'feature';
    } elsif ($sequence_type eq '5prime UTR') {
	$feattype = 'utr';
	$type = 'feature';
	$utr = '5prime';
    } elsif ($sequence_type eq '3prime UTR') {
	$feattype = 'utr';
	$type = 'feature';
	$utr = '3prime';
    }
}

### all transcripts
my $all_transcripts = 1;
#if (lc($query->param('alltranscripts')) eq 'on') {
#    $all_transcripts = 1;
#}

### unique sequence
my $uniqseqs = 0;
if (lc($query->param('uniqseqs')) eq 'on') {
    $uniqseqs = 1;
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

### header organism name
my $header_org = '';
if ($query->param('header_org')) {
    $header_org = $query->param('header_org');
}

### Line width
my $lw = 60;

#print  "<PRE><B>Command :</B> $command $parameters</PRE><P>" if ($ENV{rsat_echo} >= 1);


my %args = ('output' => 'ticket',
            'organism' => $organism,
#            'query' => \@gene_selection,
            'tmp_infile' => $gene_list_file,
            'noorf' => $noorf,
            'nogene' => $nogene,
            'from' => $from,
            'to' => $to,
            'feattype' => $feattype,
            'all_transcripts' => $all_transcripts,
            'unique_sequences' => $uniqseqs,
            'first_intron' => $first_intron,
            'non_coding' => $non_coding,
            'type' => $type,
            'utr' => $utr,
            'repeat' => $rm,
            'mask_coding' => $mask_coding,
            'line_width' => $lw,
            'ortho' => $ortho,
            'taxon' => $taxon,
            'homology_type' => $homology_type,
            'header_organism' => $header_org
    );


# $ENV{rsat_echo}=2;


### print the header
#&RSA_header("retrieve-ensembl-seq result", 'results');

#### execute the command #####
if (($query->param('output') =~ /display/i) ||
    ($query->param('output') =~ /server/i)) {

    ## Service call
    my $soap=MyInterfaces::RSATWebServices::RSATWSPortType->new();

    ## Avoid timeout
    $soap->get_transport()->timeout(600); 

    ## Send the request to the server
    my $som = $soap->retrieve_ensembl_seq({'request' => \%args});
    
    ## get the ticket
    my $results = $som->get_response();

    ## Report the remote command
    my $command = $results -> get_command();
    my $ticket = $results -> get_server();

#    my $submit_time = &AlphaDate();
    
	@months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
	@weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
	($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	$year = 1900 + $yearOffset;
	if ($second <10) {
		$second = "0".$second;	
	}
	if ($minute <10) {
		$minute = "0".$minute;	
	}
	$submit_time = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";

    ### transfer to the waiting page via a hidden form### 
    print $query->header();
    print $query->start_html( -title => 'Sending form');
    my @params = $query->param();
    print '<FORM name="send2wait" id="send2wait" method="POST" action="ws_async.cgi">';
    foreach my $parameter (sort @params) {
        print ' <INPUT type="hidden" name="'.$parameter.'" value="'.$query->param($parameter).'">';  
      }
     print ' <INPUT type="hidden" name="ticket" value="'.$ticket.'">';  
     print ' <INPUT type="hidden" name="command" value="'.$command.'">';  
     print ' <INPUT type="hidden" name="title" value="retrieve-ensembl-seq">';
     print ' <INPUT type="hidden" name="submit_time" value="'.$submit_time.'">';    
     
     print '</FORM>';
      print '<script type="text/javascript" language="JavaScript">
		//submit form
		document.send2wait.submit();
	</script>';
     print $query->end_html . "\n";    
    }

exit(0);
