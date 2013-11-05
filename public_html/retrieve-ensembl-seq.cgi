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

### print the header
&RSA_header("retrieve-ensembl-seq result", 'results');
#print $query->header();

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

# $ENV{rsat_echo}=2;
&ListParameters() if ($ENV{rsat_echo} >= 2);

$prefix = "retrieve-ensembl-seq";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); $tmp_file_name = &ShortFileName($tmp_file_path);
#$tmp_file_name = sprintf "retrieve-ensembl-seq.%s", &AlphaDate();
#@result_files = ();

################################################################
## Single or multi-genome query

#### organism
my $organism_name = $query->param('organism');

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
    s/\r//g;
    print QUERY $_, "\n";
#    print QUERY "\n";
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
#&RSAT::message::Debug("Gene list", $gene_list_file, "\n", `cat $gene_list_file`);


## The temp file is on the local server, but the queries have to be passed to the RSATWS server
## We thusxtract the queries from the temp file and send them as an array
my @gene_queries = ();
my ($gene_fh) = &OpenInputFile($gene_list_file);
while (<$gene_fh>) {
    if (/(\S+)/) {
	push @gene_queries, $1;
    }
}
#&RSAT::message::Debug("Queries", @gene_queries);

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
    } elsif ($sequence_type eq 'transcript') {
	$feattype = 'transcript';
	$type = 'feature';
    } elsif ($sequence_type eq 'mRNA') {
	$feattype = 'mrna';
	$type = 'feature';
    } elsif ($sequence_type eq 'CDS') {
	$feattype = 'cds';
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

my %args = (
            'organism' => $organism_name,
#            'query' => \@gene_selection,
#            'tmp_infile' => $gene_list_file,
            'query' => \@gene_queries,
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


## Execute the command
if (($query->param('output') =~ /display/i) ||
    ($query->param('output') =~ /server/i)) {

    $args{'output'} = 'ticket';
#    $args{'output'} = 'both';

    my ($ticket, $command, $results) = &Retrieve(%args);

    &ReportWebCommand($command) if ($ENV{rsat_echo} >= 1);

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
    my $submit_time = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";

#    &RSAT::message::Debug("ticket", $ticket);
#    &RSAT::message::Debug("command", $command);

    ### Print the result page
#    my $css_body_class = "form";
#    my $title = $query->param('title');
##    print $query->header();
##    print $query->start_html(
##	-title=>"RSA-tools : $title",
#	-class => "$css_body_class",
##	-author=>'jvanheld@bigre.ulb.ac.be',
#	-style => { 	-src => "$ENV{rsat_www}/main.css",
#			-type => 'text/css',
#			-media => 'screen' }
##	);

    ## Transfer to the result page via a hidden form
    print '<FORM name="send2result" id="send2result" method="POST" action="ws_async.cgi">';
    print ' <INPUT type="hidden" name="output" value="'.$query->param('output').'">';
    print ' <INPUT type="hidden" name="ticket" value="'.$ticket.'">';
    print ' <INPUT type="hidden" name="command" value="'.$command.'">';
    print ' <INPUT type="hidden" name="title" value="retrieve-ensembl-seq">';
    print ' <INPUT type="hidden" name="submit_time" value="'.$submit_time.'">';
    print ' <INPUT type="hidden" name="organism_name" value="'.$organism_name.'">';
    print '</FORM>';
    print '<script type="text/javascript" language="JavaScript">
		//submit form
		document.send2result.submit();
	</script>';

} elsif ($query->param('output') =~ /email/i) {
    ### Print the notification page
#    my $css_body_class = "form";
#    my $title = $query->param('title');
#    print $query->header();
#    print $query->start_html(
#	-title=>"RSA-tools : $title",
#	-class => "$css_body_class",
#	-author=>'jvanheld@bigre.ulb.ac.be',
#	-style => { 	-src => "$ENV{rsat_www}/main.css",
#			-type => 'text/css',
#			-media => 'screen' }
#	);
    $args{'output'} = $query->param('user_email');
    my $notification = &Retrieve(%args);
    $notification =~ s|(http://\S+)|<a href=$1>$1</a>|gm;
    &Info($notification);
}

print $query->end_html . "\n";

exit(0);


sub Retrieve {
    my %args = @_;

    ## Service call
    my $soap=MyInterfaces::RSATWebServices::RSATWSPortType->new();

    ## Send the request to the server
    my $som = $soap->retrieve_ensembl_seq({'request' => \%args});

    ## get the ticket
    my $results = $som->get_response();

    if ($args{'output'} eq 'ticket') {
	## Report the remote command
	my $command = $results -> get_command();
	my $ticket = $results -> get_server();
	return ($ticket, $command, $results);
    } else {
	my $notification = $results -> get_client();
	return $notification;
    }
}

# sub PrintPage {
#     my ($title,$ticket,$submit_time) = @_;

#     print $query->h2('Processing your ' . $title .' request'). "\n";
#     print "<hr/>";
#     print '<TABLE>' . "\n";
#     print "<tr><th>Job ID</th><td>" . $ticket . "</td></tr>\n";
#     print "<tr><th>Status</th><td>Running</td><td><img src='images/loader.gif'/></tr>\n";
#     print "<tr><th>Submitted at </th><td>".$submit_time."</td></tr>\n";
#     print "</TABLE>\n"; 
#     print "<hr/>";
#     print "<h2>The RSAT servers are working for you, take a break with the latest strips from 
#      <a href='http://www.phdcomics.com/' target=_blank> PHD Comics </a> !</h2>";
#     print '<iframe width="800" height="400" src="http://www.rss-info.com/rss2.php?integration=if&windowopen=1&rss=http%3A%2F%2Fwww.phdcomics.com%2Fgradfeed_justcomics.php&number=10&width=800&ifbgcol=FFFFFF&bordercol=D0D0D0&textbgcol=F0F0F0&rssbgcol=F0F0F0&showrsstitle=1&showtext=1" frameborder=0></iframe>';
# }

# sub JobStatus {
#     my $ticket = shift;
#     my $soap=MyInterfaces::RSATWebServices::RSATWSPortType->new();

#     my %args = (
# 	'ticket' => $ticket
# 	);

#     my $som = $soap->monitor({'request' => \%args});

# #    print "<h1>SOM= $som</h1>";

#     unless ($som) {
# 	printf "A fault (%s) occured: %s\n", $som->get_faultcode(), $som->get_faultstring();
#     } else {
# 	my $response = $som->get_response();

# 	## Report the status
# 	my $status = $response -> get_status();

# 	return $status;
#     }
# }

