#!/usr/bin/perl
#### this cgi script fills the HTML form for the program crer-scan
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

################################################################
### default values for filling the form
$default{sites} = '';
$default{in_format} = "ft";
$default{limits} = "None";
@supported_feature_types = qw(bed ft);

## Threshold fields
@threshold_fields = qw(site_score site_pval crer_size crer_sites crer_sites_distance crer_sig);

$descr{site_pval} = "P-value of input sites";
$default{lth_site_pval} = "None";
$default{uth_site_pval} = "1e-3";

$descr{site_score} = "Score of input sites";
$default{lth_site_score} = "None";
$default{uth_site_score} = "None";

$descr{crer_size} = "CRER size";
$default{lth_crer_size} = 100;
$default{uth_crer_size} = 500;

$descr{crer_sites} = "Sites per CRER";
$default{lth_crer_sites} = 2;
$default{uth_crer_sites} = "None";

$descr{crer_sites_distance} = "Inter-site distance (bp)";
$default{lth_crer_sites_distance} = 1;
$default{uth_crer_sites_distance} = 100;

$descr{crer_sig} = "CRER significance";
$default{lth_crer_sig} = 2;
$default{uth_crer_sig} = "None";


### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
} 

################################################################
### print the form ###


################################################################
### header
&RSA_header("crer-scan");
&PrintDescription();

print $query->start_multipart_form(-action=>"crer-scan.cgi");

print "<FONT FACE='Helvetica'>";


################################################################
## Sites
print "<b><a href='help.crer-scan.html#sites'>Sites</a></b>&nbsp;";

## Site format
print "&nbsp;"x10, "<B><A HREF='help.crer-scan.html#in_format'>Format</A></B>&nbsp;";
print $query->radio_group(-name=>'in_format',
			  -values=>[@supported_feature_types],
			  -default=>$default{in_format});

## Enter sites in a text area
print "<br>\n";
print $query->textarea(-name=>'sites',
		       -default=>$default{sites},
		       -rows=>6,
		       -columns=>80);

## Upload a file with the sites
print "<BR>Upload sites file from your computer<BR>\n";
print $query->filefield(-name=>'uploaded_file',
			-default=>'',
			-size=>45,
			-maxlength=>200);
print "<br>\n";


print "<hr>\n";

################################################################
## Return or not the sequence limits
print "<h2>Thresholds</h2>";


print "<blockquote>\n";
&PrintThresholdTable(@threshold_fields);
print "</blockquote>\n";

print "<h2>Output options</h2>";

print "<b><a href='help.crer-scan.html#limits'>Report sequence limits (only if provided in site file)</a></b>&nbsp;";
print $query->radio_group(-name=>'limits',
			  -values=>['None','Filtered','All'],
			  -default=>$default{limits});
print "<br>\n";

print "<hr>\n";

################################################################
## Send results by email or display on the browser
&SelectOutput();

################################################################
### action buttons
print "<UL><UL><TABLE>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

################################################################
## Data for the demo
print $query->start_multipart_form(-action=>"crer-scan_form.cgi");
$demo_file = "demo_files/Drosophila_melanogaster_eve_segmentation_sites_pval0.001.ft";
#$demo_sites = `grep -v '^;' $demo_file`;
$demo_sites = `cat $demo_file`;
print "<TD><B>";
print $query->hidden(-name=>'sites',-default=>$demo_sites);
print $query->hidden(-name=>'in_format',-default=>"ft");
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;

print "<TD><B><A HREF='help.crer-scan.html'>MANUAL</A></B></TD>\n";
#print "<TD><B><A HREF='tutorials/tut_crer-scan.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE class='formbutton'></UL></UL>\n";

print "</BLOCKQUOTE>\n";
print "</FONT>\n";

print $query->end_html;

exit(0);


################################################################
## Print a threshold table
sub PrintThresholdTable {
  my @threshold_fields = @_;
  print "<table>\n";
  print "<tr>\n";
  print "<th align='left'>", "Field","</th>";
  print "<th>", "Lower threshold", "</th>";
  print "<th>", "Upper threshold", "</th>";
  print "<th>", $default{"uth_".$field},"</th>";
  print "</tr>\n";  
  foreach my $field (@threshold_fields) {
    print "<tr>\n";
    my $field_description = $field;
    if (defined($descr{$field})) {
      $field_description = $descr{$field};
    }
    print "<td align='left'>",$field_description,"</td>";
    for my $side ("lth", "uth") {
      my $threshold_name = $side."_".$field;
      my $threshold_value;
      if (defined($default{$threshold_name})) {
	$threshold_value = $default{$threshold_name};
      } else {
	$threshold_value = "None";
      }
      print "<td align='center'>\n";
      print $query->textfield(-name=>$threshold_name,
			      -default=>$threshold_value,
			      -size=>4);
      print "</td>\n";
    }
    print "</tr>\n";  
  }
  print "</table>\n";
}

################################################################
## Print a HTML the description of the tool
sub PrintDescription {
print <<end_part_1;
<center>
  <p>Detect cis-regulatory element enriched regions (CRERs).</p>
  <p>Conception<sup>c</sup>, implementation<sup>i</sup> and testing<sup>t</sup>:
    <a href='mailto:marie.artufel\@etu.univ-amu.fr'>Marie Artufel</a>,
    <a href='mailto:lucie.khamvongsa\@gmail.com'>Lucie Khamvongsa</a>,
    <a target='_blank' href='http://www.bigre.ulb.ac.be/Users/jvanheld/'>Jacques van Helden</a><sup>cit</sup>,
</center>

<div class=\"menu_heading_closed\" onclick=\"toggleMenu(\'105\')\" id=\"heading105\">
  <font color='#0D73A7'>Information on the methods used in <i>crer-scan</i></font> 
</div>

<div id=\"menu105\" class=\"menu_collapsible\">

  <p>This tool takes as input a set of "sites" (genomic coordinates),
    and reports cis-regulatory enriched regions (<b>CRER</b>),
    i.e. genomic intervals containing a higher number of sites than
    expected by chance.</p>

  <p>The tool can be used to predict cis-regulatory modules (CRM) from
    collections of transcription factor binding sites obatined by
    various methods.
    <ol>


		<li>Predictions obtained by scanning sequences with
		  position-specific scoring matrices (for example the
		  output
		  of <a href="matrix-scan_form.cgi"><i>matrix-scan</i></a></li>
		
		<li>Peaks from ChIP-seq experiments.</li>

		<li>Annotated binding sites imported from a transcription factor
		  database.</li>

		<li>Any other data source that procudes a set of
		genomic features.</li>
    </ol>
  </p>

  <hr>	    
</div></p>
end_part_1
}
