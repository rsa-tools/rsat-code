#!/usr/bin/env perl
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
@supported_feature_formats = qw(bed ft);

## Threshold fields
@threshold_fields = qw(site_score site_pval crer_size crer_sites crer_sites_distance crer_sig overlap);

$descr{site_pval} = "P-value of input sites";
$default{lth_site_pval} = "Disabled";
$default{uth_site_pval} = "1e-3";
$default{demo_descr} = "";
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

$descr{overlap} = "Overlap between sites";
$default{lth_overlap} = "Disabled";
$default{uth_overlap} = 1;

$descr{crer_sig} = "CRER significance";
$default{lth_crer_sig} = 1;
$default{uth_crer_sig} = "Disabled";


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

## Demo description
#print $default{demo_descr} if ($default{demo_descr});
print "<textarea id='demo' style='display:none'></textarea>";
print "<div id='demo_descr'></div>";

print $query->start_multipart_form(-action=>"crer-scan.cgi");

print "<FONT FACE='Helvetica'>";


################################################################
## Sites
print "<b><a class='iframe' href='help.crer-scan.html#sites'>Sites</a></b>&nbsp;";

## Site format
print "&nbsp;"x10, "<B><A class='iframe' HREF='help.crer-scan.html#in_format'>Format</A></B>&nbsp;";
foreach $format (@supported_feature_formats){
    print $query->radio_group(-name=>'in_format',-id=>'in_format_'.$format,
    -values=>$format,
    -default=>$default{in_format});
}

## Enter sites in a text area
print "<br>\n";
print $query->textarea(-name=>'sites',-id=>'sites',
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
&FieldsThresholdsTable("help.crer-scan.html", \@threshold_fields, "", 1);
#&PrintThresholdTable(@threshold_fields);
print "</blockquote>\n";

print "<h2>Output options</h2>";

print "<b><a class='iframe' href='help.crer-scan.html#limits'>Report sequence limits (only if provided in site file)</a></b>&nbsp;";
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
print "<TD>", $query->reset(-id=>"reset"), "</TD>\n";
print $query->end_form;

################################################################
## Data for the demo
$demo_file = "demo_files/Drosophila_melanogaster_eve_segmentation_sites_pval0.001_nocomments.ft";
$demo_sites = "";
open ($fh, $demo_file);
while($row = <$fh>){
    chomp $row;
    $demo_sites .= $row;
    $demo_sites .= "\\n";
}

print '<script>
function setDemo(demo_sites){
    $("#reset").trigger("click");
    descr = "<H4>Comment on the demonstration :</H4>\n<blockquote class =\'demo\'><p>In this demonstration, we detect cis-regulatory enriched regions (CRER) in the 5kb upstream region of the Drosophila gene even-skipped, in order to detect putative cis-regulatory modules (CRM).\nThe program <i>crer-scan</i> takes as input a set of binding sites, and returns regions presenting a significant enrichment for these.\nBinding sites were obtained by scanning the even-skipped upstream region with 12 PSSM corresponding to transcrition factors involvd in embryonic segmentation.\nFor this demo, we intently choose a lenient significance threshold, in order to increase the sensitivity of the scan. The result contains groups of mutually overlapping CRERs, which can be displayed with <i>feature-map</i>.\n</blockquote>";
    
    demo_descr.innerHTML = descr;
    demo.value = descr;
    sites.value = demo_sites;
    $("#lth_crer_sig").val("0.1");
    $("#in_format_ft").prop("checked",true);
    
}
</script>';


print "<TD><B>";
print '<button type="button" onclick="setDemo('. "'$demo_sites'".')">DEMO</button>';
print "</B></TD>\n";

print "<TD><B><A class='iframe' HREF='help.crer-scan.html'>MANUAL</A></B></TD>\n";
#print "<TD><B><A HREF='tutorials/tut_crer-scan.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:Jacques.van-Helden\@univ-amu.fr'>MAIL</A></B></TD>\n";
print "</TR></TABLE class='formbutton'></UL></UL>\n";

print "</BLOCKQUOTE>\n";
print "</FONT>\n";

print $query->end_html;

exit(0);


################################################################
## Print a threshold table
# sub PrintThresholdTable {
#   my @threshold_fields = @_;
#   print "<table>\n";
#   print "<tr>\n";
#   print "<th align='left'>", "Field","</th>";
#   print "<th>", "Lower threshold", "</th>";
#   print "<th>", "Upper threshold", "</th>";
#   print "<th>", $default{"uth_".$field},"</th>";
#   print "</tr>\n";  
#   foreach my $field (@threshold_fields) {
#     print "<tr>\n";
#     my $field_description = $field;
#     if (defined($descr{$field})) {
#       $field_description = $descr{$field};
#     }
#     print "<td align='left'>",$field_description,"</td>";
#     for my $side ("lth", "uth") {
#       my $threshold_name = $side."_".$field;
#       my $threshold_value;
#       if (defined($default{$threshold_name})) {
# 	$threshold_value = $default{$threshold_name};
#       } else {
# 	$threshold_value = "None";
#       }
#       print "<td align='center'>\n";
#       print $query->textfield(-name=>$threshold_name,
# 			      -default=>$threshold_value,
# 			      -size=>4);
#       print "</td>\n";
#     }
#     print "</tr>\n";  
#   }
#   print "</table>\n";
# }

################################################################
## Print a HTML the description of the tool
sub PrintDescription {
print <<end_part_1;
<center>
  <p>Detect cis-regulatory element enriched regions (CRERs).</p>
  <p>Conception<sup>c</sup>, implementation<sup>i</sup> and testing<sup>t</sup>:
    <a href='mailto:marie.artufel\@etu.univ-amu.fr'>Marie Artufel</a>,
    <a href='mailto:lucie.khamvongsa\@gmail.com'>Lucie Khamvongsa</a>,
    <a target='_blank' href='http://jacques.van-helden.perso.luminy.univ-amu.fr/'>Jacques van Helden</a><sup>cit</sup>,
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

		<li>Any other data source that produces a set of
		genomic features.</li>
    </ol>
  </p>

  <hr>	    
</div></p>
end_part_1
}
