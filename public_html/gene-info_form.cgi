#!/usr/bin/env perl
#### this cgi script fills the HTML form for the program gene-info
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
$default{organism} = "";
$default{queries} = '';
$default{full} = '';
$default{match_description} = '';
$default{feattype} = "gene";

## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## TEMPORARY (2015-09): RESTRICT SUPPORTED FEATURE TYPES until the switch from NCBI
## to EnsemblGenomes as genome source is completely checked.
@supported_feature_types = qw(gene mRNA CDS);
##
## END TEMPORARY
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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
&RSA_header("gene-info");
print "<center>";
print "<p>Get information about genes (CDS, mRNA, ...) specified by identifiers, names or  descriptions. Regular expressions are supported.</p>\n";
print "</center>";
print "<blockquote>\n";

print $query->start_multipart_form(-action=>"gene-info.cgi", -id=>"form");

print "<FONT FACE='Helvetica'>";

################################################################
## Choice of the organism
&OrganismPopUp();

################################################################
## Queries
print "<B><A class='iframe' HREF='help.gene-info.html#queries'>Gene queries</A></B>&nbsp;";
print "<BR>\n";
print $query->textarea(-id=>'queries',
              -name=>'queries',
		       -default=>$default{queries},
		       -rows=>6,
		       -columns=>40);


### option to upload a file with the gene list from the client machine 
print "<BR>Upload gene list from file<BR>\n";
print $query->filefield(-name=>'uploaded_file',
			-default=>'',
			-size=>45,
			-maxlength=>200);
print "<BR>\n";

################################################################
## Feature type
print "<B><A class='iframe' HREF='help.retrieve-seq.html#feattype'>Feature type</A></B>&nbsp;";
print $query->radio_group(-name=>'feattype',
			  -values=>[@supported_feature_types],
			  -default=>$default{feattype});
print "<BR>\n";


################################################################
## Full match
print $query->checkbox(-name=>'full',
  		       -checked=>$default{full},
  		       -label=>'');
print "&nbsp;<A class='iframe' HREF='help.retrieve-seq.html#full'><B>Full string matching</B></A>";
print "<P>\n";

################################################################
## Match queries against description
print $query->checkbox(-name=>'match_description',
  		       -checked=>$default{match_description},
  		       -label=>'');
print "&nbsp;<A class='iframe' HREF='help.retrieve-seq.html#match_description'><B>Match queries against description</B></A>";
print "<P>\n";

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
$demo_queries = "ARG3\n";
$demo_queries .= "PHO\n";
$demo_queries .= "YBL05[\\d]W\n";
print "<TD><B>";
print '<script>
function setDemo(){
    $("#reset").trigger("click");
    queries.value = "ARG3\nPHO\nYBL05[\\\d]W\n";
    $("#organism").val("Saccharomyces_cerevisiae");
    $("#organism_name").val("Saccharomyces cerevisiae");
}
</script>';

print '<button type="button" onclick="setDemo();">DEMO</button>';

print "</B></TD>\n";


print "<TD><B><A class='iframe' HREF='help.gene-info.html'>MANUAL</A></B></TD>\n";
#print "<TD><B><A HREF='tutorials/tut_gene-info.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:Jacques.van-Helden\@univ-amu.fr'>MAIL</A></B></TD>\n";
print "</TR></TABLE class='formbutton'></UL></UL>\n";

print "</BLOCKQUOTE>\n";
print "</FONT>\n";

print $query->end_html;

exit(0);

