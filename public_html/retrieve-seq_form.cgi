#!/usr/bin/perl
#### this cgi script fills the HTML form for the program upstream-region
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib.pl";
require "RSA.cgi.lib.pl";

### Read the CGI query
$query = new CGI;

### default values for filling the form
$default{seq_format} = "fasta";
$default{organism} = "Saccharomyces cerevisiae";
$default{from} = "-800";
$default{to} = "-1";
$default{genes} = "";

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
} 

### print the form ###
&RSA_header("upstream-region");

### head
print "<CENTER>";
print "Returns the upstream sequences for a list of genes<P>\n";
print "</CENTER>";

print $query->start_multipart_form(-action=>"upstream-region.cgi");

print "<FONT FACE='Helvetica'>";

&OrganismPopUp;

### query (gene list)
print "<B><A HREF='help.upstream-region.html#genes'>Genes</A></B>&nbsp;";
print $query->radio_group(-name=>'genes',
			  -values=>['all','selection'],
			  -default=>'selection');

print "<BR>\n";
print "<UL>\n";

print $query->textarea(-name=>'gene_selection',
		       -default=>$default{genes},
		       -rows=>6,
		       -columns=>40);
print "</UL>\n";
print "<BR>\n";

### from to

print "<B><A HREF='help.upstream-region.html#from_to'>From</A></B>&nbsp;\n";
print $query->textfield(-name=>'from',
			-default=>$default{from},
			-size=>10);

print "&nbsp;&nbsp;";
print "<B><A HREF='help.upstream-region.html#from_to'>To</A></B>&nbsp;\n";
print $query->textfield(-name=>'to',
			-default=>$default{to},
			-size=>10);
print "<BR>\n";

### allow ORF overlap
### temporarily inactivated because it does not work with all organisms
#  print $query->checkbox(-name=>'orf_overlap',
#  		       -checked=>'checked',
#  		       -label=>'');
#  print "&nbsp;<A HREF='help.upstream-region.html#noorf'><B>allow overlap with upstream ORFs</B></A>";
#  print "<BR>\n";

print $query->hidden(-name=>'orf_overlap',-default=>'on');


### sequence format 
print "<B><A HREF='help.upstream-region.html#formats'>Sequence format</A></B>&nbsp;";
print $query->popup_menu(-name=>'format',
			 -Values=>['fasta', 
				   'IG',
				   'wconsensus',
				   'raw',
				   'multi'],
			 -default=>$default{seq_format});
print "<BR>\n";

### send results by e-mail or display on the browser
&SelectOutput;

### action buttons
print "<UL><UL><TABLE>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

### data for the demo 
print $query->start_multipart_form(-action=>"upstream-region_form.cgi");
$demo_genes = "PHO5\n";
$demo_genes .= "PHO8\n";
$demo_genes .= "PHO11\n";
$demo_genes .= "PHO81\n";
$demo_genes .= "PHO84\n";
print "<TD><B>";
print $query->hidden(-name=>'genes',-default=>$demo_genes);
print $query->hidden(-name=>'organism',-default=>"Saccharomyces cerevisiae");
print $query->hidden(-name=>'from',-default=>"-800");
print $query->hidden(-name=>'to',-default=>"-1");
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


#print "<TD><B><A HREF='demo.upstream-region.html'>DEMO</A></B></TD>\n";
print "<TD><B><A HREF='help.upstream-region.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@ucmb.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);

