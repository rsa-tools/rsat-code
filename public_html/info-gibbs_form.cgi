#!/usr/bin/perl
############################################################
#
# $Id: info-gibbs_form.cgi
#
# 
#
############################################################
#### this cgi script fills the HTML form for the program info-gibbs
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

### default values for filling the form
$default{sequence} = "";
$default{sequence_format} = "fasta";
$default{sequence_file} = "";
$default{upload_file} = "";
$default{length} = 12;
$default{expected} = 2.0;
$default{motifs} = 1;
$default{nrun} = 5;
$default{iter} = 1000;
$default{two_strands} = "checked";
$default{bg_order} = 3;
$default{background} = "upstream-noorf";
$default{bg_level} = "organism";
$default{organism} = "";
$default{taxon} = "Fungi";

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
} 

### print the form ###
&RSA_header("info-gibbs", "form");
print '<style><!-- textarea {height: 100px; width: 510px;}--></style>';
print "<center>\n";
print "<p>An enhanced gibbs sampler, based on a stochastic optimization of the information content of position-specific scoring matrices.</p>\n";
print "</center>\n";
print "<p><b>Reference:</b> Defrance M, van Helden J. (2009) <i>Info-gibbs</i>: a motif discovery algorithm that directly optimizes information content during sampling. Bioinformatics 25(20):2715-22.\n";
print "[<a target='_blank' href='http://www.ncbi.nlm.nih.gov/pubmed/19689955'>Pubmed 19689955</a>]\n";
print "[<a target='_blank' href='http://bioinformatics.oxfordjournals.org/content/25/20/2715.long'>Open acces</a>]\n";
print "<hr>\n";

print $query->start_multipart_form(-action=>"info-gibbs.cgi");

#print "<FONT FACE='Helvetica'>";

#### input sequence
&DisplaySequenceChoice;

### add reverse complement strand
print $query->checkbox(-name=>'two_strands',
		       -checked=>$default{two_strands},
		       -label=>'');
print "<a class='iframe' href=\"help.info-gibbs.html#strand\">search both strands</a>\n";


print "<hr width=\"550\" align=\"left\"></hr>\n";
print "<b>Options</b><br />\n";

### matrix length
print "<a class='iframe' href=\"help.info-gibbs.html#length\">Matrix length</a>\n";
print $query->textfield(-name=>'length',-id=>'length',
		  -default=>$default{length},
		  -size=>5);
print "<br />\n";

### expected number of sites per sequence
print "<a class='iframe' href=\"help.info-gibbs.html#expected\">Expected number of sites per sequence</a>\n";
print $query->textfield(-name=>'expected',-id=>'expected',
		  -default=>$default{expected},
		  -size=>5);
#print "&nbsp;&nbsp;";
print "<br />\n";

### motifs
print "<a class='iframe' href=\"help.info-gibbs.html#motifs\">Number of motifs to extract</a>\n";
print $query->popup_menu(-name=>'motifs',
			 -Values=>[1,2,3,4,5],
			 -default=>$default{motifs});
#print $query->textfield(-name=>'motifs',
#			-default=>$default{motifs},
#			-size=>5);
print "<br />\n";


### iterations

print "<a class='iframe' href=\"help.info-gibbs.html#iter\">Maximum number of iterations</a>\n";
print $query->textfield(-name=>'iter',
		  -default=>$default{iter},
		  -size=>5);
#print "&nbsp;&nbsp;";
print "<br />\n";

### nrun
print "<a class='iframe' href=\"help.info-gibbs.html#runs\">Number of runs</a>\n";
print $query->popup_menu(-name=>'nrun',
			 -Values=>[1..10],
			 -default=>$default{nrun});
#print $query->textfield(-name=>'nrun',
#		  -default=>$default{nrun},
#		  -size=>5);
print "<br />\n";


## Background model
print "<br />\n";
print "<a class='iframe' href=\"help.info-gibbs.html#background\">Background model</a>\n";
print "<br />\n";
print '<input type="radio" checked="checked" value="freq_estimate" name="bg_method"/><b>Estimated from input sequences</b><br />';
&PrintGenomeSubsetBgOptions();
print "<ul>";
print "&nbsp;&nbsp;<b>Markov background order</b> \n";
print $query->popup_menu(-name=>'bg_order',
			 -Values=>[0,1,2,3,4,5],
			 -default=>$default{bg_order});
print "</ul>";
### send results by email or display on the browser
print "<hr width=\"550\" align=\"left\"></hr>\n";
&SelectOutput;

### action buttons
print "<UL><UL><TABLE class='formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset(-id=>"reset"), "</TD>\n";
print $query->end_form;

### data for the demo
$demo_sequence = "";
open($fh, "demo_files/info-gibbs_demo_seq.fa");
while(my $row = <$fh>){
    chomp $row;
    $demo_sequence .= $row;
    $demo_sequence .= "\\n";
}

print '<script>
function setDemo(demo_sequence){
    $("#reset").trigger("click");
    $("input[name=bg_method][value=background]").prop("checked",true);
    $("#organism_bg_name").val("Saccharomyces cerevisiae");
    $("#organism_bg").val("Saccharomyces_cerevisiae");
    sequence.value = demo_sequence;
    sequence_format.value = "fasta";
    $("#length").val("20");
    expected.value = "1.0";
}
</script>';
print "<TD><B>";
print '<button type="button" onclick="setDemo('. "'$demo_sequence'" .')">DEMO</button>';
print "</B></TD>\n";

print "<TD><B><A class='iframe' HREF='help.info-gibbs.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A class='iframe' HREF='tutorials/tut_info-gibbs.html'>TUTORIAL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);


