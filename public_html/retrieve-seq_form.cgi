#!/usr/bin/perl
#### this cgi script fills the HTML form for the program retrieve-seq
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

### default values for filling the form
$default{sequence_format} = "fasta";
$default{seq_label} = "gene name";
$default{organism} = "Saccharomyces cerevisiae";
$default{noorf} = "checked";
$default{from} = "default";
$default{to} = "default";
$default{genes} = "selection";
$default{gene_selection} = "";
$default{sequence_type} = "upstream";

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
} 

### print the form ###
&RSA_header("retrieve sequence");

### head
print "<CENTER>";
print "Returns upstream, downstream or ORF sequences for a list of genes<P>\n";
print "</CENTER>";


print $query->start_multipart_form(-action=>"retrieve-seq.cgi");

print "<FONT FACE='Helvetica'>";

&OrganismPopUp;

### query (gene list)
print "<B><A HREF='help.retrieve-seq.html#genes'>Genes</A></B>&nbsp;";
print $query->radio_group(-name=>'genes',
			  -values=>['all','selection'],
			  -default=>$default{genes});

print "<BR>\n";
print "<UL>\n";

print $query->textarea(-name=>'gene_selection',
		       -default=>$default{gene_selection},
		       -rows=>6,
		       -columns=>40);

### option to upload a file with the gene list from the client machine 
print "<BR>Upload gene list from file<BR>\n";
print $query->filefield(-name=>'uploaded_file',
			-default=>'',
			-size=>45,
			-maxlength=>200);

print "</UL>\n";
print "<BR>\n";

### sequence type
print "<B><A HREF='help.retrieve-seq.html#sequence_type'>Sequence type</A></B>&nbsp;";
print $query->popup_menu(-name=>'sequence_type',
			 -Values=>['upstream','downstream','ORFs (unspliced)'],
			 -default=>$default{sequence_type});

### from to

print "<B><A HREF='help.retrieve-seq.html#from_to'>From</A></B>&nbsp;\n";
print $query->textfield(-name=>'from',
			-default=>$default{from},
			-size=>5);

print "&nbsp;&nbsp;";
print "<B><A HREF='help.retrieve-seq.html#from_to'>To</A></B>&nbsp;\n";
print $query->textfield(-name=>'to',
			-default=>$default{to},
			-size=>5);
print "<BR>\n";

### prevent ORF overlap
print $query->checkbox(-name=>'noorf',
  		       -checked=>$default{noorf},
  		       -label=>'');
print "&nbsp;<A HREF='help.retrieve-seq.html#noorf'><B>Prevent overlap with upstream ORFs</B></A>";
print "<BR>\n";


### sequence format 
print "<B><A HREF='help.retrieve-seq.html#formats'>Sequence format</A></B>&nbsp;";
print $query->popup_menu(-name=>'format',
			 -Values=>['fasta', 
				   'IG',
				   'wconsensus',
				   'multi'],
			 -default=>$default{sequence_format});
print "<BR>\n";

### sequence label
print "<B><A HREF='help.retrieve-seq.html#seq_label'>Sequence label</A></B>&nbsp;";
print $query->popup_menu(-name=>'seq_label',
			 -Values=>['ORF identifier', 
				   'gene name',
				   'ORF id + gene name',
				   'full identifier'
				   ],
			 -default=>$default{seq_label});
print "<BR>\n";


### send results by e-mail or display on the browser
&SelectOutput();

### data for the demo 
#@demo_genes = qw ( MET1 MET2 MET3 MET6 MET14 MET19 MET25 MET30 MUP3
#		    SAM1 SAM2);
@demo_genes = qw (DAL5 GAP1 MEP1 MEP2 PUT4 MEP3 DAL80);
$demo_genes = join "\n", @demo_genes;


### action buttons
print "<UL><UL><TABLE>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

print $query->start_multipart_form(-action=>"retrieve-seq_form.cgi");
print "<TD><B>";
print $query->hidden(-name=>'gene_selection',-default=>$demo_genes);
print $query->hidden(-name=>'organism',-default=>"Saccharomyces cerevisiae");
print $query->hidden(-name=>'from',-default=>"-800");
print $query->hidden(-name=>'to',-default=>"-1");
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


#print "<TD><B><A HREF='demo.retrieve-seq.html'>DEMO</A></B></TD>\n";
print "<TD><B><A HREF='help.retrieve-seq.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='tutorials/tut_retrieve-seq.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@ucmb.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);

