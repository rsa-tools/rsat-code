#!/usr/bin/perl
################################################################
## this cgi script fills the HTML form for the program get-orthologs

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

################################################################
## Initialize parameters

## Supported taxons (TEMPORARY, FIXED LIST)
@supported_taxons = qw(Bacteria Proteobacteria Gammaproteobacteria Enterobacteriales Enterobacteriaceae Escherichia);


## Output fields
my @output_fields = qw(ident
			   ali_len
			   mismat
			   gap_open
			   e_value
			   bit_sc
			   q_rank
			   s_rank);
my %field_description = ();
$field_description{ident} = "Percentage of identity";
$field_description{ali_len} = "Alignment length";
$field_description{mismat} = "Number of mismatches";
$field_description{gap_open} = "Number of gap openings";
$field_description{e_value} = "E-value";
$field_description{bit_sc} = "Bit score";
$field_description{q_rank} = "query rank";
$field_description{s_rank} = "target rank";


################################################################
### default values for filling the form
$default{organism} = "Escherichia_coli_K12";
$default{taxon} = "Gammaproteobacteria";
$default{queries} = '';
$default{full} = '';
$default{match_description} = '';
$default{feattype} = "CDS";

$default{lth_ali_len} = 50;
$default{uth_e_value} = 1e-5;
$default{uth_q_rank} = 1;
$default{uth_s_rank} = 1;

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
&RSA_header("get-orthologs (VERSION IN CONSTRUCTION)");
print "<CENTER>";
print "Given a list of genes from a query organism and a taxon of interest, <br>return genes coding for similar proteins in each genome of the taxon.<P>\n";
print "</CENTER>";
print "<BLOCKQUOTE>\n";

print $query->start_multipart_form(-action=>"get-orthologs.cgi");

################################################################
#### choice of the organism
&OrganismPopUp();

################################################################
### gene queries
print "<B><A HREF='help.get-orthologs.html#queries'>Query genes</A></B>&nbsp;";
print "<BR>\n";
print $query->textarea(-name=>'queries',
		       -default=>$default{queries},
		       -rows=>3,
		       -columns=>40);

### option to upload a file with the gene list from the client machine 
print "<BR>Upload gene list from file<BR>\n";
print $query->filefield(-name=>'uploaded_file',
			-default=>'',
			-size=>45,
			-maxlength=>200);
print "<p>\n";

################################################################
## Taxon of interest
print "<B><A HREF='help.get-organisms.html#taxon'>Taxon of interest</A>&nbsp;</B>\n";
#&RSAT::error::FatalError($default{taxon});
print $query->popup_menu(-name=>'taxon',
			 -default=>$default{taxon},
			 -Values=>[@supported_taxons]
			);
print "<p>\n";

################################################################
## Return fields + thresholds

print "<B><A HREF='help.get-organisms.html#return'>Return fields</A>&nbsp;</B>\n";
print "<ul>\n";
print "<table cellpadding=3 border=1>\n";
print "<tr><th>Field</th><th>Lower<br>Threshold</th><th>Upper<br>Threshold</th></tr>";
foreach my $field (@output_fields) {
    print "\n<tr>\n";
    print "<th align=left><a href='help.get-organisms.html#",$field,"'>";
    print $field_description{$field};
    print "</a></th>\n";

    foreach my $th ("lth", "uth") {
	my $param = $th."_".$field;
	my $default_param = "none";
	if (defined($default{$param})) {
	    $default_param = $default{$param};
	}
	print "<td align=center>";
	print $query->textfield(-name=>$param,
				-default=>$default_param,
				-size=>5);
	print "</td>";
    }
    print "</tr>\n";
}
print "</table>\n";
print "</ul>\n";


################################################################
### send results by email or display on the browser
print "<p>\n"; 
&SelectOutput();

################################################################
### action buttons
print "<UL><UL><TABLE>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

################################################################
### data for the demo 
print $query->start_multipart_form(-action=>"get-orthologs_form.cgi");
$demo_queries = "lexA\n";
$demo_queries .= "recA\n";
$demo_queries .= "uvrB\n";
print "<TD><B>";
print $query->hidden(-name=>'queries',-default=>$demo_queries);
print $query->hidden(-name=>'organism',-default=>"Escherichia_coli_K12");
print $query->hidden(-name=>'taxon',-default=>"Gammaproteobacteria");
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


print "<TD><B><A HREF='help.get-orthologs.html'>MANUAL</A></B></TD>\n";
#print "<TD><B><A HREF='tutorials/tut_get-orthologs.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@scmbb.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</BLOCKQUOTE>\n";

print $query->end_html;

exit(0);

