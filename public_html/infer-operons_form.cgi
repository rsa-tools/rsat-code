#!/usr/bin/perl
#### this cgi script fills the HTML form for the program retrieve-seq
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
$default{genes} = "selection";
$default{organism} = "Escherichia coli K12";
#$default{organism} = "Escherichia_coli_K_12_substr__MG1655_uid57779";
$default{dist_thr} = 55;
$default{min_gene_nb} = 2;
$default{return_leader} = "checked";
$default{return_trailer} = "";
$default{return_operon} = "checked";
$default{return_query} = "checked";
$default{return_name} = "checked";
$default{return_upstr_dist} = "checked";
$default{return_q_info} = "";
$default{return_up_info} = "";
$default{return_down_info} = "";
$default{return_gene_nb} = "checked";

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
}

## Output fields
my @output_fields = qw(query
		       name
		       leader
		       trailer
		       operon
		       upstr_dist
		       q_info
		       up_info
		       down_info
		       gene_nb
		       );
my %field_description = ();
$field_description{leader} = "Predicted operon leader gene";
$field_description{trailer} = "Predicted operon trailer gene";
$field_description{operon} = "Composition of the operon";
$field_description{query} = "Query";
$field_description{name} = "Query gene name";
$field_description{upstr_dist} = "Distance to upstream neighbour (negative for overlapping genes)";
$field_description{q_info} = "Detailed info on the query gene";
$field_description{up_info} = "Detailed info on the gene located upstream the query";
$field_description{down_info} = "Detailed info on the gene located downstream the query";
$field_description{gene_nb} = "Number of genes in the predicted operon";

### print the form ###
&RSA_header("infer operon", 'form');

### head
print "<CENTER>";
print "Infers the operon to which each coding gene of a given list belongs in a prokaryotic genome.";
print "<br>This program was developed by <a target=_blank href='http://www.bigre.ulb.ac.be/people/Members/rekins'>Rekins Janky</a> and <a target=_blank href=http://www.bigre.ulb.ac.be/Users/jvanheld/>Jacques van Helden</a>.</center>";
print "</CENTER>";


print $query->start_multipart_form(-action=>"infer-operons.cgi");


#### Single organism
# if ($default{single_multi_org} eq 'single') {
#     $CHECKED = "checked";
# } else {
#     $CHECKED = "";
# }
# print ("<INPUT TYPE='radio' NAME='single_multi_org' VALUE='single' $CHECKED>", 
#        "<A HREF=help.infer-operons.html#single_org>",
#        "<b>Single organism</b>",
#        "</A>\n");
# print "&nbsp;"x4, &OrganismPopUpString();
# print "<p>\n";

# #### Multiple organisms
# if ($default{single_multi_org} eq 'multi') {
#     $CHECKED = "checked";
# } else {
#     $CHECKED = "";
# }
# print ("<INPUT TYPE='radio' NAME='single_multi_org' VALUE='multi' $CHECKED>", 
#        "<A HREF=help.infer-operons.html#multi_org>",
#        "<b>Multiple organisms</b>",
#        "</a>\n"
#       );

&OrganismPopUp();

### query (gene list)
print "<p>";
print "<B><A HREF='help.infer-operons.html#genes'>Genes</A></B>&nbsp;";
print $query->radio_group(-name=>'genes',
			  -values=>['all','selection'],
			  -default=>$default{genes});

print "<BR>\n";
print "<UL>\n";

print $query->textarea(-name=>'gene_selection',
		       -default=>$default{gene_selection},
		       -rows=>6,
		       -columns=>65);
### option to upload a file with the gene list from the client machine 
print "<BR>Upload gene list from file<BR>\n";
print $query->filefield(-name=>'uploaded_file',
			-default=>'',
			-size=>45,
			-maxlength=>200);

### distance threshold
print "</UL><BR><HR>\n";
print "<B><A HREF='help.infer-operons.html#dist_thr'>Distance threshold (bp)</A></B>&nbsp;\n";
print $query->textfield(-name=>'dist_thr',
			-default=>$default{dist_thr},
			-size=>5);

print "&nbsp;"x5;

print "<B><A HREF='help.infer-operons.html#min_gene_nb'>Minimum number of genes</A></B>&nbsp;\n";
print $query->textfield(-name=>'min_gene_nb',
			-default=>$default{min_gene_nb},
			-size=>5);
print "<BR><HR>\n";

################################################################
#### Return fields
print "<p><B><A HREF='help.infer-operons.html#return'>Return fields</A></B>&nbsp;<br>\n";
my $i = 0;
foreach my $field (@output_fields) {
  my $return_field = "return_".$field;
  print $query->checkbox(-name=>$return_field,
			 -checked=>$default{$return_field},
			 -label=>' ');
  print join "", "<a href='help.infer-operons.html#",$field,"'>", $field_description{$field}, "</a>\n";
  print "<br>\n";
}

### send results by email or display on the browser
print "<BR><HR>\n";
&SelectOutput();

### data for the demo 
@demo_genes = qw (bioD bioA trpB trpE hisB metB);
$demo_genes = join "\n", @demo_genes;


### action buttons
print "<UL><UL><TABLE class='formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

## Demo 1: selected genes
print $query->start_multipart_form(-action=>"infer-operons_form.cgi");
print "<TD><B>";
print $query->hidden(-name=>'gene_selection',-default=>$demo_genes);
print $query->hidden(-name=>'organism',-default=>"Escherichia coli K12");
print $query->hidden(-name=>'dist_thr',-default=>"55");
print $query->hidden(-name=>'leader',-default=>"checked");
print $query->hidden(-name=>'query',-default=>"checked");
print $query->hidden(-name=>'operon',-default=>"checked");
print $query->submit(-label=>"DEMO 1 (selected genes)");
print "</B></TD>\n";
print $query->end_form;

## Demo 2: all genes
print $query->start_multipart_form(-action=>"infer-operons_form.cgi");
print "<TD><B>";
print $query->hidden(-name=>'genes',-default=>'all');
print $query->hidden(-name=>'organism',-default=>"Escherichia coli K12");
print $query->hidden(-name=>'dist_thr',-default=>"55");
print $query->hidden(-name=>'leader',-default=>"checked");
print $query->hidden(-name=>'query',-default=>"checked");
print $query->hidden(-name=>'operon',-default=>"checked");
print $query->submit(-label=>"DEMO 2 (all genes)");
print "</B></TD>\n";
print $query->end_form;


#print "<TD><B><A HREF='demo.infer-operons.html'>DEMO</A></B></TD>\n";
print "<TD><B><A HREF='help.infer-operons.html'>MANUAL</A></B></TD>\n";
#print "<TD><B><A HREF='tutorials/tut_infer-operons.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

#print "</FONT>\n";

print $query->end_html;

exit(0);

