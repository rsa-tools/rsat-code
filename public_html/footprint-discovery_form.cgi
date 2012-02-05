#!/usr/bin/perl
################################################################
## this cgi script fills the HTML form for the program footprint-discovery
#!/usr/bin/perl
BEGIN {
    if ($0 =~ /([^(\/)]+)$/) {
	push @INC, "$`lib/";
    }
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";
use RSAT::Tree;


$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

################################################################
## Initialize parameters


## Output fields
my @output_fields = qw(ident
		       ali_len
		       mismat
		       gap_open
		       e_value
		       bit_sc
		       rank);
my %field_description = ();
$field_description{ident} = "Percentage of identity";
$field_description{ali_len} = "Alignment length";
$field_description{mismat} = "Number of mismatches";
$field_description{gap_open} = "Number of gap openings";
$field_description{e_value} = "E-value";
$field_description{bit_sc} = "Bit score";
$field_description{rank} = "Rank";
#$field_description{s_rank} = "target rank";


################################################################
### default values for filling the form

## Default values for dyad-analysis
%default = ();
&LoadDyadDefault(\%default);

## Default parameters for get-orthologs
&LoadGetOrthoDefault(\%default);

## Other default parameters
$default{dyads_filter} = 'checked';
$default{bg_model} = 'taxfreq';
$default{leaders} = '';
$default{uth_rank} = 50;
$default{to_matrix} = 1;

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
&RSA_header("footprint-discovery", "form");
print "<CENTER>";
print "Given one or several genes from a query organism, collect all orthologous genes for a given taxonomical level <br>and discover conserved elements in their promoters.<br>\n";
print "(Program developed by <A href='mailto:rekins\@bigre.ulb.ac.be'>Rekin's Janky</A>\n";
print "and <a href='mailto:jvanheld\@bigre.ulb.ac.be'>Jacques van Helden</A>).\n";
print "<br>Reference: <a target='_blank' href=\"http://www.biomedcentral.com/1471-2105/9/37\">Janky & van Helden, BMC Bioinformatics 2008, 9:37.</a>";
print "</CENTER>";
print "<BLOCKQUOTE>\n";

&ListDefaultParameters() if ($ENV{rsat_echo} >= 2);

print $query->start_multipart_form(-action=>"footprint-discovery.cgi");

################################################################
## Print the options for the selection of orthologs
print "<hr>";
&PrintOrthoSelectionSection();

### use predicted leader genes
print "<br>";
print $query->checkbox(-name=>'leaders',
		       -checked=>$default{leaders},
		       -label=>'');
print "<A HREF='help.footprint-discovery.html#leader'><B>\n";
print "predict operon leader genes";
print "</B></A>\n";

################################################################
#### Footprint-discovery-specific options
print "<hr>";
print "<p><b>Options for <i>dyad-analysis</i></b></p>";

print "<ul>";

### filtering dyads
print $query->checkbox(-name=>'dyads_filter',
		       -checked=>$default{dyads_filter},
		       -label=>'');
print "<a href='help.footprint-discovery.html#filtering'><b>\n";
print "dyad filtering</b></a>\n";
print "(only accept dyads having at least one occurrence in the promoter of the query gene)";

### use background model
print "<br>";
print "<b><a href='help.footprint-discovery.html#bg_model'>Background model</A>&nbsp;</b>\n";
print $query->popup_menu(-name=>'bg_model',
			 -Values=>['taxfreq', ## taxon-wise background model (taxfreq)
				   'monads'], ## background model estimated from the input sequence set
			 -default=>$default{bg_model});


################################################################
## Print dyad return fields
&PrintDyadReturnFields(no_matrix=>1);


# #### Convert patterns to matrix
# print $query->checkbox(-name=>'to_matrix',
# 		       -checked=>$default{to_matrix},
# 		       -label=>'');
# print "&nbsp;Convert assembled patterns to Position-Specific Scoring Matrices";
# print "<BR>";

print "</ul>";

################################################################
### send results by email or display on the browser
print "<hr>";
&SelectOutput('email', email_only=>1);

################################################################
### action buttons
print "<UL><UL><TABLE class = 'formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

################################################################
### data for the demo 
print $query->start_multipart_form(-action=>"footprint-discovery_form.cgi");
$demo_queries = "lexA\n";
#$demo_queries .= "recA\n";
#$demo_queries .= "uvrB\n";
print "<TD><B>";
print $query->hidden(-name=>'queries',-default=>$demo_queries);
print $query->hidden(-name=>'organism',-default=>"Escherichia_coli_K12");
print $query->hidden(-name=>'taxon',-default=>"Enterobacteriales");
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


print "<TD><B><A HREF='help.footprint-discovery.html'>MANUAL</A></B></TD>\n";
#print "<TD><B><A HREF='tutorials/tut_footprint-discovery.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</BLOCKQUOTE>\n";

print "<hr class=solid>";

print $query->end_html;

exit(0);

