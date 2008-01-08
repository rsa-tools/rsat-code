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
require "RSA.cgi.lib";
require "RSA2.cgi.lib";
use RSAT::Tree;


$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

################################################################
## Initialize parameters

@selected_organisms = qw(
			 Bacillus_subtilis
			 Brucella_melitensis
			 Escherichia_coli_O157H7
			 Escherichia_coli_K12
			 Pseudomonas_putida_KT2440
			 Porphyromonas_gingivalis_W83
			 Salmonella_typhimurium_LT2
			 Streptococcus_mutans
			 Streptomyces_coelicolor
			 Treponema_denticola_ATCC_35405
			 );

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
$default{organism} = "Escherichia_coli_K12";
$default{taxon} = "Gammaproteobacteria";
$default{queries} = '';
$default{full} = '';
$default{match_description} = '';
$default{feattype} = "CDS";
$default{lth_ali_len} = 50;
$default{uth_e_value} = 1e-5;
$default{uth_rank} = 1;
#$default{uth_s_rank} = 1;

$default{return_ident} = "checked";
$default{return_ali_len} = "";
$default{return_mismat} = "";
$default{return_gap_open} = "";
$default{return_e_value} = "checked";
$default{return_bit_sc} = "";
$default{return_rank} = "";
#$default{return_s_rank} = "";

$default{dyads_filter} = 'checked';
$default{bg_model} = 'taxfreq';
$default{leaders} = '';

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
print "Program developed by <A HREF='mailto:rekins\@scmbb.ulb.ac.be'>Rekin's Janky</A>\n";
print "and <A HREF='mailto:jvanheld\@scmbb.ulb.ac.be'>Jacques van Helden</A>).\n";
print "<br>Reference: Janky & van Helden (2008), BMC Bioinformatics, in press. ";
print "</CENTER>";
print "<BLOCKQUOTE>\n";

print $query->start_multipart_form(-action=>"footprint-discovery.cgi");

################################################################
#### choice of the organism
&OrganismPopUp(@selected_organisms);

################################################################
### gene queries
print "<B><A HREF='help.footprint-discovery.html#queries'>Query genes</A></B>&nbsp;";
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
&TaxonPopUp();

################################################################
#### Options
print "<BR>";
print "<HR width=550 align=left>\n";
print "<B>Options</B>&nbsp;";
print "<BR>";

### filtering dyads
print "&nbsp;" x 5;
print $query->checkbox(-name=>'dyads_filter',
		       -checked=>$default{dyads_filter},
		       -label=>'');
print "<A HREF='help.footprint-discovery.html#filtering'><B>\n";
print "filter dyads having occurrences in the query gene\n";
print "</B></A>\n";

### use predicted leader genes
print "<BR>";
print "&nbsp;" x 5;
print $query->checkbox(-name=>'leaders',
		       -checked=>$default{leaders},
		       -label=>'');
print "<A HREF='help.footprint-discovery.html#leader'><B>\n";
print "use prediction of upstream leader genes (slow process)\n";
print "</B></A>\n";

### use background model
print "<BR>";
print "&nbsp;" x 5;
print "<B><A HREF='help.footprint-discovery.html#bg_model'>Background model</A>&nbsp;</B>\n";
print $query->popup_menu(-name=>'bg_model',
			 -Values=>['taxfreq',#taxon-wise background model (taxfreq)',
				  'monad'],
			 -default=>$default{bg_model});

################################################################
### send results by email or display on the browser
print "<p>\n"; 
print "<HR width=550 align=left>\n";
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
print "<TD><B><A HREF='mailto:jvanheld\@scmbb.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</BLOCKQUOTE>\n";

print $query->end_html;

exit(0);

