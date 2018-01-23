#!/usr/bin/perl
################################################################
## this cgi script fills the HTML form for the program get-orthologs
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

## Read the CGI query
$query = new CGI;

################################################################
## Initialize parameters

my @homology_types = qw( ortholog paralog homeolog all );

## Output fields
my @output_fields = qw(
    target_id
    ref_organism
    subtype
    query_id
    query_organism
    ident_target
    ident_query
    );

my %field_description = ();
$field_description{target_id} = "Target gene identifier";
$field_description{ref_organism} = "Reference organism";
$field_description{subtype} = "Compara homology subtype";
$field_description{query_id} = "Query gene identifier";
$field_description{query_organism} = "Query organism";
$field_description{subtype} = "Compara homology subtype";
$field_description{ident_target} = "%identity with respect to target length";
$field_description{ident_query} = "%identity with respect to query length";

################################################################
## default values for get-orthologs
%default = ();
&LoadGetOrthoComparaDefault(\%default);

## Replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
}

################################################################
## Print the form
################################################################

################################################################
## header
&RSA_header("get-orthologs-compara", "form");
print "<CENTER>";
print "Returns orthologues, plus optionally paralogues and homeologues, for a
    set of genes in one or more organisms. <br>Relies on primary data from
    Ensembl Compara.<br><br>\n";
print "Program developed by <A HREF='mailto:bcontreras\@eead.csic.es'>Bruno Contreras-Moreira</A>\n";
print "and <A HREF='mailto:Jacques.van-Helden\@univ-amu.fr'>Jacques van Helden</A>.\n";
print "</CENTER>";
print "<BLOCKQUOTE>\n";

&ListParameters() if ($ENV{rsat_echo} >= 2); ## For debugging

################################################################
## Display the form only if it is relveant for the organisms supported
## on this RSAT instance.
&check_compara_tools();

################################################################
## Print the demo description if specified
if ($default{demo_descr}) {
    print "<font color='darkgreen'>\n";
    print "<p><b>Demo:</b> ", $default{demo_descr}, "</p>\n";
    print "</font>\n";
}


################################################################
## Form header

print $query->start_multipart_form(-action=>"get-orthologs-compara.cgi");

&ListDefaultParameters() if ($ENV{rsat_echo} >= 2);

################################################################
## Print the options for the selection of orthologs
print "<hr/>\n"; 
&PrintOrthoComparaSelectionSection();

## homology type
print "<B><A HREF='help.get-orthologs-compara.html#type'>Homology type</A></B>&nbsp;<br>";
print $query->radio_group(-name=>'type',
             -values=>[@homology_types],
             -default=>$default{type});
print "<BR>\n";


################################################################
## thresholds
print "<hr/>\n";

print "<B><A HREF='help.get-orthologs-compara.html#return'>Return fields</A>&nbsp;</B>\n";
print "<ul>\n";
print "<table cellpadding=3>\n";
print ("<tr>",
       "<th>Field</th>",
       "<th>Lower</th>",
       "</tr>\n");
my $field_number = 0;
foreach my $field (@output_fields) {
    print "\n<tr>\n";
    my $return_field = "return_".$field;
    print "<th align=left>";
    printf("%d. ",++$field_number);
    print join "", "<a href='help.get-orthologs-compara.html#",$field,"'>", 
        $field_description{$field}, "</a>\n";
    print "</th>\n";
    if($field eq "ref_org" || $field eq "query_id" || $field eq "target_id" || $field eq "subtype"){
	    print "<td></td>";
    } 
    elsif($field eq "ident_target") {
        my $default_param = 0;
        if (defined($default{$field})) {
            $default_param = $default{$field};
        }
        print "<td align=center>";
        print $query->textfield(-name=>$field,
                    -default=>$default_param,
                    -size=>5);
        print "</td>";
    } 
    elsif($field eq "ident_query") {
        my $default_param = 0;
        if (defined($default{$field})) {
            $default_param = $default{$field};
        }
        print "<td align=center>";
        print $query->textfield(-name=>$field,
                    -default=>$default_param,
                    -size=>5);
        print "</td>";
    }
    print "</tr>\n";
}
print "</table>\n";
print "</ul>\n";


################################################################
## Send results by email or display on the browser
print "<hr/>\n"; 
&SelectOutput();

################################################################
## Action buttons
print "<UL><UL><TABLE class = 'formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

################################################################
## Data for the demo on ortholog searches
print $query->start_multipart_form(-action=>"get-orthologs-compara_form.cgi");
my $demo_descr = "Search orthologs for gene FT1 (Bradi1g48830) from <i>Brachypodium distachyon</i> in several grasses.";
print "<TD><B>";
print $query->hidden(-name=>'queries',-default=>"BRADI4G31367.1");
print $query->hidden(-name=>'type',-default=>"ortholog");
print $query->hidden(-name=>'demo_descr',-default=>$demo_descr);
print $query->hidden(-name=>'organism',
-default=>"brachypodium_distachyon\nhordeum_vulgare\noryza_indica\noryza_sativa\nsetaria_italica\nsorghum_bicolor\ntriticum_aestivum\ntriticum_urartu\zea_mays");
print $query->submit(-label=>"DEMO 1 (FT1 orthologs)");
print "</B></TD>\n";
print $query->end_form;

################################################################
## Data for the demo on ortholog searches
print $query->start_multipart_form(-action=>"get-orthologs-compara_form.cgi");
my $demo_descr2 = "Search paralogs for a gene from <i>Arabidopsis thaliana</i> in its own genome.";
print "<TD><B>";
print $query->hidden(-name=>'queries',-default=>"AT5G45730.1");
print $query->hidden(-name=>'type',-default=>"paralog");
print $query->hidden(-name=>'ident_target',-default=>"0");
print $query->hidden(-name=>'ident_query',-default=>"0");
print $query->hidden(-name=>'demo_descr',-default=>$demo_descr2);
print $query->hidden(-name=>'organism',-default=>"arabidopsis_thaliana");
print $query->submit(-label=>"DEMO 2 (inparalogs)");
print "</B></TD>\n";
print $query->end_form;

print "<TD><B><A HREF='help.get-orthologs-compara.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:Jacques.van-Helden\@univ-amu.fr'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</BLOCKQUOTE>\n";

print $query->end_html;

exit(0);

