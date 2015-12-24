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

### Read the CGI query
$query = new CGI;

################################################################
## Initialize parameters

my @homology_types = qw( ortholog paralog homeolog all );

## Output fields
my @output_fields = qw(
    target_id
    ref_org
    subtype
    query_id
    ident_target
    ident_query
    );

my %field_description = ();
$field_description{target_id} = "Target gene identifier";
$field_description{query_id} = "Query gene identifier";
$field_description{ref_org} = "Reference organism";
$field_description{subtype} = "Compara homology subtype";
$field_description{ident_target} = "%identity with respect to target length";
$field_description{ident_query} = "%identity with respect to query length";

################################################################
### default values for get-orthologs
%default = ();
&LoadGetOrthoComparaDefault(\%default);

### Replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
}

################################################################
### print the form ###


################################################################
### header
&RSA_header("get-orthologs-compara", "form");
print "<CENTER>";
print "Returns orthologues, plus optionally paralogues and homoeologues, for a
    set of genes in one or more organisms. <br>Relies on primary data from
    Ensembl Compara.<br><br>\n";
print "Program developed by <A HREF='bcontreras\@eead.csic.es'>Bruno Contreras-Moreira</A>\n";
print "and <A HREF='mailto:Jacques.van-Helden\@univ-amu.fr'>Jacques van Helden</A>.\n";
print "</CENTER>";
print "<BLOCKQUOTE>\n";

################################################################
## Display the form only if it is relveant for the organisms supported
## on this RSAT instance.
&check_compara_tools();

################################################################
## Form header

print $query->start_multipart_form(-action=>"get-orthologs-compara.cgi");

&ListDefaultParameters() if ($ENV{rsat_echo} >= 2);

################################################################
## Print the options for the selection of orthologs
print "<hr/>\n"; 
&PrintOrthoComparaSelectionSection();

## homology type
print "<B><A HREF='help.get-orthologs-compara.html#feattype'>Homology type</A></B>&nbsp;<br>";
print $query->radio_group(-name=>'type',
             -values=>[@homology_types],
             -default=>$default{feattype});
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
### send results by email or display on the browser
print "<hr/>\n"; 
&SelectOutput();

################################################################
### action buttons
print "<UL><UL><TABLE class = 'formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

################################################################
### data for the demo 
print $query->start_multipart_form(-action=>"get-orthologs-compara_form.cgi");
$demo_queries = "BRADI4G31367.1\n";
print "<TD><B>";
print $query->hidden(-name=>'queries',-default=>$demo_queries);
print $query->hidden(-name=>'organism',-default=>"triticum_aestivum");
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;

print "<TD><B><A HREF='help.get-orthologs-compara.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:Jacques.van-Helden\@univ-amu.fr'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</BLOCKQUOTE>\n";

print $query->end_html;

exit(0);

