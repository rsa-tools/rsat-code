#!/usr/bin/env perl
################################################################
## this cgi script fills the HTML form for the program get-orthologs
#!/usr/bin/env perl
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
my @output_fields = qw(ref_name
		       query_name
		       query_organism
		       ident
		       ali_len
		       mismat
		       gap_open
		       e_value
		       bit_sc
		       rank
		       s_rank
		     );
my %field_description = ();
$field_description{ref_name} = "Reference gene name";
$field_description{query_name} = "Query gene name";
$field_description{query_organism} = "Query organism";
$field_description{ident} = "Percentage of identity";
$field_description{ali_len} = "Alignment length";
$field_description{mismat} = "Number of mismatches";
$field_description{gap_open} = "Number of gap openings";
$field_description{e_value} = "E-value";
$field_description{bit_sc} = "Bit score";
$field_description{rank} = "Rank";
$field_description{s_rank} = "Reciprocal rank";


################################################################
### default values for get-orthologs
%default = ();
&LoadGetOrthoDefault(\%default);

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
&RSA_header("get-orthologs", "form");
print "<CENTER>";
print "Given a list of genes from a query organism and a taxon of interest, <br>return genes coding for similar proteins in each genome of the taxon.<br>\n";
print "Program developed by <A HREF='https://www.kuleuven.be/wieiswie/en/person/u0076845'>Rekin's Janky</A>\n";
print "and <A HREF='mailto:Jacques.van-Helden\@univ-amu.fr'>Jacques van Helden</A>).\n";
print "</CENTER>";
print "<BLOCKQUOTE>\n";

################################################################
## Display the form only if it is relveant for the organisms supported
## on this RSAT instance.
&check_phylo_tools();

################################################################
## Form header

print $query->start_multipart_form(-action=>"get-orthologs.cgi");

&ListDefaultParameters() if ($ENV{rsat_echo} >= 2);

################################################################
## Print the options for the selection of orthologs
print "<hr/>\n"; 
&PrintOrthoSelectionSection();

################################################################
## Return fields + thresholds
print "<hr/>\n";

print "<B><A HREF='help.get-orthologs.html#return'>Return fields</A>&nbsp;</B>\n";
print "<ul>\n";
print "<table cellpadding=3>\n";
print ("<tr>",
       "<th>Field</th>",
       "<th>Lower<br>Threshold</th>",
       "<th>Upper<br>Threshold</th>",
       "</tr>\n");
foreach my $field (@output_fields) {
    print "\n<tr>\n";
    my $return_field = "return_".$field;
    print "<th align=left>";
    print $query->checkbox(-name=>$return_field,
			   -checked=>$default{$return_field},
			   -label=>' ');
    print join "", "<a href='help.get-orthologs.html#",$field,"'>", $field_description{$field}, "</a>\n";
    print "</th>\n";
    if (($field eq "ref_name")||($field eq "query_name")||($field eq "query_organism")){
	print "<td></td>";
    } else {
	foreach my $th ("ortho_lth", "ortho_uth") {
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

print "<TD><B>";

print '<script>
function setDemo(){
    $("#reset").trigger("click");
    queries.value = "lexA\n";
    $("#organism").val("Escherichia_coli_K_12_substr__MG1655_uid57779");
    $("#organism_name").val("Escherichia coli K 12 substr  MG1655 uid57779");
    $("#taxon").val("Enterobacteriales");
    $("#taxon_name").val("Enterobacteriales");
}
</script>';

print '<button type="button" onclick="setDemo();">DEMO</button>';
print "</B></TD>\n";


print "<TD><B><A HREF='help.get-orthologs.html'>MANUAL</A></B></TD>\n";
#print "<TD><B><A HREF='tutorials/tut_get-orthologs.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:Jacques.van-Helden\@univ-amu.fr'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</BLOCKQUOTE>\n";

print $query->end_html;

exit(0);

