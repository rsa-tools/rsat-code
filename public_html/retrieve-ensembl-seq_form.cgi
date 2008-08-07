#!/usr/bin/perl
#### this cgi script fills the HTML form for the program retrieve-ensembl-seq
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

use DBI();

### Read the CGI query
$query = new CGI;

### default values for filling the form
#$default{sequence_format} = "fasta";
#$default{seq_label} = "gene identifier + organism + gene name";
#$default{seq_label} = "gene name";
$default{organism} = "Homo sapiens";
$default{rm} = "";
$default{noorf} = "";
#$default{imp_pos} = "checked";
$default{from} = "-2000";
$default{to} = "-1";
$default{genes} = "selection";
$default{gene_selection} = "";
$default{sequence_type} = "upstream";
$default{feattype} = "mRNA";
$default{feature} = "gene";
$default{alltranscripts} = "";
$default{single_multi_org} = "single";
# $default{ids_only} = "";
# $default{gene_col} = 1;
# $default{org_col} = 2;
$default{taxon_selection} = "";
$default{homology_selection} = "";
$default{header_org} = "scientific";

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
}

### print the form ###
&RSA_header("retrieve-ensembl-seq", 'form');

### head
print "<CENTER>";
print "Returns upstream, downstream, intronic, exonic or UTR sequences for a list of genes from the EnsEMBL database. <P>\n";
print "Program developed by <A HREF='mailto:oly\@bigre.ulb.ac.be (Olivier Sand)'>Olivier Sand</A> with the help of <A HREF='mailto:morgane\@bigre.ulb.ac.be (Morgane Thomas-Chollier)'>Morgane Thomas-Chollier</A><P>";
print "</CENTER>";


print $query->start_multipart_form(-action=>"retrieve-ensembl-seq.cgi");


#print "<FONT FACE='Helvetica'>";

#### Query organism list
my @selected_organisms;

eval {
	my $dbh = DBI->connect("DBI:mysql:host=ensembldb.ensembl.org;mysql_connect_timeout=10", "anonymous", "", {'RaiseError' => 1});
};

if ($@) {
	print "<p><font size=2 color=red>No answer from the EnsEMBL database ; server may be down. Try again later...</font></p>".$query->end_html;
	exit(0);
} else {

my $dbh = DBI->connect("DBI:mysql:host=ensembldb.ensembl.org", "anonymous", "", {'RaiseError' => 0});
my $sth = $dbh->prepare("SHOW DATABASES");
$sth->execute();
my $previous_org = "bogus";
while (my $ref = $sth->fetchrow_hashref()) {
    if ($ref->{'Database'} =~ /_core_\d+/) {
	$dbversion = $ref->{'Database'};
	$dbversion =~ s/.+_core_//;
	$dbversion =~ s/_.+//;
	$ref->{'Database'} =~s/_core_.+//;
	if ($ref->{'Database'} ne $previous_org) {
	    push @selected_organisms, $ref->{'Database'};
	    $previous_org = $ref->{'Database'};
	}

    }
}
$sth->finish();
$dbh->disconnect();

my $organismPopup = "";
$organismPopup .=  "<BR/><B>Query organism</B>&nbsp;";
$organismPopup .=  "<SELECT NAME='organism'>\n";
foreach my $org (@selected_organisms) {
    $name = ucfirst($org);
    $name =~s/_/ /;
    if ((lc($org) eq lc($default{organism})) ||
	(lc($name) eq lc($default{organism}))) {
	$organismPopup .=  "<OPTION SELECTED VALUE='$org'>$name\n";
    } else {
	$organismPopup .=  "<OPTION VALUE='$org'>$name\n";
    }
}
$organismPopup .=  "</SELECT>";

print $organismPopup;

#### Database version
print "&nbsp;"x4;
print "<A HREF='http://www.ensembl.org'>EnsEMBL</A> database version: <B>";
print $dbversion;
print "</B><BR/>\n";

#### Single organism
if ($default{single_multi_org} eq 'single') {
    $CHECKED = "checked";
} else {
    $CHECKED = "";
}
print ("<INPUT TYPE='radio' NAME='single_multi_org' VALUE='single' $CHECKED>", 
       "<A HREF=help.retrieve-ensembl-seq.html#single_org>",
       "<b>Single organism</b>",
       "</A>\n");
print "<BR/>";

#### Multiple organisms
if ($default{single_multi_org} eq 'multi') {
    $CHECKED = "checked";
} else {
    $CHECKED = "";
}
print ("<INPUT TYPE='radio' NAME='single_multi_org' VALUE='multi' $CHECKED>", 
       "<A HREF=help.retrieve-ensembl-seq.html#multi_org>",
       "<b>Multiple organisms</b>",
       "</a>\n"
      );

print "&nbsp;"x10;
print "<B>Optional filters</B>";
print "&nbsp;"x4;

#### Taxon
#print "<p>";
print "<B><A HREF='help.retrieve-ensembl-seq.html#taxon'>Taxon</A></B>&nbsp;";
#print "<BR/>\n";
print $query->textfield(-name=>'taxon_selection',
		       -default=>$default{taxon_selection},
		       -size=>20);

#### Homology type
#print "<p>";
print "&nbsp;"x4;
print "<B><A HREF='help.retrieve-ensembl-seq.html#homology'>Homology type</A></B>&nbsp;";
#print "<BR/>\n";
print $query->textfield(-name=>'homology_selection',
		       -default=>$default{homology_selection},
		       -size=>20);

### query (gene list)
print "<p>";
print "<B><A HREF='help.retrieve-ensembl-seq.html#genes'>EnsEMBL IDs (gene, transcript or protein IDs)</A></B>&nbsp;";
#print $query->radio_group(-name=>'genes',
#			  -values=>['all','selection'],
#			  -default=>$default{genes});

print "<BR/>\n";
#print "<UL>\n";

print $query->textarea(-name=>'gene_selection',
		       -default=>$default{gene_selection},
		       -rows=>6,
		       -columns=>65);
### option to upload a file with the gene list from the client machine 
print "<BR>Upload gene list (EnsEMBL IDs) from file<BR>\n";
print $query->filefield(-name=>'uploaded_file',
			-default=>'',
			-size=>45,
			-maxlength=>200);

## IDs only
#print "<br>", $query->checkbox(-name=>'ids_only',
#			       -checked=>$default{ids_only},
#			       -label=>'');
#print "<a href=help.retrieve-seq.html#ids_only>Query contains only IDs (no synonyms)</a>";

print "</UL>\n";
print "<BR/>\n";

print "<P><B>Choose the type of sequence to retrieve:</B></P>";

print "<TABLE border=1>\n";
print "<TR align='center'>\n";
print "<TD>\n";
print "<B>Upstream or downstream sequence</B>";
print "</TD>\n";
print "<TD>\n";
print "<B>Feature sequence</B>";
print "</TD>\n";
print "</TR>\n";
print "<TR valign='top'>\n";
print "<TD>\n";


### sequence type
print "<B><A HREF='help.retrieve-ensembl-seq.html#sequence_type'>Sequence type</A></B>&nbsp;";
print $query->popup_menu(-name=>'sequence_type',
			 -Values=>['upstream','downstream'],
			 -default=>$default{sequence_type});
print "&nbsp;"x4;

### from to
print "<B><A HREF='help.retrieve-ensembl-seq.html#from_to'>From</A></B>&nbsp;\n";
print $query->textfield(-name=>'from',
			-default=>$default{from},
			-size=>5);

print "&nbsp;&nbsp;";
print "<B><A HREF='help.retrieve-ensembl-seq.html#from_to'>To</A></B>&nbsp;\n";
print $query->textfield(-name=>'to',
			-default=>$default{to},
			-size=>5);
print "<BR>\n";

#### Reference feature
print "<B><A HREF='help.retrieve-seq.html#feattype'>Relative to feature</A></B>&nbsp;";
print $query->radio_group(-name=>'feattype',
			  -values=>['gene','mRNA','CDS'],
			  -default=>$default{feattype});
print "<BR/>\n";

### All transcripts
print $query->checkbox(-name=>'alltranscripts',
  		       -checked=>$default{alltranscripts},
  		       -label=>'');
print "&nbsp;<A HREF='help.retrieve-ensembl-seq.html#alltranscripts'><B>Retrieve sequence relative to each alternative transcript (with mRNA feature)</B></A>";
print "<BR>\n";

### prevent overlap
print "<B><A HREF='help.retrieve-ensembl-seq.html#prevent_overlap'>Prevent overlap with neighbouring</A></B>&nbsp;";
print $query->popup_menu(-name=>'prevent_overlap',
			 -Values=>['none','ORF','gene'],
			 -default=>$default{prevent_overlap});
print "<BR>\n";

### prevent ORF overlap
#print $query->checkbox(-name=>'noorf',
#  		       -checked=>$default{noorf},
#  		       -label=>'');
#print "&nbsp;<A HREF='help.retrieve-ensembl-seq.html#noorf'><B>Prevent overlap with neighbour orfs (noorf)</B></A>";
#print "<BR>\n";

### prevent gene overlap
#print $query->checkbox(-name=>'nogene',
#  		       -checked=>$default{nogene},
#  		       -label=>'');
#print "&nbsp;<A HREF='help.retrieve-ensembl-seq.html#nogene'><B>Prevent overlap with neighbour genes (nogene)</B></A>";
#print "<BR>\n";

print "</TD>\n";

print "<TD>\n";

#### Feature to retrieve
print "<B><A HREF='help.retrieve-seq.html#feature'>Feature to retrieve</A></B>&nbsp;";
print $query->popup_menu(-name=>'feature',
			  -values=>['gene','introns','first intron','exons','non-coding exons','UTR'],
			  -default=>$default{feature});
print "<BR/>\n";
print "</TD>\n";

print "</TR>\n";
print "</TABLE>\n";

### Repeat masking
print $query->checkbox(-name=>'rm',
  		       -checked=>$default{rm},
  		       -label=>'');
print "&nbsp;<A HREF='help.retrieve-ensembl-seq.html#rm'><B>Mask repeats</B></A>";
print "&nbsp;<A HREF='help.retrieve-ensembl-seq.html#rm_list'><B>(only valid for organisms with annotated repeats)</B></A>";
print "<BR>\n";

### Mask coding
print $query->checkbox(-name=>'maskcoding',
  		       -checked=>$default{maskcoding},
  		       -label=>'');
print "&nbsp;<A HREF='help.retrieve-ensembl-seq.html#maskcoding'><B>Mask coding sequences</B></A>";
print "<BR>\n";

### Organism in header
print "<B><A HREF='help.retrieve-ensembl-seq.html#header_org'>Organism name in sequence header</A></B>&nbsp;";
print $query->popup_menu(-name=>'header_org',
			 -Values=>['scientific','common','none'],
			 -default=>$default{prevent_overlap});
print "<BR>\n";

### send results by email or display on the browser
&SelectOutput("server");

### data for the demo 
@demo_genes = qw (ENSG00000139618 ENSG00000138411);
$demo_genes = join "\n", @demo_genes;


### action buttons
print "<UL><UL><TABLE class='formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

print $query->start_multipart_form(-action=>"retrieve-ensembl-seq_form.cgi");
print "<TD><B>";
print $query->hidden(-name=>'gene_selection',-default=>$demo_genes);
print $query->hidden(-name=>'organism',-default=>"Homo sapiens");
print $query->hidden(-name=>'from',-default=>"-2000");
print $query->hidden(-name=>'to',-default=>"-1");
# $ENV{rsat_www} = 	'http://rsat.scmbb.ulb.ac.be/rsat/';
print $query->hidden(-name=>'noorf',-default=>"");
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


#print "<TD><B><A HREF='demo.retrieve-ensembl-seq.html'>DEMO</A></B></TD>\n";
print "<TD><B><A HREF='help.retrieve-ensembl-seq.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='tutorials/tut_retrieve-ensembl-seq.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

#print "</FONT>\n";

print $query->end_html;

exit(0);
}
