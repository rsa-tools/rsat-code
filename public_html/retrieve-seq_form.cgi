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
$default{sequence_format} = "fasta";
#$default{seq_label} = "gene identifier + organism + gene name";
$default{seq_label} = "gene name";
$default{organism} = "Saccharomyces cerevisiae";
$default{rm} = "";
$default{noorf} = "checked";
$default{imp_pos} = "checked";
$default{from} = "default";
$default{to} = "default";
$default{genes} = "selection";
$default{gene_selection} = "";
$default{sequence_type} = "upstream";
$default{feattype} = "CDS";
$default{single_multi_org} = "single";
$default{ids_only} = "";
# $default{gene_col} = 1;
# $default{org_col} = 2;

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
}

### print the form ###
&RSA_header("retrieve sequence", 'form');

### head
print "<CENTER>";
print "Returns upstream, downstream or ORF sequences for a list of genes<P>\n";
print "</CENTER>";
print "<b>Remark: If you want to retrieve sequences from an organism that is in the <a href='http://www.ensembl.org'>EnsEMBL</a> database, we recommand to use the <a href='retrieve-ensembl-seq_form.cgi'>retrieve-ensembl-seq</a> program instead</b><p>\n";

print $query->start_multipart_form(-action=>"retrieve-seq.cgi");


#print "<FONT FACE='Helvetica'>";

#### Single organism
if ($default{single_multi_org} eq 'single') {
    $CHECKED = "checked";
} else {
    $CHECKED = "";
}
print ("<INPUT TYPE='radio' NAME='single_multi_org' VALUE='single' $CHECKED>", 
       "<A HREF=help.retrieve-seq.html#single_org>",
       "<b>Single organism</b>",
       "</A>\n");
print "&nbsp;"x4, &OrganismPopUpString();
print "<p>\n";

#### Multiple organisms
if ($default{single_multi_org} eq 'multi') {
    $CHECKED = "checked";
} else {
    $CHECKED = "";
}
print ("<INPUT TYPE='radio' NAME='single_multi_org' VALUE='multi' $CHECKED>", 
       "<b>Multiple organisms</b>",
       " (2-column input, check <A HREF=help.retrieve-seq.html#multi_org>help</a> for format)",
       "\n"
      );

# ### Gene/organism columns
# print "&nbsp;"x10;
# print "<B><A HREF='help.retrieve-seq.html#gene_col'>Gene column</A></B>&nbsp;\n";
# print $query->textfield(-name=>'gene_col',
# 			-default=>$default{gene_col},
# 			-size=>5);

# print "&nbsp;&nbsp;";
# print "<B><A HREF='help.retrieve-seq.html#org_col'>Organism column</A></B>&nbsp;\n";
# print $query->textfield(-name=>'org_col',
# 			-default=>$default{org_col},
# 			-size=>5);
# print "<BR>\n";


## &OrganismPopUp;

### query (gene list)
print "<p>";
print "<B><A HREF='help.retrieve-seq.html#genes'>Genes</A></B>&nbsp;";
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

## IDs only
print "<br>", $query->checkbox(-name=>'ids_only',
			       -checked=>$default{ids_only},
			       -label=>'');
print "<a href=help.retrieve-seq.html#ids_only>Query contains only IDs (no synonyms)</a>";

print "</UL>\n";
print "<BR>\n";

#### feature type
print "<B><A HREF='help.retrieve-seq.html#feattype'>Reference feature type (reference coordinate for positions)</A></B>&nbsp;<br>";
print $query->radio_group(-name=>'feattype',
			  -values=>[@supported_feature_types],
			  -default=>$default{feattype});
print "<BR>\n";

### sequence type
print "<B><A HREF='help.retrieve-seq.html#sequence_type'>Sequence type</A></B>&nbsp;";
print $query->popup_menu(-name=>'sequence_type',
			 -Values=>['upstream','downstream','ORFs (unspliced)'],
			 -default=>$default{sequence_type});

### from to
print "&nbsp;&nbsp;";
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
print "&nbsp;<A HREF='help.retrieve-seq.html#noorf'><B>Prevent overlap with neighbour genes (noorf)</B></A>";
print "<BR>\n";

### Repeat masking
print $query->checkbox(-name=>'rm',
  		       -checked=>$default{rm},
  		       -label=>'');
print "&nbsp;<A HREF='help.retrieve-seq.html#rm'><B>Mask repeats</B></A>";
print "&nbsp;<A HREF='help.retrieve-seq.html#rm_list'><B>(only valid for organisms with annotated repeats)</B></A>";
print "<BR>\n";

### allows for imprecise postions
print $query->checkbox(-name=>'imp_pos',
  		       -checked=>$default{imp_pos},
  		       -label=>'');
print "&nbsp;<A HREF='help.retrieve-seq.html#imp_pos'><B>Admit imprecise positions</A></B>";
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
			 -Values=>['gene identifier', 
				   'gene name',
				   'gene identifier + name',
				   'gene identifier + organism + gene name',
				   'full identifier'
				   ],
			 -default=>$default{seq_label});
print "<BR>\n";

## Pass the taxon from get-orthologs for the further programs
if ($query->param('taxon')) {
  print $query->hidden(-name=>'taxon',-default=>$query->param('taxon'));
}

### send results by email or display on the browser
&SelectOutput("server");

### data for the demo 
@demo_genes = qw (DAL5 GAP1 MEP1 MEP2 PUT4 MEP3 DAL80);
$demo_genes = join "\n", @demo_genes;


### action buttons
print "<UL><UL><TABLE class='formbutton'>\n";
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
# $ENV{rsat_www} = 	'http://rsat.ulb.ac.be/rsat/';
print $query->hidden(-name=>'noorf',-default=>"");
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


#print "<TD><B><A HREF='demo.retrieve-seq.html'>DEMO</A></B></TD>\n";
print "<TD><B><A HREF='help.retrieve-seq.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='tutorials/tut_retrieve-seq.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

#print "</FONT>\n";

print $query->end_html;

exit(0);

