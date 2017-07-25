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
$default{imp_pos} = "";
$default{from} = "default";
$default{to} = "default";
$default{genes} = "selection";
$default{gene_selection} = "";
$default{sequence_type} = "upstream";
$default{feattype} = "gene";
$default{single_multi_org} = "single";
$default{ids_only} = "";
# $default{gene_col} = 1;
# $default{org_col} = 2;

## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## TEMPORARY (2015-09): RESTRICT SUPPORTED FEATURE TYPES until the switch from NCBI
## to EnsemblGenomes as genome source is completely checked.
@supported_feature_types = qw(gene mRNA CDS);
##
## END TEMPORARY
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

#print $query->start_multipart_form(-action=>"retrieve-seq.cgi");
print "<form method='post' enctype='multipart/form-data' action='retrieve-seq.cgi'>";

print "<input type='text' id='menu_open' name='menu_open' style='display:none'>";

&ListParameters() if ($ENV{rsat_echo} >= 2);

#print "<FONT FACE='Helvetica'>";

#### Single organism
if ($default{single_multi_org} eq 'single') {
    $CHECKED = "checked";
} else {
    $CHECKED = "";
}
print ("<INPUT TYPE='radio' NAME='single_multi_org' VALUE='single' $CHECKED>", 
       "<A class='iframe' HREF=help.retrieve-seq.html#single_org>",
       "<b>Single organism</b>",
       "</A>\n");
#print "&nbsp;"x4, &OrganismPopUpSelectable();
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
       " (2-column input, check <A class='iframe' HREF=help.retrieve-seq.html#multi_org>help</a> for format)",
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
print "<B><A class='iframe' HREF='help.retrieve-seq.html#genes'>Genes</A></B>&nbsp;";
print $query->radio_group(-name=>'genes',
			  -values=>['all','selection'],
			  -default=>$default{genes});

print "<BR>\n";
print "<UL>\n";

#print $query->textarea(-name=>'gene_selection',
#		       -default=>$default{gene_selection},
#		       -rows=>6,
#		       -columns=>65);

print "<textarea id='gene_selection' name='gene_selection' rows='6' cols='65'>$default{gene_selection}</textarea>";

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
print "<a class='iframe' href=help.retrieve-seq.html#ids_only>Query contains only IDs (no synonyms)</a>";

print "</UL>\n";
print "<BR>\n";

#### feature type
print "<B><A class='iframe' HREF='help.retrieve-seq.html#feattype'>Reference feature type (reference coordinate for positions)</A></B>&nbsp;<br>";
print $query->radio_group(-name=>'feattype',
			  -values=>[@supported_feature_types],
			  -default=>$default{feattype});
print "<BR>\n";

### sequence type
print "<B><A class='iframe' HREF='help.retrieve-seq.html#sequence_type'>Sequence type</A></B>&nbsp;";
print $query->popup_menu(-name=>'sequence_type',
			 -Values=>['upstream','downstream','ORFs (unspliced)'],
			 -default=>$default{sequence_type});

### from to
print "&nbsp;&nbsp;";
print "<B><A class='iframe' HREF='help.retrieve-seq.html#from_to'>From</A></B>&nbsp;\n";
#print $query->textfield(-name=>'from',
#			-default=>$default{from},
#			-size=>5);

print "<input type='text' id='from' name='from' value=$default{from} size='5'/>";

print "&nbsp;&nbsp;";
print "<B><A class='iframe' HREF='help.retrieve-seq.html#from_to'>To</A></B>&nbsp;\n";
#print $query->textfield(-name=>'to',
#			-default=>$default{to},
#			-size=>5);
print "<input type='text' id='to' name='to' value=$default{to} size='5'/>";


print "<BR>\n";

### prevent ORF overlap
#print $query->checkbox(-name=>'noorf',
#  		       -checked=>$default{noorf},
#  		       -label=>'');
print "<input type='checkbox' name='noorf' id='noorf' checked='$default{noorf}' />";

print "&nbsp;<A class='iframe' HREF='help.retrieve-seq.html#noorf'><B>Prevent overlap with neighbour genes (noorf)</B></A>";
print "<BR>\n";

### Repeat masking
print $query->checkbox(-name=>'rm',
  		       -checked=>$default{rm},
  		       -label=>'');
print "&nbsp;<A class='iframe' HREF='help.retrieve-seq.html#rm'><B>Mask repeats</B></A>";
print "&nbsp;<A class='iframe' HREF='help.retrieve-seq.html#rm_list'><B>(only valid for organisms with annotated repeats)</B></A>";
print "<BR>\n";

################################################################
## Allows for imprecise postions
##
## 2014-03-03: JvH temporarily inactivates this option becauses it provokes a
## bug with Arabidopsis thaliana genome.
##
## Besides, I should revise the utility of this option.
##
# print $query->checkbox(-name=>'imp_pos',
#   		       -checked=>$default{imp_pos},
#   		       -label=>'');
# print "&nbsp;<A HREF='help.retrieve-seq.html#imp_pos'><B>Admit imprecise positions</A></B>";
# print "<BR>\n";

### sequence format 
print "<B><A class='iframe' HREF='help.retrieve-seq.html#formats'>Sequence format</A></B>&nbsp;";
print $query->popup_menu(-name=>'format',
			 -Values=>['fasta', 
				   'IG',
				   'wconsensus',
				   'multi'],
			 -default=>$default{sequence_format});
print "<BR>\n";

### sequence label
print "<B><A class='iframe' HREF='help.retrieve-seq.html#seq_label'>Sequence label</A></B>&nbsp;";
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
$demo_genes = join "\\n", @demo_genes;


### action buttons
print "<UL><UL><TABLE class='formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset(-id=>"reset"), "</TD>\n";



print $query->end_form;



print "<TD><B>";

print "<script>
function setDemo(demo_genes){
    \$('#reset').trigger('click');
    \$('#gene_selection').val(demo_genes);
    \$('#from').val('-800');
    \$('#to').val('-1');
    \$('#noorf').removeAttr('checked');
}
</script>";
print '<button type="button" onclick="setDemo('. "'$demo_genes'" .')">DEMO</button>';
print "</B></TD>\n";


#print "<TD><B><A HREF='demo.retrieve-seq.html'>DEMO</A></B></TD>\n";
print "<TD><B><A class='iframe' HREF='help.retrieve-seq.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='htmllink.cgi?title=RSAT : Tutorials&file=tutorials/tut_retrieve-seq.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:Jacques.van-Helden\@univ-amu.fr'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

#print "</FONT>\n";


print $query->end_html;

exit(0);

