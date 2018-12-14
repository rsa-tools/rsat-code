#!/usr/bin/env perl
#### this cgi script fills the HTML form for the program matrix-clustering
BEGIN {
    if ($0 =~ /([^(\/)]+)$/) {
        push (@INC, "$`lib/");
    }
    require "RSA.lib";
}
#if ($0 =~ /([^(\/)]+)$/) {
#    push (@INC, "$`lib/");
#}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";
use RSAT::matrix;
use RSAT::MatrixReader;

$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

### default values for filling the form
$default{sequence_format} = "fasta";
#$default{seq_label} = "gene identifier + organism + gene name";
$default{seq_label} = "gene name";
$default{organism} = "";
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

################################################################
### print the form ###

################################################################
### header
&RSA_header_bootstrap("retrieve sequence", "form");

print $query->start_multipart_form(-action=>"retrieve-seq.cgi");

print '
<!-- Form with bootstrap -->
<div class="container">
<div class="row">
<div class="col-lg-9 col-md-5 col-sm-8 col-xs-9 bhoechie-tab-container">
<div class="col-lg-2 col-md-3 col-sm-3 col-xs-3 bhoechie-tab-menu">
<div class="list-group">
<a href="#" class="list-group-item active text-center">
<h4 class="glyphicon"><i class="fa fa-info-circle fa-2x"></i></h4><br/>retrieve sequence
</a>
<a href="#" class="list-group-item text-center">
<h4 class="glyphicon"><i class="fa fa-tag fa-2x"></i></h4><br/>Mandatory options
</a>
<a href="#" class="list-group-item text-center">
<h4 class="glyphicon"><i class="fa fa-tasks fa-2x"></i></h4><br/>Advanced options
</a>
<a href="#" class="list-group-item text-center">
<h4 class="glyphicon"><i class="fa fa-play-circle fa-2x"></i></h4><br/>Run analysis
</a>
</div>
</div>
<div class="col-lg-9 col-md-9 col-sm-9 col-xs-9 bhoechie-tab">


<!-- ################################################################ -->
<!-- ### info ### -->

<div class="bhoechie-tab-content active">

<h2> <img src="images/RSAT_logo.jpg" style="max-width:150px;max-height:60px;padding-bottom:10px" alt="RSAT server" border="0"></img>
retrieve sequence</h2>
<span class="fa-stack fa-lg">
<i class="fa fa-info-circle fa-stack-1x"></i>
</span>
Starting from a list of genes, returns upstream, downstream or ORF sequences. Dedicated to genomes <b>locally-installed</b> in RSAT.<br>
To retrieve sequences from an organism that is in the EnsEMBL database, we recommend to use the <a href="retrieve-ensembl-seq_form.cgi">retrieve-ensembl-seq</a> program instead. <br>
<span class="fa-stack fa-lg">
<i class="fa fa-user fa-stack-1x"></i>
</span>

<a target="_blank" href="http://jacques.van-helden.perso.luminy.univ-amu.fr/ ">Jacques van Helden</a>.<br>

<span class="fa-stack fa-lg">
<i class="fa fa-folder-open fa-stack-1x"></i>
</span>
<a href="sample_outputs/retrieve-seq_demo20180226.fasta">Sample output</a><br>

<span class="fa-stack fa-lg">
<i class="fa fa-book fa-stack-1x"></i>
</span>
<a class="iframe" href="help.retrieve-seq.html">User Manual</a><br>

<span class="fa-stack fa-lg">
<i class="fa fa-graduation-cap fa-stack-1x"></i>
</span>
<a class="iframe" href="tutorials/tut_retrieve-seq.html">Tutorial</a><br>

<span class="fa-stack fa-lg">
<i class="fa fa-twitter fa-stack-1x"></i>
</span>
<a href="https://twitter.com/rsatools" target="_blank">Ask a question to the RSAT team</a><br>

<span class="fa-stack fa-lg">
<i class="fa fa-pencil fa-stack-1x"></i>
</span>
Cite the publication: <a href="https://twitter.com/rsatools" target="_blank"></a><br>
<div class="panel panel-default">
<div class="panel-body">
Medina-Rivera A*, Defrance M*, Sand O*, Herrmann C, Castro-Mondragon J, Delerce J, Jaeger S, Blanchet C, Vincens P, Caron C, Staines DM, Contreras-Moreira B, Artufel M, Charbonnier - Khamvongsa L, Hernandez C, Thieffry D, Thomas-Chollier M., van Helden J. <i>"RSAT 2015 : Regulatory Sequence Analysis Tools"</i>, Nucleic Acid Research 43(W1):W50-W56 (2015) [<a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4489296/" target="_blank">Pubmed</a>][<a href="https://academic.oup.com/nar/article/43/W1/W50/2467887" target="_blank">Full text</a>]
<br>
van Helden, J., Andre, B. & Collado-Vides, J. (2000). <i>A web site for the computational analysis of yeast regulatory sequences</i>. Yeast 16(2), 177-187. [<a href="https://www.ncbi.nlm.nih.gov/pubmed/10641039" target="_blank">Pubmed 10641039</a>]

</div>
</div>
</div>

<!-- ################################################################ -->
<!-- ### mandatory options ### -->
<div class="bhoechie-tab-content">
<!-- Organism -->
<div class="panel panel-danger">
<div class="panel-heading">Organism <i class="fa fa-info-circle" data-container="body" data-trigger="hover" rel="popover" data-placement="right" data-content="Select the organism from which to retrieve the sequences, or select the option to retrieve sequences from multiple organisms. The organism name is then specified below"></i></div>

<div class="panel-body">
<div class="form-group">
';

#### Single organism
if ($default{single_multi_org} eq 'single') {
    $CHECKED = "checked";
} else {
    $CHECKED = "";
}
print ("<INPUT TYPE='radio' NAME='single_multi_org' VALUE='single' $CHECKED>",
"<b>Single organism </b>");
print '<a class="badge badge-primary iframe" HREF="help.retrieve-seq.html#single_org">Info</a><br/><br/>';

print "&nbsp;"x4, &OrganismPopUpString();
print "<p><br/>";

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

print'

</div>
</div>
</div>

<!-- Gene -->
<div class="panel panel-danger">
<div class="panel-heading">Gene
<i class="fa fa-info-circle" rel="popover" data-container="body" data-trigger="hover" data-placement="right" data-content="Specify the genes relative to which the sequences will be retrieved"></i>
</div>
<div class="panel-body">
<div class="form-group">';
if($default{genes} eq "all"){
    $CHECKED = "checked";
}else{
    $CHECKED = "";
}
print "<input type='radio' name='genes' value='all' $CHECKED>all genes of this organism<br/>";
if($default{genes} eq "selection"){
    $CHECKED = "checked";
}else{
    $CHECKED = "";
}
print "<input type='radio' name='genes' value='selection' $CHECKED>genes specified below ";
print '<a class="badge badge-primary iframe" HREF="help.retrieve-seq.html#genes">Info</a>';

print "<BR><br>\n";
print $query->textarea(-id=>"gene_selection", -name=>"gene_selection", -default=>$default{gene_selection}, -rows=>5, -columns=>60, class=>"form-control");

### option to upload a file with the gene list from the client machine
print "<BR>Upload gene list from file<BR>\n";
print $query->filefield(-name=>'uploaded_file',
-default=>'',
-size=>45,
-maxlength=>200);

print '</div></div></div>';


print '
<!-- Sequence type -->
<div class="panel panel-danger">
<div class="panel-heading">Sequence type
<i class="fa fa-info-circle" data-container="body" data-trigger="hover" rel="popover" data-placement="right" data-content="Precise the type of sequence and its localisation relative to the gene"></i>
</div>
<div class="panel-body">
<div class="form-group">
<a class="badge badge-primary iframe" HREF="help.retrieve-seq.html#sequence_type">Info</a>';

print $query->popup_menu(-name=>'sequence_type', class=>'form-control',
-Values=>['upstream','downstream','ORFs (unspliced)'],
-default=>$default{sequence_type});

### from to
print "&nbsp;"x10;
print "<B><A class='iframe' HREF='help.retrieve-seq.html#from_to'>From</A></B>&nbsp;\n";
print "<input type='text' id='from' name='from' value=$default{from} class='form-control' style='width:80px'/>";
print "&nbsp;&nbsp;";
print "<B><A class='iframe' HREF='help.retrieve-seq.html#from_to'>To</A></B>&nbsp;\n";
print "<input type='text' id='to' name='to' value=$default{to} style='width:80px' class='form-control'/>";
print '</div></div></div></div>
<!-- ################################################################-->
<!-- ### advanced options ###-->

<!-- ADVANCED OPTIONS -->
<div class="bhoechie-tab-content">

<div class="panel panel-warning">
<div class="panel-heading">Advanced options</div>
<div class="panel-body"> <br>';

print $query->checkbox(-name=>'ids_only',
-checked=>$default{ids_only},
-label=>'');

print 'Query contains only IDs (no synonyms) <a class="badge badge-primary iframe" HREF="help.retrieve-seq.html#ids_only">Info</a><br/><br/>';

print "Reference feature type (reference coordinate for positions) <a class='badge badge-primary iframe' HREF='help.retrieve-seq.html#feattype'>Info</a><br>";
print $query->radio_group(-name=>'feattype',
-values=>[@supported_feature_types],
-default=>$default{feattype});
print "<BR><br/>";
print "<input type='checkbox' name='noorf' id='noorf' checked='$default{noorf}' />";
print "&nbsp;Prevent overlap with neighbour genes (noorf) <a class='badge badge-primary iframe' HREF='help.retrieve-seq.html#noorf'>Info</a>";
print "<BR><br/>";

### Repeat masking
print $query->checkbox(-name=>'rm',
-checked=>$default{rm},
-label=>'');
print "&nbsp;Mask repeats";
print "&nbsp;(only valid for organisms with annotated repeats) <a class='badge badge-primary iframe' HREF='help.retrieve-seq.html#rm_list'>Info</a>";

print '</div>
</div>
</div>

<!--################################################################-->
<!--### output & run ###-->

<div class="bhoechie-tab-content">

<!-- ## Specific options for output files-->
<div class="panel panel-info">
<div class="panel-heading">Output options</div>

<div class="panel-body">
<div class="form-group">';

print "Sequence format&nbsp;";
print $query->popup_menu(-name=>'format', class=>'form-control',
-Values=>['fasta',
'IG',
'wconsensus',
'multi'],
-default=>$default{sequence_format});
print " <a class='badge badge-primary iframe' HREF='help.retrieve-seq.html#formats'>Info</a><BR>\n";

### sequence label
print "Sequence label&nbsp;";
print $query->popup_menu(-name=>'seq_label', class=>'form-control',
-Values=>['gene identifier',
'gene name',
'gene identifier + name',
'gene identifier + organism + gene name',
'full identifier'
],
-default=>$default{seq_label});
print " <a class='badge badge-primary iframe' HREF='help.retrieve-seq.html#seq_label'>Info</a><BR>\n";

## Pass the taxon from get-orthologs for the further programs
if ($query->param('taxon')) {
    print $query->hidden(-name=>'taxon',-default=>$query->param('taxon'));
}

### send results by email or display on the browser
print '<hr/>';
&SelectOutput("server");
print " </div>
</div> </div>";


################################################################
## Action buttons

print '<script> function formreset(){
demo_descr.innerHTML = "";
} </script>';
print $query->submit(-label=>"GO", -class=>"btn btn-success", -type=>"button");
print " ";
print "<input type='reset' id='reset' class='btn btn-warning' onclick='formreset()' value='RESET'>";
print $query->end_form;
print ' </div>
</div>
</div>
';

################################################################
## Demo area
print "<textarea id='demo' style='display:none'></textarea>";
print "<div id='demo_descr' class='col-lg-9 col-md-5 col-sm-8 col-xs-9 demo-buttons-container'></div>";
print "</div> ";

### data for the demo
@demo_genes = qw (DAL5 GAP1 MEP1 MEP2 PUT4 MEP3 DAL80);
$demo_genes = join "\\n", @demo_genes;

print '<script>
function setDemo1(demo_genes){
    $("#reset").trigger("click");
    descr_1 = "<H4>Demonstration</H4>\n \
    <blockquote class =\'blockquote text-justify small\'>\
    In this demo, we will retrieve the 800bp sequence upstream of a list of genes from the organism <i>Saccharomyces cerevisiae</i>. Check the panel <b>Mandatory inputs</b> and then <b>Run analysis</b></blockquote>";
    
    demo_descr.innerHTML = descr_1;
    $("#organism_name").val("Saccharomyces cerevisiae");
    $("#organism").val("Saccharomyces_cerevisiae");
    $("#gene_selection").val(demo_genes);
    $("#from").val("-800");
    $("#to").val("-1");
    $("#noorf").removeAttr("checked");
}
</script>';

print ' <div class="col-lg-9 col-md-5 col-sm-8 col-xs-9 demo-buttons-container">

<button type="button" class="btn btn-info" onclick="setDemo1('. "'$demo_genes'" .')">DEMO</button> ';
print "</div>";

print $query->end_html;

exit(0);


