#!/usr/bin/env perl
############################################################
#
# $Id: RSAT_home.cgi,v 1.80 2013/10/09 07:04:55 jvanheld Exp $
#
# Time-stamp: <2003-10-22 11:53:22 jvanheld>
#
############################################################
#### this cgi script fills the HTML form for the program dna-pattern
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";
require RSAT::server;

$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

$query = new CGI;

&RSA_header("RSAT home","info");

#print $query->header;
#print $query->start_html(-class => "info",
#			 -author=>'Jacques.van-Helden\@univ-amu.fr',
#			 -style => {
#			   -src => ["main.css","font-awesome.min.css"],
#			   -type => 'text/css',
#			   -media => 'screen'
#			 });


print "<blockquote>";

&RSAT::server::DetectDeniedIP();

print <<EndText;

<table class="title" cellpadding="10" width="100%">
  <tr>

    <td align=center>
      <h1>Regulatory Sequence Analysis Tools</h1>
    </td>
    


  </TR>
</table>

 <div align="center">
	
	<br/>
	Welcome to <b>Regulatory Sequence Analysis Tools</b>
	(<b>RSAT</b>).  <p> <a href="RSAT_portal.html" target="_top"><img src="images/RSAT_logo_w150.png"  alt="RSAT"></a>
	  <br/>This web site provides a series of modular computer
	  programs specifically designed for the detection of regulatory
	  signals in non-coding sequences.</p>
	
	<p>RSAT servers have been up and running since 1997. The
	  project was initiated
	  by <a href="http://jacques.van-helden.perso.luminy.univ-amu.fr/"><b>Jacques
	  van Helden</b></a>, and is now pursued by
	  the <a href="people.html"><b>RSAT team</b></a>.
	</p>
	<p><i>This website is free and open to all users.</i></p>

<div class="hr-like"> </div>

<table class='simple' border='0'>
  <tr>
    <td>
    
    <div id="display" >
    <p ><span style="color: #cc6600;"><i  class="fa fa-question-circle fa-lg"></i><b> Which program to use ?</b></span> A guide to our main tools for new users. </p>

    <div class="holder">
    <b>1 - Choose your type of data to analyse</b></br>
        <select id="datatypes">
            <option value="">Choose Data type</option>
            <option value="data1">ChIP-seq</option>
            <option value="data2">List of gene names</option>
            <option value="data3">Sequences</option>
            <option value="data4">Matrices (PSSM)</option>
            <option value="data5">Coordinates (BED)</option>
            <option value="data6">List of variants</option>
        </select>
    <br/>
            </div>
        
        <div class="holder">
         <b>2 - Choose your biological question / analysis to perform </b></br>
            <select id="questions" disabled="true">
                <option value="">Choose selection</option>
                <option value="peak-motifs" class="data1">Which TF motifs are over-represented in this data set ?</option>
                
                <option value="retrieve-seq-programs" class="data2">I want to extract their promoter sequences (or other sequence features) </option>
                <option value="footprint-programs" class="data2">Which regulatory elements are conserved in promoters of orthologs ? (only for prokaryotes and fungi)</option>
                
                <option value="disco-programs" class="data3">Are there over-represented motifs in these sequences ? </option>
                <option value="scan-programs" class="data3">I want to scan these sequences with a motif </option>
                
                
                <option value="matrix-scan" class="data4">I want to scan sequences with these matrices </option>
                <option value="matrix-compa-programs" class="data4">I want to compare matrices (with known collections) </option>
                <option value="matrix-compa-programs" class="data4">I want to cluster and align matrices </option>
                <option value="convert-matrix" class="data4">I want to convert the matrix format </option>
                
                <option value="fetch-sequences" class="data5">I want to extract the sequences corresponding to these coordinates </option>
                
                 <option value="retrieve-variation-seq" class="data6">Obtain the variants and their flanking sequences </option>
                 <option value="scan-variations" class="data6">Which transcription factor binding sites are affected by these variants ? </option>
            </select>
        </div>
    
     <div class="holder">
         <b>3 - Relevant RSAT programs </b></br>
            <select id="tools" disabled="true">
                <option value="">Choose selection</option>
                <option value="peak-motifs_form.cgi" class="peak-motifs">peak-motifs</option>
                
                <option value="retrieve-seq_form.cgi" class="retrieve-seq-programs">retrieve sequences</option>
                <option value="retrieve-ensembl-seq_form.cgi" class="retrieve-seq-programs">retrieve Ensembl sequences (only for organisms supported at Ensembl.org)</option>
                
                 <option value="footprint-discovery_form.cgi" class="footprint-programs">footprint discovery</option>
                 <option value="footprint-scan_form.cgi" class="footprint-programs">footprint-scan</option>
                 
                 <option value="oligo-analysis_form.cgi" class="disco-programs">oligo-analysis (words)</option>
                 <option value="dyad-analysis_form.cgi" class="disco-programs">dyad-analysis (spaced pairs)</option>
                 
                 <option value="dna-pattern_form.cgi" class="scan-programs">dna-pattern</option>
                 <option value="matrix-scan-quick_form.cgi" class="scan-programs">matrix-scan (quick)</option>
                 <option value="matrix-scan_form.cgi" class="scan-programs">matrix-scan (full options)</option>
                 
                 <option value="matrix-scan-quick_form.cgi" class="matrix-scan">matrix-scan (quick)</option>
                 <option value="matrix-scan_form.cgi" class="matrix-scan">matrix-scan (full options)</option>
                 
                 <option value="compare-matrices_form.cgi" class="matrix-compa-programs">compare matrices</option>
                 <option value="matrix-clustering_form.cgi" class="matrix-compa-programs">matrix clustering</option>
                 
                 <option value="convert-matrix_form.cgi" class="convert-matrix">convert matrix</option>
                 
                 <option value="fetch-sequences_form.php" class="fetch-sequences">fetch sequences from UCSC (only for organisms supported at genome.ucsc.edu)</option>
                 
                 <option value="retrieve-variation-seq_form.cgi" class="retrieve-variation-seq">retrieve variation sequences</option>
                 <option value="variation-scan_form.cgi" class="scan-variations">scan variations</option>
            </select>
        </div>
        <br/>
             
        <i class="fa fa-hand-o-left fa-lg"></i> Complete list of online tools is in the left menu
    </div>
</div>

      
    </td>

    <td>
    <p><span style="color:red"><i class="fa fa-book fa-lg"></span></i> Check <span style="color:red"><b>latest RSAT paper</b></span> in <b><a target='_blank' href="http://nar.oxfordjournals.org/content/early/2015/04/21/nar.gkv362.full" target="_blank"> NAR web software issue 2015</a>
<pre>
	 </pre>
	<p><i class="fa fa-graduation-cap fa-lg"></i> Check <b>RSAT tutorial</b> at <b><a target='_blank' href="http://rsa-tools.github.io/tutorial_eccb14/index.html" target="tools">ECCB'14</a></b> and <a href="http://rsa-tools.github.io/teaching/index.html" target="tools"><b>all training material</b></a> 
	
	<p><i class="fa fa-graduation-cap fa-lg"></i> Learn how to use <b>Peak-motifs</b> with a <b>Nature Protocol</b> <a href='http://www.nature.com/nprot/journal/v7/n8/full/nprot.2012.088.html' target=_blank>[view article]</a></font>
		  
	<!--p><li> Latest features of RSAT presented in the <b>2011 NAR Web server issue</b> <br/> <a href='http://www.ncbi.nlm.nih.gov/pubmed/21715389' target=_blank>[Pubmed 21715389]</a></li-->
	
	<!--p><i class="fa fa-rss-square fa-lg"></i> Stay Tuned !! <a href="http://www.bigre.ulb.ac.be/forums/feed.php" target="_top"><b>RSS feed</b></a> to all RSAT news. -->
	  
	 <p><pre>
	 </pre><i class="fa fa-star fa-lg"></i> Try our <b>new programs</b> <img src="images/onebit_49.png" class="newprograms">
    </td>
  </tr>


</table>

</div>




<div class="hr-like"> </div>
<div align="center">
    <table class='serverlist' cellspacing='3' style="font-size:14px;">

  	<tr align=center valign=bottom>
  	  
  	  <td>
  	    <a href="http://rsat.france-bioinformatique.fr/fungi/" target="_top">
  	      <img src="images/logo_fungi.jpg" height='100' border='0'
  		   alt="fungi"></a>
  	    <br>hosted by <a target='_blank' href="https://www.france-bioinformatique.fr/">Institut Fran&Ccedil;ais de Bioinformatique (IFB), France </a>
  	  </td>
  
  	  <td>
  	    <a href="http://prokaryotes.rsat.eu/" target="_top">
  	      <img src="images/logo_prokaryotes.jpg" height='100' border='0' alt="prokaryotes"></a>
  	   <br>hosted by <a target='_blank' href="https://www.ccg.unam.mx/genomica-computacional/">Centro de Ciencias Genomicas (CCG), Cuernavaca, Mexico</a>
  	  </td>
  
  	<tr align=center valign=bottom>
  	  
  	  <td>
  	    <a href="http://rsat.france-bioinformatique.fr/fungi/" target="_top">
  	      <img src="images/logo_fungi.jpg" height='100' border='0'
  		   alt="fungi"></a>
  	    <br>hosted by <a target='_blank' href="https://www.france-bioinformatique.fr/">Institut Fran&Ccedil;ais de Bioinformatique (IFB), <br>Paris, France</a>
  	  </td>
  
  	  <td>
  	    <a href="http://prokaryotes.rsat.eu/" target="_top">
  	      <img src="images/logo_prokaryotes.jpg" height='100' border='0' alt="prokaryotes"></a>
  	   <br>hosted by <a target='_blank' href="https://www.ccg.unam.mx/genomica-computacional/">Centro de Ciencias Genomicas (CCG), <br>Cuernavaca, Mexico</a>
  	  </td>
  
      <td>
  	    <a href="http://rsat.france-bioinformatique.fr/metazoa/" target="_top">
  	      <img src="images/logo_metazoa.jpg" height='100' border='0' alt="metazoa"></a>
  	   <br>hosted by <a target='_blank' href="https://www.france-bioinformatique.fr/">Institut Fran&Ccedil;ais de Bioinformatique (IFB), <br>Paris, France</a>
  	  </td>
  
  	</tr>
  	<tr align=center valign=bottom>
  
  	  <td>
  	    <a href="http://protists.rsat.eu/" target="_top">
  	      <img src="images/logo_protists.jpg" height='100' border='0' alt="protists"></a>
  	    <br>hosted by <a target='_blank' href="https://www.france-bioinformatique.fr/">Institut Fran&Ccedil;ais de Bioinformatique (IFB), <br>Paris, France</a>
  	  </td>
  
  	  <td>
  	    <a href="http://plants.rsat.eu/" target="_top">
  	      <img src="images/logo_plants.jpg" height='100' border='0' alt="plants"></a>
  	    <br>hosted by <a target='_blank' href="https://eead-csic-compbio.github.io/">Lab Computational & structural biology, EEAD-CSIC, <br>Zaraagoza, Spain</a>
  	  </td>
  
      <td>
  	    <a href="http://teaching.rsat.eu" target="_top" >
  	      <img src="images/logo_teaching.jpg" height='100' border='0' alt="teaching"></a>
  	    <br>hosted by <a target='_blank' href="https://www.france-bioinformatique.fr/">Institut Fran&Ccedil;ais de Bioinformatique (IFB), <br>Paris, France</a>
  	  </td>
  	  
  	</tr>

    
    </table>
</div>

<p/>
    
EndText

@orgs =  &RSAT::OrganismManager::get_supported_organisms_web();
print "<div align='center'>";
print "<i class='fa fa-bar-chart  fa-lg'></i> <b>", scalar(@orgs), "</b> <i>organisms supported";
print " on $ENV{rsat_site} (<a href='$ENV{rsat_www}' target=_top>",$ENV{rsat_www},"</a>)</i>\n";
# print &ListSupportedOrganisms("html_list");
if ($group_specificity = $ENV{group_specificity}) {
   print ". Group specificity: <b>",$group_specificity,"</b></p>";
}
print "</div>\n";

&UpdateLogFile();
$count = &UpdateCounterFile();


print <<EndAddress;
<div class="hr-like"> </div>
<p/>
   <span style="color: #cc6600;"><i class="fa fa-pencil fa-lg"></i> <b>Citing RSAT complete suite of tools:</b></span>
	<div align="left">
      <ul>
<li><font color="red"><b>New !</b></font>  Medina-Rivera A, Defrance M, Sand O, Herrmann C, Castro-Mondragon J, Delerce J, Jaeger S, Blanchet C, Vincens P, Caron C, Staines DM, Contreras-Moreira B, Artufel M, Charbonnier-Khamvongsa L, Hernandez C, Thieffry D, Thomas-Chollier M, van Helden J (2015)
	  <b>RSAT 2015: Regulatory Sequence Analysis Tools </b>. Nucleic Acids Res. 2015 (Web Server issue) in press.
	  <!--a href=http://www.ncbi.nlm.nih.gov/pubmed/21715389>[Pubmed 21715389]</a-->
	  <a href='http://nar.oxfordjournals.org/content/early/2015/04/21/nar.gkv362.full'>[Full text]</a>
	</li>
	<li>Thomas-Chollier M, Defrance M, Medina-Rivera A, Sand O, Herrmann C, Thieffry D, van Helden J. (2011)
	  <b>RSAT 2011: regulatory sequence analysis tools</b>. Nucleic Acids Res. 2011 Jul;39(Web Server issue):W86-91.
	  <a href=http://www.ncbi.nlm.nih.gov/pubmed/21715389>[Pubmed 21715389]</a>
	  <a href='http://nar.oxfordjournals.org/content/39/suppl_2/W86.long'>[Full text]</a>
	</li>
	<li>Thomas-Chollier, M., Sand, O., Turatsinze, J. V., Janky, R.,
	  Defrance, M., Vervisch, E., Brohee, S. & van Helden, J. (2008). <b>RSAT:
	    regulatory sequence analysis tools</b>. Nucleic Acids
	  Res.
	  <a href=http://www.ncbi.nlm.nih.gov/pubmed/18495751>[Pubmed 18495751]</a>
	  <a href='http://nar.oxfordjournals.org/cgi/content/full/36/suppl_2/W119'>[Full text]</a>
	</li>
	<li>van Helden, J. (2003). <b>Regulatory sequence analysis
	    tools</b>. Nucleic Acids Res. 2003 Jul 1;31(13):3593-6. [<a target=_blank
									href="http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=12824373&dopt=Abstract">Pubmed
	    12824373</a>] [<a target=_blank
			      href=http://nar.oupjournals.org/cgi/content/full/31/13/3593>Full
	    text</a>] [<a target=_blank
			  href=http://nar.oupjournals.org/cgi/reprint/31/13/3593.pdf>pdf</a>]
	</li>
      </ul>
</div>

<i class="fa fa-pencil fa-lg"></i> For citing <b>individual tools</b>: the reference of each tool is indicated on top of their query form.</br>
<i class="fa fa-pencil fa-lg"></i> <a href="publications.html">All RSAT publications</a> and all publications for <a href=neat_publications.html> the Network Analysis Tools</a>. 


</p>

      

<hr  class = 'portal'/>
<h4 class='footer'>
      <div align="center">RSAT logos designed by Mauricio Guzman (<a target="_blank" href="http://www.altamirastudio.com.mx/">http://www.altamirastudio.com.mx/</a>)</div>

<i>
For suggestions or information request, please contact :
<br>
<a href="mailto:Jacques.van-Helden\@univ-amu.fr (Jacques van Helden)">
Jacques van Helden (Jacques.van-Helden\@univ-amu.fr)</i>
</a>
</h4>
EndAddress

my $version_git = &RSAT::server::GetGitLastCommitDate();

print "<div align=center> Running with RSAT code from : ".$version_git."</div>";

print "</blockquote>\n";
&google_analytics_tag();

print "<script type=\"text/javascript\">";

 print  `cat $ENV{RSAT}/perl-scripts/lib/js/DataTables-1.10.4/media/js/jquery.js`;

print "</script>";

print "<script type=\"text/javascript\">";

 print  `cat $ENV{RSAT}/perl-scripts/lib/js/chained/jquery.chained.min.js`;

print "</script>";


#print '
#<script src="http://www.appelsiini.net/download/jquery.chained.js" type="text/javascript" charset="utf-8"></script> ';
 
print '<script type="text/javascript" charset="utf-8"> 
$(function(){ 
	$("#questions").chained("#datatypes"); 
	$("#tools").chained("#questions"); 
	
	document.getElementById("tools").onchange = function() {
        if (this.selectedIndex!==0) {
            window.location.href = this.value;
        }        
    };
	
	}); 
	
	
</script> ';



print $query->end_html, "\n";

exit(0);

