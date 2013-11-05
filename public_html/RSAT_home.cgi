#!/usr/bin/perl
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
print $query->header;
print $query->start_html(-class => "info",
			 -author=>'jacques.van.helden@ulb.ac.be',
			 -style => { 	-src => "$ENV{rsat_www}/main.css",
                             	       	-type => 'text/css',
                             		-media => 'screen' });

print "<blockquote>";

&RSAT::server::DetectDeniedIP();

print <<EndText;

<table class="title" cellpadding="10" width="100%">
  <tr>
    <td align=center valign = top>
      <a href="http://www.bigre.ulb.ac.be" target=_blank>
	<img src="images/bigre_logo.png" alt="BiGRE lab" border=1 height=75>
	<br>
	<font color="#0000dd" size=-2>
	  BiGRe</a> - <a href="http://www.ulb.ac.be" target=_blank>ULB</a>
      </font>
    </td>
    
    <td align=center>
      <h1>Regulatory Sequence Analysis Tools</h1>
    </td>
    
    
    <td align=center valign = top width="160">
      <a target=_blank href="http://embnet.ccg.unam.mx/">
	<img src="images/ccg_logo_lr.jpg" alt="CCG" border=0  width=75></a>
      <br>
      <font color="#006600" size=-2>Laboratorio de Biologia Computacional<br>UNAM/CCG</font>
    </td>
    
  </TR>
</table>

<div class = 'serverlist'>
  
  <table border="0" cellspacing="3" cellpadding="7">
    <TR>
    <td><a href="tool_map.html"><b>
    Tool Map</b></a></td>

    <td><a href="intro.html"><b>
    Introduction</b></a></td>

    <td><a href="http://www.bigre.ulb.ac.be/forums/"><b>
    Forum</b></a></td>

    <td><a href="tutorials/tutorials.html"><b>
    Tutorials</b></a></td>

    <td><a href="publications.html"><b>
    Publications</b></a></td>

    <td><a href="credits.html"><b>
    Credits</b></a></td>

    <td><a href="people.html"><b>
    People</b></a></td>

    <td><a href="data/"><b>
    Data</b></a></td>

<!--    <td><a href="links.html"><b>
    Links</b></a></td>
-->

    <td><a href="distrib/index.html"><b>
    Download</b></a></td>
    
    <td><a href="http://www.bigre.ulb.ac.be/forums/"><b>
    [Forum]</b></a></td>

    </tr>
  </table>
</div>


<table border=0 cellpadding=3>
  <tr>
    <td COLSPAN=3>
      <p>Welcome to <b>Regulatory Sequence Analysis Tools</b>
	(<b>RSAT</b>). This web site provides a series of modular computer
	programs specifically designed for the detection of regulatory
	signals in non-coding sequences.</p>
    </td>
  </tr>
  <tr>
    <td  WIDTH="50%" bgcolor=#F6E6CA>
      <ul> 
	<li><font color="red"><b>New !</b></font> Learn how to use <b>Peak-motifs</b> with a <b>Nature Protocol</b> <a href='http://www.nature.com/nprot/journal/v7/n8/full/nprot.2012.088.html' target=_blank>[view article]</a></font></li>
	
	<p><li> <b>Peak-motifs</b> is now published in <b>NAR</b> <a href='http://nar.oxfordjournals.org/content/early/2011/12/08/nar.gkr1104.full' target=_blank>[view article]</a></li>
	  
	<p><li> Latest features of RSAT presented in the <b>2011 NAR Web server issue</b> <br/> <a href='http://www.ncbi.nlm.nih.gov/pubmed/21715389' target=_blank>[Pubmed 21715389]</a></li>
      </ul>
    </td>

    <td WIDTH="30%" bgcolor=#F6E6CA>
      <ul>
      <li>Try our <b>new programs</b> <img src="images/onebit_49.png" class="newprograms"></li> <p>
	<li> Check the <b>latest news</b> in <a
						href="http://www.bigre.ulb.ac.be/forums/viewforum.php?f=25&sid=a2f000bd4cc6f7f8260ae9547db9a72d"><b>our [Forum] </b></a></li>
	
	<p><li> Stay Tuned !! <b>RSS feed</b> to all RSAT news <a href="http://www.bigre.ulb.ac.be/forums/feed.php" target="_top"><IMG class="rss" SRC="images/feed.png" BORDER=0></a></li>
	  
	  
	<p><li> <b><a href="citing_rsat.html"> How to cite RSAT ? </a> </b></li></p>
      </ul>
    </td>
    <td WIDTH="*"></td>
  </tr>

  <tr>
    <td COLSPAN=3>
      <p> This website is free and open to all users.
      </p>
    </td>
  </tr>
</table>

<hr  align='left'></hr>
<div id = 'serverlist'>
  <table class='serverlist'>
    
    
    <tr>
      <td align=center colspan=3><h3>Regulatory Sequence Analysis Tools - Web servers</h3></td>
    </tr>
    <tr align=center valign=bottom>
      
      <td>
	<a href="http://rsat.ulb.ac.be/rsat/" target="_top">
	  <b>Brussels - Belgium</b><br>
	  <img src="images/manneken_pis.jpg" height=80 border=0><br>
	  http://rsat.ulb.ac.be/rsat/</a>
      </td>

<!--
      <td>
	<a href="http://wwwsup.scmbb.ulb.ac.be/rsat/" target="_top">
	  <b>Brussels (2) - Belgium</b><br>
	  <img src="images/Atomium_icon_80.png" height=80 border=0><br>
	  http://wwwsup.scmbb.ulb.ac.be/rsat/</a>
      </td>
-->


      <td align=center>
	<a href='http://www.rsat.fr/' target='_top'>
	  <b>Marseille TAGC - France</b><br>
	  <img src="images/marseille_calanques.jpg" height=80 border=0><br>
	  http://www.rsat.fr/</a>
      </td>

      <td ALIGN=CENTER>
	<a href="http://embnet.ccg.unam.mx/rsa-tools/" target="_top">
	  <b>Cuernavaca - Mexico</b><br>
	  <img src="images/zapata.jpg" height=80 border=0><br>
	  http://embnet.ccg.unam.mx/rsa-tools/</a>
      </td>

    </tr><tr>

      <td align=center>
	<a href='http://rsat01.biologie.ens.fr/rsa-tools/' target='_top'>
	  <b>ENS Paris - France</b><br>
	  <img src="images/paris.jpg" height=80 border=0><br>
	  http://rsat01.biologie.ens.fr/rsa-tools/</a>
      </td>


      <td align=center>
	<a href="http://bongcam1.hgen.slu.se/rsat/" target="_top">
	  <!--  <a href="http://liv.hgen.slu.se/rsa-tools/" target="_top">-->
	  <b>Uppsala - Sweden</b><br>
	  <img src="images/uppsala_lcb.jpg" height=80 border=0><br>
	  http://bongcam1.hgen.slu.se/rsat/</a>
	<!--http://liv.hgen.slu.se/rsa-tools/</a>-->
      </td>

      <td align=center>
	<a href="http://lmmge-webserver.unisa.it/rsat/" target="_top">
	  <b>Salerno - Italia</b><br>
	  <img src="images/salerno_icon.jpg" height=80 border=0><br>
	  http://lmmge-webserver.unisa.it/rsat//</a>
      </td>


    </tr><tr>
      <td align=center>
	<a href="http://anjie.bi.up.ac.za/rsa-tools/" target="_top">
	  <b>Pretoria - South-Africa</b><br>
	  <img src="images/pretoria_icon.jpg" height=80 border=0><br>
	  http://anjie.bi.up.ac.za/rsa-tools/</a>
      </td>
      
      <td>
	<a href="http://rsat.sb-roscoff.fr/" target="_top">
	  <b>Roscoff - France</b><br>
	  <img src="images/roscoff_sb.jpg" height=80 border=0><br>
	  http://rsat.sb-roscoff.fr/</a>
      </td>


      <!--
	  <td align=center>
	    <a href="http://rsat.ccb.sickkids.ca/" target="_top">
	      <b>Toronto - Canada</b><br>
	      <img src="images/toronto.jpg" height=80 border=0><br>
	      http://rsat.ccb.sickkids.ca/</a>
	  </td>
      </tr><tr>
      
      
  <td align=center>
    <a href="http://af.boku.ac.at:4080/rsa-tools/" target="_top">
      <b>Vienna - Austria</b><br>
      <img src="http://www.wien.gv.at/english/cityhall/images/cityhall.jpg" height=80 border=0><br>
      http://af.boku.ac.at:4080/rsa-tools/</a>
  </td>
  
  <td></td>
  
	  </tr><tr>
	  
  <td ALIGN=CENTER>
    <a href="http://www.flychip.org.uk/rsa-tools/" target="_top">
      <b>Cambridge - UK</b><br>
      <img src="images/cambridge.jpg" height=80 border=0><br>
      http://www.flychip.org.uk/rsa-tools/</a>
  </td>
  -->

    </tr>
  </table>
</div>
EndText

@orgs =  &RSAT::OrganismManager::get_supported_organisms();
print "<h4 align ='center'>", scalar(@orgs) ," organisms supported on <a href='$ENV{rsat_www}' target=_top>",$ENV{rsat_www},"</a></h4>\n";
# print &ListSupportedOrganisms("html_list");

&UpdateLogFile();
$count = &UpdateCounterFile();

print <<EndAddress;
<p>
<hr  align = 'left'>
<h4 class='footer'><i>
For suggestions or information request, please contact :
<br>
<a href="mailto:jvanheld\@bigre.ulb.ac.be (Jacques van Helden)">
Jacques van Helden (jvanheld\@bigre.ulb.ac.be)</i>
</a>
</h4>
EndAddress

print "</blockquote>\n";
&google_analytics_tag();
print $query->end_html, "\n";

exit(0);

