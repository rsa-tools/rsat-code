#!/usr/bin/perl
############################################################
#
# $Id: RSAT_home.cgi,v 1.61 2011/04/05 14:44:32 jvanheld Exp $
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
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";


$query = new CGI;
print $query->header;
print $query->start_html(-class => "info",
			 -author=>'jacques.van.helden@ulb.ac.be',
			 -style => { 	-src => "$ENV{rsat_www}/main.css",
                             	       	-type => 'text/css',
                             		-media => 'screen' });

print "<blockquote>";

print <<EndText;

<table class = 'title' cellpadding=10 width = 100%>
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
    <font color="#006600" size=-2>Laboratorio de Biologia Computacional<br>UNAM/CCG
    </FONT>
    </td>

    </TR>
    </TABLE>
    <div class = 'serverlist'>

    <table border=0 cellspacing=3 cellpadding=7>
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

    </tr>
    </table>
    </div>

    <table  cellpadding=10>
    <tr >
    <td>

    <p>Welcome to <b>Regulatory Sequence Analysis Tools</b>
    (<b>RSAT</b>). This web site provides a series of modular computer
    programs specifically designed for the detection of regulatory
    signals in non-coding sequences.</p>

    <ul>
	<p><li>Try our <b>new programs</b> <img src="images/onebit_49.png" class="newprograms"></li>

     <p><li><font color="red">New !</font> <a target='_top' href='http://rsat.ulb.ac.be/rsat/tutorial_ECCB_2010/'>Supporting
     material</a> of the 3rd Tutorial presented at </a> of <a target='_blank' href='http://www.eccb10.org/'>ECCB
     2010</a> (Sunday Sept 26).</li></p>

     <p><li> <font color='red'>Check the latest news in </font> <a
    href="http://www.bigre.ulb.ac.be/forums/viewforum.php?f=25&sid=a2f000bd4cc6f7f8260ae9547db9a72d"><b>our forum </b></a></li></p>

	<p><li> Stay Tuned !! RSS feed to all RSAT news <a href="http://www.bigre.ulb.ac.be/forums/feed.php" target="_top"><IMG class="rss" SRC="images/feed.png" BORDER=0></a></li>


    <li> <b><a href="citing_rsat.html"> How to cite RSAT ? </a> </b></li>
	</ul>

    <p> This website is free and open to all users.
   </p>
</td>


    <td>
    <table class='title' cellpadding=0 cellspacing=0><tr><td align = 'center'>
    <h3>Warnings<br><br>

    <a href=warnings.html>Vertebrate genomes</a>

    </h3>
    </table>
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
  <IMG SRC="images/manneken_pis.jpg" HEIGHT=80 BORDER=0><br>
  http://rsat.ulb.ac.be/rsat/</a>
  </td>

  <td>
  <a href="http://wwwsup.scmbb.ulb.ac.be/rsat/" target="_top">
  <b>Brussels (2) - Belgium</b><br>
  <IMG SRC="images/Atomium_icon_80.png" HEIGHT=80 BORDER=0><br>
  http://wwwsup.scmbb.ulb.ac.be/rsat/</a>
  </td>

  <td ALIGN=CENTER>
  <a href="http://embnet.ccg.unam.mx/rsa-tools/" target="_top">
  <b>Cuernavaca - Mexico</b><br>
  <IMG SRC="images/zapata.jpg" HEIGHT=80 BORDER=0><br>
  http://embnet.ccg.unam.mx/rsa-tools/</a>
  </td>

</tr>
<tr>

  <td align=center>
  <a href="http://liv.bmc.uu.se/rsa-tools/" target="_top">
  <b>Uppsala - Sweden</b><br>
  <IMG SRC="images/uppsala_lcb.jpg" HEIGHT=80 BORDER=0><br>
  http://liv.bmc.uu.se/rsa-tools/</a>
  </td>

  <td align=center>
  <a href='http://tagc.univ-mrs.fr/rsa-tools/' target='_top'>
  <b>Marseille TAGC - France</b><br>
  <IMG SRC="images/calanques.jpg" HEIGHT=80 BORDER=0><br>
  http://tagc.univ-mrs.fr/rsa-tools/</a>
<br><font size=-2>(photo by <a target=_blank href=http://www.lim.univ-mrs.fr/~guenoche/Walk1.html>Alain Gu&eacute;noche</a>)</font>
  </td>

  <td align=center>
  <a href='http://rsat01.biologie.ens.fr/rsa-tools/' target='_top'>
  <b>ENS Paris - France</b><br>
  <IMG SRC="images/paris.jpg" height=80 border=0><br>
  http://rsat01.biologie.ens.fr/rsa-tools/</a></font>
  </td>

  </tr><tr>

  <td align=center>
  <a href="http://anjie.bi.up.ac.za/rsa-tools/" target="_top">
  <b>Pretoria - South-Africa</b><br>
  <img src="images/pretoria_icon.jpg" height=80 border=0><br>
  http://anjie.bi.up.ac.za/rsa-tools/</a>
  </td>

<!--
  <td align=center>
  <a href="http://rsat.ccb.sickkids.ca/" target="_top">
  <b>Toronto - Canada</b><br>
  <IMG SRC="images/toronto.jpg" HEIGHT=80 BORDER=0><br>
  http://rsat.ccb.sickkids.ca/</a>
  </td>
</tr><tr>


  <td align=center>
  <a href="http://af.boku.ac.at:4080/rsa-tools/" target="_top">
  <b>Vienna - Austria</b><br>
  <IMG SRC="http://www.wien.gv.at/english/cityhall/images/cityhall.jpg" HEIGHT=80 BORDER=0><br>
  http://af.boku.ac.at:4080/rsa-tools/</a>
  </td>

  <td></td>

-->


<!--
</tr><tr>

  <td ALIGN=CENTER>
  <a href="http://www.flychip.org.uk/rsa-tools/" target="_top">
  <b>Cambridge - UK</b><br>
  <IMG SRC="images/cambridge.jpg" HEIGHT=80 BORDER=0><br>
  http://www.flychip.org.uk/rsa-tools/</a>
  </td>
-->

    </tr>
    </TABLE>
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

