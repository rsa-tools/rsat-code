#!/usr/bin/perl
############################################################
#
# $Id: RSAT_home.cgi,v 1.28 2008/02/12 20:19:49 rsat Exp $
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
			   -author=>'jvanheld@scmbb.ulb.ac.be',
			   
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
    <H1>Regulatory Sequence Analysis Tools</H1>
    </td>
    
    
    <td align=center valign = top width="160">
    <a target=_blank href="http://embnet.ccg.unam.mx/">
    <img src="images/lablogo.gif" alt="lab logo" border=0 height=48 width=75></a>
    <br>
    <font color="#006600" size=-2>LABORATORIO DE BIOLOGIA COMPUTACIONAL
    </FONT>
    </TD>
    
    </TR>
    </TABLE>
    <div class = 'serverlist'>
    
    <table border=0 cellspacing=3 cellpadding=7>
    <TR>
    <TD><A HREF="tool_map.html"><B>
    Tool Map</B></A></TD>
    
    <TD><A HREF="intro.html"><B>
    Introduction</B></A></TD>
    
    <TD><A HREF="FAQ.html"><B>
    FAQ</B></A></TD>
    
    <TD><A HREF="tutorials/tutorials.html"><B>
    Tutorials</B></A></TD>
    
    <TD><A HREF="publications.html"><B>
    Publications</B></A></TD>
    
    <TD><A HREF="credits.html"><B>
    Credits</B></A></TD>
    
    <TD><A HREF="data/"><B>
    Data</B></A></TD>
    
    <TD><A HREF="change_history.html"><B>
    Change history</B></A></TD>
    
    <TD><A HREF="links.html"><B>
    Links</B></A></TD>
    
    </tr>
    </table>
    </div>
    
    <table  cellpadding=10> 
    <tr >
    <td>
    <P>
    Welcome to <b>Regulatory Sequence Analysis Tools</b>
    (<B>RSAT</B>). This web site provides a series of modular computer
    programs specifically designed for the detection of regulatory
    signals in non-coding sequences.

    <P> <font color='red'>New !</font> The <b>RSAT lookup</b> has been
      reshaped by <a target=_blank
      href=http://www.scmbb.ulb.ac.be/Users/morgane/>Morgane
      Thomas-Chollier</a> and <a target=_blank
      href=http://www.scmbb.ulb.ac.be/Users/sylvain/>Sylvain
      Broh&eacute;e</a>.

    </p>

    <P>
      <font color='red'>New !</font> For bioinformaticians, RSAT tools
      are now also available as <a href=web_services.html><b>Web Services</b></a>. 
    </p>

    <P> This website is free and open to all users.
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
  <A HREF="http://rsat.scmbb.ulb.ac.be/rsat/" target="_top">
  <B>Brussels - Belgium</B><BR>
  <IMG SRC="images/manneken_pis.jpg" HEIGHT=80 BORDER=0><BR>
  http://rsat.scmbb.ulb.ac.be/rsat/</A>
  </TD>

  <TD ALIGN=CENTER>
  <A HREF="http://embnet.ccg.unam.mx/rsa-tools/" target="_top">
  <B>Cuernavaca - Mexico</B><BR>
  <IMG SRC="images/zapata.jpg" HEIGHT=80 BORDER=0><BR>
  http://embnet.ccg.unam.mx/rsa-tools/</A>
  </TD>
  
  <TD ALIGN=CENTER>
  <A HREF="http://liv.bmc.uu.se/rsa-tools/" target="_top">
  <B>Uppsala - Sweden</B><BR>
  <IMG SRC="images/uppsala_lcb.jpg" HEIGHT=80 BORDER=0><BR>
  http://liv.bmc.uu.se/rsa-tools/</A>
  </TD>

</tr>
<tr>

  <td align=center>
  <A HREF="http://crfb.univ-mrs.fr/rsaTools/" target="_top">
  <B>Marseille - France</B><BR>
  <IMG SRC="images/calanques.jpg" HEIGHT=80 BORDER=0><BR>
  http://crfb.univ-mrs.fr/rsaTools/</a>
<br><font size=-2>(photo by <a target=_blank href=http://www.lim.univ-mrs.fr/~guenoche/Walk1.html>Alain Gu&eacute;noche</a>)</font>
  </TD>

  <td align=center>
  <A HREF="http://rsat.ccb.sickkids.ca/" target="_top">
  <B>Toronto - Canada</B><BR>
  <IMG SRC="images/toronto.jpg" HEIGHT=80 BORDER=0><BR>
  http://rsat.ccb.sickkids.ca/</A>
  </TD>

  <td align=center>
  <a href="http://www.bi.up.ac.za/rsa-tools/" target="_top">
  <B>Pretoria - South-Africa</b><br>
  <img src="images/pretoria_icon.jpg" height=80 border=0><br>
  http://www.bi.up.ac.za/rsa-tools/</a>
  </td>


<!--
</tr><tr>

  <TD ALIGN=CENTER>
  <A HREF="http://www.flychip.org.uk/rsa-tools/" target="_top">
  <B>Cambridge - UK</B><BR>
  <IMG SRC="images/cambridge.jpg" HEIGHT=80 BORDER=0><BR>
  http://www.flychip.org.uk/rsa-tools/</A>
  </TD>
-->

    </tr>
    </TABLE>
    </div>
EndText

@orgs =  &ListSupportedOrganisms("keys");
print "<H4 align ='center'>", scalar(@orgs) ," organisms supported on <A HREF='$ENV{rsat_www}' target=_top>",$ENV{rsat_www},"</A></H4>\n";
# print &ListSupportedOrganisms("html_list");

&UpdateLogFile();
$count = &UpdateCounterFile();

print <<EndAddress;
<P>
<HR  align = 'left'>
<H4 class='footer'><i>
For suggestions or information request, please contact :
<BR>
<A HREF="mailto:jvanheld\@scmbb.ulb.ac.be (Jacques van Helden)">
Jacques van Helden (jvanheld\@scmbb.ulb.ac.be)</i>
</A>
</H4>
EndAddress

print "</blockquote>";
print $query->end_html, "\n";

exit(0);

