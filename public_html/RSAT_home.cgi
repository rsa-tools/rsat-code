#!/usr/bin/perl
############################################################
#
# $Id: RSAT_home.cgi,v 1.9 2005/11/25 07:01:16 rsat Exp $
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
require "RSA.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";


$query = new CGI;
print $query->header;
print $query->start_html(-title=>'Regulatory Sequence Analysis Tools',
			 -author=>'jvanheld@scmbb.ulb.ac.be',
			 -BGCOLOR=>'#FFEEDD');

print <<EndText;

  
  <table cellpadding=10 align=center>
  <tr>

  
  <td valign=top width="160">
  <center>
  <img src="images/ULB_Blue_215.gif" alt="ULB_Blue_215.gif" border=0 height=75 width=75>
  <br>
  <font color="#0000dd" size=-2>
  <a href="http://www.scmbb.ulb.ac.be" target=_blank>SCMBB</a> - <a href="http://www.ulb.ac.be" target=_blank>ULB</a>
  </font>
  </center>
  </td>

  <td align=center valign=center>
  <font face="Helvetica,Arial" SIZE=+1><B>Regulatory Sequence Analysis Tools</B></font>
  </td>
  
  
  <TD VALIGN=TOP WIDTH="160">
  <CENTER>
  <A HREF="http://www.ccg.unam.mx/Computational_Genomics/"
  target=_blank> <IMG SRC="images/lablogo.gif" ALT="lab logo" BORDER=0
  HEIGHT=48 WIDTH=75></A>
  <BR>
  <FONT COLOR="#006600" SIZE=-2>
  LABORATORIO DE BIOLOGIA COMPUTACIONAL
  </FONT>
  </CENTER>
  </TD>
  
  </TR>
  </TABLE>
  
  <HR SIZE=4 WIDTH="100%">
  
  <TABLE BORDER=0 CELLSPACING=3 CELLPADDING=7 ALIGN=CENTER BGCOLOR="#FFDDCC">
  <TR>
  <TD><FONT FACE="Helvetica" SIZE=-1><A HREF="intro.html"><B>
  Introduction</B></A></FONT></TD>
  
  <TD><FONT FACE="Helvetica" SIZE=-1><A HREF="FAQ.html"><B>
  FAQ</B></A></FONT></TD>
  
  <TD><FONT FACE="Helvetica" SIZE=-1><A HREF="tutorials/tutorials.html"><B>
  Tutorials</B></A></FONT></TD>
  
  <TD><FONT FACE="Helvetica" SIZE=-1><A HREF="publications.html"><B>
  Publications</B></A></FONT></TD>
  
  <TD><FONT FACE="Helvetica" SIZE=-1><A HREF="credits.html"><B>
  Credits</B></A></FONT></TD>
  
  <TD><FONT FACE="Helvetica" SIZE=-1><A HREF="data/"><B>
  Data</B></A></FONT></TD>
  
  <TD><FONT FACE="Helvetica" SIZE=-1><A HREF="change_history.html"><B>
  Change history</B></A></FONT></TD>
  
  <TD><FONT FACE="Helvetica" SIZE=-1><A HREF="links.html"><B>
  Links</B></A></FONT></TD>
  
  </TR>
  </TABLE>

<font face=Arial,Helvetica size=-1>  
  
<table> 
<tr  valign=top>
<td>
  <P>
  Welcome to Regulatory Sequence Analysis Tools (<B>RSAT</B>). This
  site provides a series of modular computer programs specifically
  designed for the detection of regulatory signals in non-coding
  sequences.  <P>

  <P>
  This web site is freely available for academic users. For users from
  commercial companies, please read our <A
  HREF="disclaimer.html">disclaimer</A>.
    </font>
</td>


<td>
<table  border=1 bgcolor='#ffcccc' cellpadding=5 cellspacing=5 border=3><tr><td>
<B>
<font face=arial,helvetica color='#ff0000'>Warnings</font>
<ul>
<li><a href=warnings.html>Mammalian genomes</a>
</ul>
</b>
</table>
</td>

</tr>
</table>


  <P>
  <h2 align=center>Regulatory Sequence Analysis Tools - Web servers</h2>
  
<table border=5 cellspacing=5 cellpadding=5 align=center>

<tr align=center valign=bottom>

  <td colspan=3>
  <font size=-1>
  <A HREF="http://rsat.scmbb.ulb.ac.be/rsat/" target="_top">
  <B>Brussels - Belgium</B><BR>
  <IMG SRC="images/manneken_pis.jpg" HEIGHT=80 BORDER=0><BR>
  http://rsat.scmbb.ulb.ac.be/rsat/</A>
  </FONT>
  </TD>

</tr></tr>


  <TD ALIGN=CENTER>
  <font size=-1>
  <A HREF="http://embnet.cifn.unam.mx/rsa-tools/" target="_top">
  <B>Cuernavaca - Mexico</B><BR>
  <IMG SRC="images/zapata.jpg" HEIGHT=80 BORDER=0><BR>
  http://embnet.cifn.unam.mx/rsa-tools/</A>
  </FONT>
  </TD>
  
  <TD ALIGN=CENTER>
  <FONT SIZE=-1>
  <A HREF="http://liv.bmc.uu.se/rsa-tools/" target="_top">
  <B>Uppsala - Sweden</B><BR>
  <IMG SRC="images/uppsala_lcb.jpg" HEIGHT=80 BORDER=0><BR>
  http://liv.bmc.uu.se/rsa-tools/</A>
  </FONT>
  </TD>

  <td align=center>
  <font size=-1>
  <A HREF="http://gin.univ-mrs.fr/~jvanheld/rsa-tools/" target="_top">
  <B>Marseille - France</B><BR>
  <IMG SRC="images/calanques.jpg" HEIGHT=80 BORDER=0><BR>
  http://gin.uiv-mrs.fr/~jvanheld/rsa-tools/</A>
  </FONT>
<br><font size=-2>(photo by <a target=_blank href=http://www.lim.univ-mrs.fr/~guenoche/Walk1.html>Alain Gu&eacute;noche</a>)</font>
  </TD>

</tr></tr>


  <TD ALIGN=CENTER>
  <font size=-1>
  <A HREF="http://www.flychip.org.uk/rsa-tools/" target="_top">
  <B>Cambridge - UK</B><BR>
  <IMG SRC="images/cambridge.jpg" HEIGHT=80 BORDER=0><BR>
  http://www.flychip.org.uk/rsa-tools/</A>
  </FONT>
  </TD>

  <td align=center>
  <font size=-1>
  <A HREF="http://rsat.ccb.sickkids.ca/" target="_top">
  <B>Toronto - Canada</B><BR>
  <IMG SRC="images/toronto.jpg" HEIGHT=80 BORDER=0><BR>
  http://rsat.ccb.sickkids.ca/</A>
  </FONT>
  </TD>

  <td align=center>
  <font size=-1>
  <a href="http://www.bi.up.ac.za/rsa-tools/" target="_top">
  <B>Pretoria - South-Africa</b><br>
  <img src="images/pretoria_icon.jpg" height=80 border=0><br>
  http://www.bi.up.ac.za/rsa-tools/</a>
  </font>
  </td>


</tr>

  </TABLE>
  
EndText

@orgs =  &ListSupportedOrganisms("keys");

print "<H4>", scalar(@orgs) ," organisms supported on <A HREF='$WWW_RSA' target=_top>",$WWW_RSA,"</A></H4>\n";
print &ListSupportedOrganisms("html_list");

&UpdateLogFile();
$count = &UpdateCounterFile();

# print ("This page has been visited\n<B>", $count,
#        "</B>\ntimes since May 1998.\n",
#        "Enjoy your visit !");


print <<EndAddress;
<P>
<HR WIDTH="100%">
<CENTER>
<FONT SIZE=-1>
For suggestions or information request, please contact :
<BR>
<A HREF="mailto:jvanheld\@scmbb.ulb.ac.be (Jacques van Helden)">
Jacques van Helden (jvanheld\@scmbb.ulb.ac.be)
</A>
</FONT>
</CENTER>
</FONT>
EndAddress

print $query->end_html, "\n";

exit(0);








#test
