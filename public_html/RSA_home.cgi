#!/usr/bin/perl
#### this cgi script fills the HTML form for the program dna-pattern
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib.pl";
require "RSA_cgi_lib.pl";


$query = new CGI;
print $query->header;
print $query->start_html(-title=>'Regulatory Sequence Analysis Tools',
			 -author=>'jvanheld@ucmb.ulb.ac.be',
			 -BGCOLOR=>'#FFEEDD');

print <<EndText;

  <FONT FACE="Helvetica">
  
  <TABLE CELLPADDING=10 ALIGN=CENTER>
  <TR>
  <TD VALIGN=TOP WIDTH="160">
  <CENTER>
  <A HREF="http://www.cifn.unam.mx/Computational_Biology/index.html" target=_top>
  <IMG SRC="lablogo.gif" ALT="lab logo" BORDER=0 HEIGHT=48 WIDTH=75></A>
  <BR>
  <FONT COLOR="#006600" SIZE=-2>
  LABORATORIO DE BIOLOGIA COMPUTACIONAL
  </FONT>
  </CENTER>
  </TD>
  
  <TD ALIGN=CENTER VALIGN=CENTER>
  <FONT FACE="Helvetica,Arial" SIZE=+1><B>
  Regulatory Sequence Analysis Tools</B></FONT>
  </TD>
  
  
  
  <TD VALIGN=TOP WIDTH="160">
  <CENTER>
  <A HREF="http://www.ucmb.ulb.ac.be" TARGET=_top>
  <IMG SRC="ULB_Blue_215.gif" ALT="ULB_Blue_215.gif" BORDER=0 HEIGHT=75 WIDTH=75>
  <BR>
  <FONT COLOR="#0000dd" SIZE=-2>
  UCMB-ULB
  </FONT>
  </A>
  </CENTER>
  </TD>
  </TR>
  </TABLE>
  
  <HR SIZE=4 WIDTH="100%">
  
  <TABLE BORDER=0 CELLSPACING=3 CELLPADDING=7 ALIGN=CENTER BGCOLOR="#FFDDCC">
  <TR>
  <TD><FONT FACE="Helvetica" SIZE=-1><A HREF="intro.html"><B>
  Introduction</B></A></FONT></TD>
  
  <TD><FONT FACE="Helvetica" SIZE=-1><A HREF="tutorial.html"><B>
  Tutorial</B></A></FONT></TD>
  
  <TD><FONT FACE="Helvetica" SIZE=-1><A HREF="data/"><B>
  Data</B></A></FONT></TD>
  
  <TD><FONT FACE="Helvetica" SIZE=-1><A HREF="publications.html"><B>
  Publications</B></A></FONT></TD>
  
  <TD><FONT FACE="Helvetica" SIZE=-1><A HREF="credits.html"><B>
  Credits</B></A></FONT></TD>
  
  <TD><FONT FACE="Helvetica" SIZE=-1><A HREF="links.html"><B>
  Links</B></A></FONT></TD>
  
  </TR>
  </TABLE>
  
  
  <FONT SIZE=-1>
  
  <CENTER>
  <P>
  <BLINK>
  <FONT COLOR=#FF6666 SIZE=+2>
  <B>New ! 16 bacteria supported</B>
  </FONT>
  </BLINK>
  </CENTER>
  
  <P>
  Welcome to Regulatory Sequence Analysis Tools. This site provides a
  series of modular computer programs specifically designed for the
  detection of regulatory signals in non-coding sequences.
  
  

  <P> Initially devoted to the yeast <I>Saccharomyces cerevisiae</I>,
  Regulatory Sequence Analysis tools have now been extended to 16
  completely sequences bacterial genomes.  Most tools present now a
  menu allowing to select the organism on which the analysis will be
  performed. The extension is ready for most tools, but some still
  have to be updated (ORF-info, neighbour-orfs, downstream sequence).


  <P>Basically, the same tools can be used for bacteria as for
  yeast. For pattern discovery in bacteria, we strongly recommend the
  <A HREF="dyad-detector_form.cgi">dyad-detector</A>, which has been
  specifically designed for the detection of HTH-like sites. The
  dyad-detector is also recommended for yeast, where it allows to
  detect the sites bound by Zn cluster-containing factors.

  <P> This web site is freely available for academic users. For users
  from commercial companies, please read our <A
  HREF="disclaimer.html">disclaimer</A>.

  <P>
  <B>RSA-tools on the web</B>
  
  <TABLE BORDER=0 CELLPADDING=10 ALIGN=CENTER>
  <TR ALIGN=CENTER VALIGN=BOTTOM>

  <TD>
  <A HREF="http://www.ucmb.ulb.ac.be/bioinformatics/rsa-tools/" target="_top">
  <B>Belgium</B><BR>
  <IMG SRC="manneken_pis.jpg" HEIGHT=80 WIDTH=53 BORDER=0><BR>
  <FONT SIZE=-2>
 http://www.ucmb.ulb.ac.be/bioinformatics/rsa-tools/</A>
  </FONT>
  </TD>

  <TD ALIGN=CENTER>
  <A HREF="http://copan.cifn.unam.mx/Computational_Biology/yeast-tools" target="_top">
  <B>Mexico (yeast only)</B><BR>
  <IMG SRC="zapata.jpg" HEIGHT=80 WIDTH=89 BORDER=0><BR>
  <FONT SIZE=-2>http://copan.cifn.unam.mx/Computational_Biology/yeast-tools</A>
  </FONT>
  </TD>
  
  </TR>
  </TABLE>
  
EndText


print "<H4>Organisms supported on <A HREF='$WWW_RSA'>",$WWW_RSA,"</A></H4>\n";
print &ListSupportedOrganisms("html_list");

&UpdateLogFile;
print "This page has been visited\n<B>";
print &UpdateCounterFile;
print "</B>\ntimes since May 1998.\n";
print "Enjoy your visit !";


print <<EndAddress;
<P>
<HR WIDTH="100%">
<CENTER>
<FONT SIZE=-1>
For suggestions or information request, please contact :
<BR>
<A HREF="mailto:jvanheld\@ucmb.ulb.ac.be (Jacques van Helden)">
Jacques van Helden (jvanheld\@ucmb.ulb.ac.be)
</A>
</FONT>
</CENTER>
</FONT>
EndAddress

print $query->end_html, "\n";

exit(0);








