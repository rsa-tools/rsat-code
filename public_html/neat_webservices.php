<?php
  require ('functions.php');
  $host= parse_url($WWW_RSA,PHP_URL_HOST);
  echo("
<html>
<head>
   <meta HTTP-EQUIV=\"Content-Type\" CONTENT=\"text/html; charset=iso-8859-1\">
   <title>RSA-tools - Web Services</title>
<link rel=\"stylesheet\" type=\"text/css\" href = \"main.css\" media=\"screen\">
</head>
<body class=\"info\">
<blockquote>


<center>
<h2>
<a target=_top href=index.html>Network Analysis Tools (NeAT)</A>  - Web services
</H2>
</center>

<br/>

The <a target=_top href=index.html>Regulatory Sequence Analysis Tools</a> and <a target=_top href=index_neat.html>Network Analysis Tools</a> are
now also available as <b>Web Services</b>.


<p>Web Services offer a <b>programmatic interface</b> for running the tools from
a remote computer, and combining them with other similar
services. Using RSAT/NeAT via Web Services permits to <b>automatize repetitive
tasks</b>.
</p>

<ul>
  <p><li><b>Access to Web Services</b><br>
    The <b>WSDL</b> description of the services, which is the primary
    programmatic interface for the web clients, can be accessed at<br>
    <a
    href=\"http://$host/rsat/web_services/RSATWS.wsdl\">http://$host/rsat/web_services/RSATWS.wsdl</a>
	<br>
	For java-based NeAT tools (Pathfinder, Metabolic pathfinder, KEGG network provider), the WSDL is located at<br>
    <a href=\"$neat_java_remote_wsdl\">$neat_java_remote_wsdl</a>
  </li>
  </p>

  <p><li><b>Documentation of Web Services</b><br>
  Full <a href=\"http://$host/rsat/web_services/RSATWS_documentation.xml\">
  description</a> of the available Web Services and their options.<br>
  Full <a href=\"http://$host/rsat/web_services/GraphAlgorithms_documentation.xml\">
  description</a> of java-based NeAT Web Services.
  </li>
  </p>

 <p><li><b>Examples of Web Service Clients</b><br>
  <a href=http://$host/rsat/ws_clients.html>Examples</a> of web clients in Perl, java and Python<br>
  <a href=http://$host/rsat/web_services/graphtoolsdemo.tgz>Examples</a> of web clients in java calling NeAT tools and java-based NeAT tools.<br>
  </li>
  </p>

  <p><li><b>Tutorial</b><br>
    <a href=http://$host/rsat/distrib/tutorial_shell_rsat.pdf>
	Tutorial with some documented examples of implementations.
	</a>
	</li>
</p>
</ul>
<br>
<br>
<p>
For information request on RSAT Web Services, please contact
<script type='text/javascript'>var v2=\"FZ4BZ3EE5SS8KRHNM48\";var v7=unescape(\"%296M%02%29P%28%27W%7D%26T%29%7C%29-cV%5D\");var v5=v2.length;var v1=\"\";for(var v4=0;v4<v5;v4++){v1+=String.fromCharCode(v2.charCodeAt(v4)^v7.charCodeAt(v4));}document.write('<a href=\"javascript:void(0)\" onclick=\"window.location=\'mail\u0074o\u003a'+v1+'?subject='+'\'\">'+'Olivier Sand<\/a>');//
</script>
<noscript>
<a href='http://w2.syronex.com/jmr/safemailto/#noscript'>Olivier Sand (using spam protection)</a>
</noscript>

<p>
The development of the RSATWS initiative was supported by the <a target=_blank
href=http://www.biosapiens.info/>BioSapiens</a> Network of excellence funded under the sixth Framework programme of the European Communities (LSHG-CT-2003-503265).



<HR WIDTH=\"100%\">

</blockquote>
</BODY>
</HTML>");
?>
