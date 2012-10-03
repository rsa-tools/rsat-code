<!-- Note by JvH: there is not a single dynamic element here, we should check the interest of using a php rather than a html page -->

<html>
<head>
<title>NeAT Home Page</title>
<link rel='stylesheet' type='text/css' href='main.css' media='screen'/>
<meta http-equiv='Content-Type' content='text/html; charset=iso-8859-1' />
</head>
<body class='info'>
<br>
<table class = 'title' cellpadding='10' width ='100%'>
    <tr>


    <td align=center valign = top>
    <font color='#0000dd' size='-2'>
    <a href='http://www.bigre.ulb.ac.be' target=_blank>
    <img src='images/bigre_logo.png' alt='BiGRe lab' border='1' height='75'>
    <br>BiGRe</a>
    </font>
    </td>

    <td align=center>
    <H1>Network analysis tools</H1>
    </td>


    <td align=center valign = top width='160'>
      <font color='#006600' size=-2><a href='http://www.ulb.ac.be' target=_blank>
	  <img src='images/logo_ULB.jpg' alt='ULB logo' border='0' width='75' height='75'>
	  <br>
	  Universit&eacute; Libre de Bruxelles</a>
      </FONT>
    </td>
    </TR>
</TABLE>

     <div class = 'serverlist'>
    <table border=0 cellspacing=3 cellpadding=7 class='serverlist'>
    <TR>
    <td><a href='images/NeAT_flowchart.png'><B>
    Tool Map</B></A></td>

    <td><a href='neat_intro.html'><B>
    Introduction</B></A></td>

    <td><a href='http://www.bigre.ulb.ac.be/forums/'><B>
    Forum</B></A></td>

    <td><a href='neat_tutorial.html'><B>
    Tutorials</B></A></td>


    <td><a href='neat_publications.html'><B>
    Publications</B></A></td>

    <td><a href='neat_credits.html'><B>
    Credits</B></A></td>

    <td><a href='data/'><B>
    Data</B></A></td>

    <td><a href='neat_links.html'><B>

    Links</B></A></td>

    <td><a href='distrib/index.html'><B>
    Download</B></A></td>

    </tr>
    </table>
    </div>
	<br>
	<br>
    Welcome to <b>Network Analysis Tools</b>
    (<B>NeAT</B>). This web site provides a series of modular computer
    programs specifically designed for the analysis of biological networks.


<h3>News</h3>

<h4>New tools</h4>

<ul>
	<li>
	  <p>In the context of the EU-funded
	    MICROME project, focused on the
	    annotation of bacterial metabolism, we
	    developed a simplified interface for the
	    pathway extraction tool, specifically
	    adapted to discover metabolic pathways
	    from sets of functionally related
	    bacterial genes (e.g. co-expression
	    clusters, operons, syntons, synteny
	    groups, ...).</p>
	</li>
</ul>

<h4>Recent publications</h4>		      
<ol>

<p><li><u>Book:</u> Jacques van Helden, Ariane Toussaint and Denis
  Thieffry (2012). Bacterial Molecular Networks. Volume in the series
  Methods in Molecular Biology 804 (28 chapters). 
  [<a target='_blank' href='http://www.springer.com/life+sciences/microbiology/book/978-1-61779-360-8'>Publisher's site</a>]
</li></p>

<p><li>
  van Helden, J., Toussaint, A. and Thieffry, D. (2012). Bacterial
  molecular networks: bridging the gap between functional genomics and
  dynamical modelling. Methods Mol Biol 804, 1-11.
  [<a target='_blank' href='http://www.ncbi.nlm.nih.gov/pubmed/22144145'>PMID 22144145</a>]
</li></p>

<p><li>Lima-Mendez, G. (2012). Reticulate Classification of Mosaic
    Microbial Genomes Using NeAT Website. Methods Mol Biol 804, 81-91.
  [<a target='_blank' href='http://www.ncbi.nlm.nih.gov/pubmed/22144149'>PMID 22144149</a>]
</li></p>

<p><li>
  Faust, K. and van Helden, J. (2012). Predicting Metabolic Pathways
  by Sub-network Extraction. Methods Mol Biol 804, 107-30.
  [<a target='_blank' href='http://www.ncbi.nlm.nih.gov/pubmed/22144151'>PMID 22144151</a>]
</li></p>

<p><li>Broh&eacute;e, S. (2012). Using the NeAT Toolbox to Compare Networks
    to Networks, Clusters to Clusters, and Network to Clusters. Methods
    Mol Biol 804, 327-42.
    [<a target='_blank' href='http://www.ncbi.nlm.nih.gov/pubmed/22144162'>PMID 22144162</a>]
</li></p>

<li>Faust, K., Croes, D. and van Helden, J. (2011). Prediction of
  metabolic pathways from genome-scale metabolic networks. Biosystems
  105, 109-21.
  [<a target='_blank' href='http://www.ncbi.nlm.nih.gov/pubmed/21645586'>PMID 21645586</a>]
  [<a target='_blank' href='http://www.sciencedirect.com/science/article/pii/S0303264711000839'>doi:10.1016/j.biosystems.2011.05.004</a>]
</li>

<li>Faust, K., Dupont, P., Callut, J. and van Helden,
  J. (2010). Pathway discovery in metabolic networks by subgraph
  extraction. Bioinformatics 26:1211-8. <a target='_blank'
					   href='http://www.ncbi.nlm.nih.gov/pubmed/20228128'>[Pubmed 20228128]</a>
</li>

<li><p><a href='publications.html'> ... other
publications </a></p></li>

</ol>

<p>This website is free and open to all users.</p>


<table align = 'center'><tr><td>
<img src='images/neat_logo2_shadow.png' border=0 height=40%>
</td></tr></table>

<HR  align = 'left'>

<H4 class='footer'><i>
For suggestions or information request, please contact :
<BR>
<a href='mailto:sylvain-at-bigre.ulb.ac.be'>
Sylvain Broh&eacute;e (sylvain-at-bigre.ulb.ac.be)</i>
</A>
</H4>

</body>

</html>
<?php
  require ('functions.php');
  echo("");
?>
