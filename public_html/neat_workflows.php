<?php
  require ('functions.php');
  $host= parse_url($WWW_RSA,PHP_URL_HOST);
  echo("
<html>
<head>
   <meta HTTP-EQUIV=\"Content-Type\" CONTENT=\"text/html; charset=iso-8859-1\">
   <title>RSA-tools - Web Services - Taverna workflows</title>
<link rel=\"stylesheet\" type=\"text/css\" href = \"main.css\" media=\"screen\">
</head>
<body class=\"info\">
<blockquote>


<center>
<h2><a href=\"NeAT_home.html\">Network Analysis Tools (NeAT)</A> - Taverna workflows
</h2>
NeAT Taverna workflows were developped with the help of Eric Vervisch.
</center>
<br>
<p/>
The <a href=\"http://taverna.sourceforge.net/\" target='_blank'> Taverna </a> workbench provides tools to build, edit and browse workflows. Its graphical interface facilitates an easy building and running of workflows over distributed computer resources.
<p/>
We present below some workflows. They use NeAT and other resources exposed as Web Services. They can be interconnected into reusable  analysis workflows.
<p/>
You can execute those workflows with your own data. Meanwhile data file examples are also available for each workflow.
<p/><b>HOW TO ?</b> To execute these workflows in the Taverna Workbench, download the
workflow (SCUFL XML format) on your machine and load it in Taverna with \"File > Open workflow...\".
Alternatively, copy and paste the URL of the workflow download link into
Taverna with \"File > Open workflow location...\". Sample input data are
provided for each workflow. To use a sample input, download  the
file on your machine, then load this file in Taverna with \"File > Run
workflow >Load input\".

<h3> How to find paths in protein-protein-interaction networks?</h3>

Yeast two-hybrid datasets are known to be very noisy. However, a
protein-protein-interaction confirmed by several different experiments might be more reliable
than one found in only one experiment.
This example demonstrates a workflow that carries out path finding on the intersection
network of two different yeast two-hybrid networks.

<table class=\"normal\" border=0 >
<tr>
<td align=\"\">
<p/>
<a href=\"web_services/taverna/Neat_example_workflow.png\">
<img src=\"web_services/taverna/Neat_example_workflow.png\" width=100>
</a>
</td>
<td valign=\"top\">
In this workflow, the following steps are performed:
<p/>	     - first, the intersection between two undirected networks obtained from yeast two-hybrid data is
               computed with <a href=\"compare_graphs_form.php\">compare-graphs</a>
<p/>	     - second, <a href=\"pathfinder_form.php\">Pathfinder</a> is launched on the intersection network
<p/>

<p/> ---
<p/><ul><li>
<a href=\"web_services/taverna/Neat_example_workflow.xml\"><b>Download Workflow File (SCUFL)</b></a></li>
<p/><li>
<a href=\"web_services/taverna/Neat_example_input.xml\"><b>Download Sample Input</b></a></li>
(The sample input consists of the two yeast two-hybrid data sets from Uetz et al, 2001 and Ito et al, 2002,
which are also available in the <a href=\"demo_files/\">DEMO data</a> section.)
</ul>
</td>
</tr>
</table>

<h3> How to identify biologically relevant clusters in biological networks?</h3>

This workflow forms a part of the Nature protocol presenting NeAT. In the first step,
we extract clusters from a yeast-specific network obtained
from <a href=\"http://string.embl.de/\" target=\"_blank\">STRING</a>.
Next, these clusters are compared to
<a href=\"http://mips.gsf.de/\" target=\"_blank\">MIPS</a> complexes,
which comprise known physical interactions between proteins.
We end up with a list that gives for each MIPS cluster the score of its overlap with
clusters obtained from the STRING network.

<table class=\"normal\" border=0 >
<tr>
<td align=\"\">
<p/>
<a href=\"web_services/taverna/NeAT_Nature_Protocol_workflow.png\">
<img src=\"web_services/taverna/NeAT_Nature_Protocol_workflow.png\" width=100>
</a>
</td>
<td valign=\"top\">
In this workflow, the following steps are performed:
<p/>	     - first, clusters are extracted from the input network using <a href=\"mcl_form.php\">MCL</a>.
<p/>	     - second, the clusters are compared to known protein complexes obtained from MIPS
               with <a href=\"compare_classes_form.php\">compare_classes</a>
<p/>

<p/> ---
<p/><ul><li>
<a href=\"web_services/taverna/NeAT_Nature_Protocol_workflow.xml\"><b>Download Workflow File (SCUFL)</b></a></li>
<p/><li>
<a href=\"web_services/taverna/NeAT_Nature_Protocol_workflow_data.xml\"><b>Download Sample Input</b></a></li>
(The sample input consists of the STRING yeast network and the MIPS classes, both available
in the  <a href=\"http://$host/rsat/data/published_data/nature_protocols/network_analysis/\">Nature protocol data</a> section.)
</ul>
</td>
</tr>
</table>

<h3> How to find paths in organism-specific metabolic networks?</h3>

Path finding in organism-specific metabolic networks allows to detect
pathway variants that are specific to an organism or to a set of related organisms.
<a href=\"http://www.genome.jp/kegg/pathway.html\" target=\"_blank\">KEGG PATHWAY</a>
is a database that offers organism-specific metabolic information.

<table class=\"normal\" border=0 >
<tr>
<td align=\"\">
<p/>
<a href=\"web_services/taverna/Neat_pathfinding_workflow.png\">
<img src=\"web_services/taverna/Neat_pathfinding_workflow.png\" width=100>
</a>
</td>
<td valign=\"top\">
In this workflow, the following steps are performed:
<p/>	     - first, an organism-specific metabolic network is constructed using the <a href=\"keggnetworkprovider_form.php\">KEGG network provider</a>
<p/>	     - second, paths are found in this network using <a href=\"pathfinder_form.php\">Pathfinder</a>
<p/>

<p/> ---
<p/><ul><li>
<a href=\"web_services/taverna/Neat_pathfinding_workflow.xml\"><b>Download Workflow File (SCUFL)</b></a></li>
<p/>
</ul>
</td>
</tr>
</table>

<HR WIDTH=\"100%\">

</blockquote>
</BODY>
</HTML>");
?>
