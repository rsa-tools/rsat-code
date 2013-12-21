<?php
  require ('functions.php');
  echo("
<HTML>
	<head>
		<meta http-equiv='Content-Type' content='text/html; charset=iso-8859-1'>
		<title>Network Analysis Tools - menu</title>
		<link rel='stylesheet' type='text/css' href='menu.css' media='screen' />
		<script src='RSAT_menu.js' type='text/javascript'></script>
	</head>
	<body>
<ul id='tabmenu' class='rsat'>
  <li><a href='index.html' target='_top'>RSAT</a></li>
  <li><a class='active' href='index_neat.html' target='_top'>NeAT</a></li>
</ul>

<div id='content' class='graph'>

<h2><img src='images/neat_logo_small.png' /><a href='NeAT_home.php' target='tools'>Network Analysis Tools</a></h2>


<div class='menu'>
<div class='menu_heading_open'
onclick='toggleMenu('201')' id='heading201'>Networks</div>
<div id='menu201' >
	<a class='menu_item' href='compare_graphs_form.php' target='tools'>Network comparison</a>
	<a class='menu_item' href='graph_topology_form.php' target='tools'>Node topology statistics</a>
	<a class='menu_item' href='graph_neighbours_form.php' target='tools'>Get node neighborhood</a>
	<a class='menu_item' href='random_graph_form.php' target='tools'>Randomize network</a>
	<a class='menu_item' href='alter_graph_form.php' target='tools'>Alter network</a>
	<a class='menu_item' href='convert_graph_form.php' target='tools'>Format conversion / layout calculation</a>
	<a class='menu_item_last' href='display_graph_form.php' target='tools'>Graph display</a>
</div>
</div>
<div class='menu'>
<div class='menu_heading_open'
onclick='toggleMenu('202')' id='heading202'>Path finding and pathway extraction</div>
<div id='menu202' >
        <a class='menu_item' target='tools' href='$neat_java_host/rsat/pathfinder_form.php'>k-shortest path finding </a>
	<a class='menu_item' target='tools' href='$neat_java_host/metabolicpathfinding/metabolicPathfinder_form.jsp'>Metabolic path finding</a>
	<a class='menu_item_last' target='tools' href='$neat_java_host/metabolicpathfinding/pathwayinference_form.jsp'>Pathway extraction</a>
	<a class='menu_item_last' target='tools' href='pathway-extractor_form.php'>Pathway extraction (microme prototype)</a>
</div>

</div>


<div class='menu'>
<div class='menu_heading_open'
onclick='toggleMenu('203')' id='heading203'>Clusters</div>
<div id='menu203' >
	<a class='menu_item' href='compare_classes_form.php' target='tools'>Compare classes/clusters</a>
	<a class='menu_item_last' href='contingency_stats_form.php' target='tools'>Contingency stats</a>
	<a class='menu_item_last' href='convert_classes_form.php' target='tools'>Convert classification to different format</a>

</div>
</div>

<div class='menu'>
<div class='menu_heading_open'
onclick='toggleMenu('204')' id='heading204'>Clusters/Networks</div>
<div id='menu204' >
<!--	<a class='menu_item' href='graph_clique_form.php' target='tools'>Find cliques</a>-->
	<a class='menu_item' href='mcl_form.php' target='tools'>MCL clustering</a>
	<a class='menu_item' href='rnsc_form.php' target='tools'>RNSC clustering</a>
	<a class='menu_item' href='graph_get_clusters_form.php' target='tools'>Map clusters onto network</a>
	<a class='menu_item_last' href='graph_cluster_membership_form.php' target='tools'>Cluster membership</a>

</div>
</div>

<div class='menu'>
<div class='menu_heading_open'
onclick='toggleMenu('206')' id='heading206'>Graphics</div>
<div id='menu206' >
	<a class='menu_item' href='draw_heatmap_form.php' target='tools'>Draw heatmap</a>
</div>
</div>

<div class='menu'>
<div class='menu_heading_open'
onclick='toggleMenu('205')' id='heading205'>Roc curves</div>
<div id='menu205' >
	<a class='menu_item' href='roc-stats_form.cgi' target='tools'>ROC curves and stats</a>

</div>
</div>

<div class='menu'>
<div class='menu_heading_open'
onclick='toggleMenu('206')' id='heading206'>Data</div>
<div id='menu206' >
	<a class='menu_separator'>Collect data from external servers</a>
	<a class='menu_item' href='string_dataset_form.php' target='tools'>Download a subgraph from STRING</a>
	<a class='menu_item' href='$neat_java_host/rsat/keggnetworkprovider_form.php' target='tools'>Download organism-specific networks from KEGG</a>
	<a class='menu_separator'>Sample data</a>
	<a class='menu_item' href='http://rsat.ulb.ac.be/rsat/demo_files/' target='tools'>DEMO data</a>
	<a class='menu_item_last' href='http://rsat.ulb.ac.be/rsat/data/published_data/nature_protocols/network_analysis/'
	target='tools'>Nature protocol data</a>
</div>
</div>


<div class='menu'>
<div class='menu_heading_open'
onclick='toggleMenu('207')' id='heading207'>Web services</div>
	<div id='menu207' >
		<a class='menu_item' href='$neat_java_host/rsat/neat_webservices.php' target='tools'>Documentation and clients</a>
		<a class='menu_item_last' href='neat_workflows.php' target='tools'>Taverna workflows</a>
	</div>
</div>

<div class='menu'>
<div class='menu_heading_open'
onclick='toggleMenu('208')' id='heading208'>Help</div>
	<div id='menu208' >
	<a class='menu_item' href='images/NeAT_flowchart.png' target='tools'>Site map</a>
	<a class='menu_item' href='neat_tutorial.html' target='tools'>Tutorial</a>
	<a class='menu_item' href='http://www.bigre.ulb.ac.be/forums/' target='tools'>Contact & Forum</a>
	</div>
</div>

<div class='menu'>
<div class='menu_heading_open'
onclick='toggleMenu('209')' id='heading209'>Information</div>
	<div id='menu209' >
	<a class='menu_item' href='neat_intro.html' target='tools'>Introduction</a>
	<a class='menu_item' href='people.html' target='tools'>People</a>
	<a class='menu_item' href='neat_publications.html' target='tools'>Publications</a>
	<a class='menu_item' href='distrib/index.html' target='tools'>Distribution</a>
	<a class='menu_item' href='neat_credits.html' target='tools'>Credits</a>
	<a class='menu_item' href='neat_links.html' target='tools'>Links</a>
	</div>
</div>

<h3>Feedback </h3>
<h3>
<a href='mailto:sylvain-at-bigre.ulb.ac.be'>Sylvain Broh&eacute;e</a>
</h3>
</h4>
</div>
</body>
</html>");
?>
