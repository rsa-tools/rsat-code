<?php echo("

<link rel='stylesheet' href='css/simple-sidebar.css'></link>
        <link rel='stylesheet' type='text/css' href='menu.css' media='screen,projection,print' />
        <link rel='stylesheet' href='css/colorbox.css'></link>
        <script src='js/jquery.js'></script>
        <script src='RSAT_menu.js' type='text/javascript'></script>
		<script src='js/jquery.colorbox-min.js'></script>

        <script>
            $(document).ready(function(){
                $('.iframe').colorbox({iframe:true, innerWidth:'70%', innerHeight:'70%'});
            });
        </script>

<div id='wrapper'> 
                   <!-- Sidebar -->
            <div id='sidebar-wrapper'>
             <div id='menubody'>
             <div id='tabmenu'>

<ul class='rsat'>
  <li><a href='RSAT_home.cgi' target='_top'>RSAT</a></li>
  <li><a class='active' href='NeAT_home.php' target='_top'>NeAT</a></li>
</ul>
</div> <!-- /#tabmenu -->
<div id='content' class='graph'>

<h2><img src='images/neat_logo_small.png' /><a href='NeAT_home.php' >Network Analysis Tools</a></h2>

<div class='menu_expand' onclick=\"expandAll('10')\" id='expand'> > view all tools</div>
<div class='menu'>
<div class='menu_heading_closed'
onclick=\"toggleMenu('1')\" id='heading1'>Networks</div>
<div id='menu1' class='menu_collapsible'>
	<a class='menu_item' href='compare_graphs_form.php' >Network comparison</a>
	<a class='menu_item' href='graph_topology_form.php' >Node topology statistics</a>
	<a class='menu_item' href='graph_neighbours_form.php' >Get node neighborhood</a>
	<a class='menu_item' href='random_graph_form.php' >Randomize network</a>
	<a class='menu_item' href='alter_graph_form.php' >Alter network</a>
	<a class='menu_item' href='convert_graph_form.php' >Format conversion / layout calculation</a>
	<a class='menu_item_last' href='display_graph_form.php' >Graph display</a>
</div>
</div>
<div class='menu'>
<div class='menu_heading_closed'
onclick=\"toggleMenu('2')\" id='heading2'>Path finding and pathway extraction</div>
<div id='menu2' class='menu_collapsible'>
        <a class='menu_item'  href='$neat_java_host/pathfinder_form.php'>k-shortest path finding </a>
	<a class='menu_item'  href='http://wwwsup.scmbb.ulb.ac.be/metabolicpathfinding/metabolicPathfinder_form.jsp'>Metabolic path finding</a>
	<a class='menu_item_last'  href='http://wwwsup.scmbb.ulb.ac.be/metabolicpathfinding/pathwayinference_form.jsp'>Pathway extraction</a>
	<a class='menu_item_last'  href='pathway-extractor_form.php'>Pathway extraction (microme prototype)</a>
<!--
        <a class='menu_item'  href='$neat_java_host/pathfinder_form.php'>k-shortest path finding </a>
	<a class='menu_item'  href='$neat_java_host/metabolicpathfinding/metabolicPathfinder_form.jsp'>Metabolic path finding</a>
	<a class='menu_item_last'  href='$neat_java_host/metabolicpathfinding/pathwayinference_form.jsp'>Pathway extraction</a>
	<a class='menu_item_last'  href='pathway-extractor_form.php'>Pathway extraction (microme prototype)</a>
-->
</div>

</div>


<div class='menu'>
<div class='menu_heading_closed'
onclick=\"toggleMenu('3')\" id='heading3'>Clusters</div>
<div id='menu3' class='menu_collapsible'>
	<a class='menu_item' href='compare_classes_form.php' >Compare classes/clusters</a>
	<a class='menu_item_last' href='contingency_stats_form.php' >Contingency stats</a>
	<a class='menu_item_last' href='convert_classes_form.php' >Convert classification to different format</a>

</div>
</div>

<div class='menu'>
<div class='menu_heading_closed'
onclick=\"toggleMenu('4')\" id='heading4'>Clusters/Networks</div>
<div id='menu4' class='menu_collapsible'>
<!--	<a class='menu_item' href='graph_clique_form.php' >Find cliques</a>-->
	<a class='menu_item' href='mcl_form.php' >MCL clustering</a>
	<a class='menu_item' href='rnsc_form.php' >RNSC clustering</a>
	<a class='menu_item' href='graph_get_clusters_form.php' >Map clusters onto network</a>
	<a class='menu_item_last' href='graph_cluster_membership_form.php' >Cluster membership</a>

</div>
</div>

<div class='menu'>
<div class='menu_heading_closed'
onclick=\"toggleMenu('5')\" id='heading5'>Graphics</div>
<div id='menu5' class='menu_collapsible'>
	<a class='menu_item' href='draw_heatmap_form.php' >Draw heatmap</a>
</div>
</div>

<div class='menu'>
<div class='menu_heading_closed'
onclick=\"toggleMenu('6')\" id='heading6'>Roc curves</div>
<div id='menu6' class='menu_collapsible'>
	<a class='menu_item' href='roc-stats_form.cgi' >ROC curves and stats</a>

</div>
</div>

<div class='menu'>
<div class='menu_heading_closed'
onclick=\"toggleMenu('7')\" id='heading7'>Data</div>
<div id='menu7' class='menu_collapsible'>
	<a class='menu_separator'>Collect data from external servers</a>
	<a class='menu_item' href='string_dataset_form.php' >Download a subgraph from STRING</a>
	<a class='menu_item' href='$neat_java_host/rsat/keggnetworkprovider_form.php' >Download organism-specific networks from KEGG</a>
	<a class='menu_separator'>Sample data</a>
	<a class='menu_item' href='http://www.rsat.eu/demo_files/' >DEMO data</a>
	<a class='menu_item_last' href='http://www.rsat.eu/data/published_data/nature_protocols/network_analysis/'
	>Nature protocol data</a>
</div>
</div>


<div class='menu'>
<div class='menu_heading_closed'
onclick=\"toggleMenu('8')\" id='heading8'>Web services</div>
	<div id='menu8' class='menu_collapsible'>
		<a class='menu_item' href='$neat_java_host/rsat/neat_webservices.php' >Documentation and clients</a>
		<a class='menu_item_last' href='neat_workflows.php' >Taverna workflows</a>
	</div>
</div>

<div class='menu'>
<div class='menu_heading_open'
onclick=\"toggleMenu('9')\" id='heading9'>Help</div>
	<div id='menu9' >
	<a class='menu_item' href='htmllink.cgi?title=NeAT : Site map&file=images/NeAT_flowchart.png' >Site map</a>
	<a class='menu_item' href='htmllink.cgi?title=NeAT : Tutorial&file=neat_tutorial.html' >Tutorial</a>
	<a class='menu_item' href='http://www.bigre.ulb.ac.be/forums/' >Contact & Forum</a>
	</div>
</div>

<div class='menu'>
<div class='menu_heading_open'
onclick=\"toggleMenu('10')\" id='heading10'>Information</div>
	<div id='menu10' >
	<a class='menu_item' href='htmllink.cgi?title=NeAT : Introduction&file=neat_intro.html' >Introduction</a>
	<a class='menu_item' href='htmllink.cgi?title=NeAT : People&file=people.html' >People</a>
	<a class='menu_item' href='htmllink.cgi?title=NeAT : Publication&file=neat_publications.html' >Publications</a>
	<a class='menu_item' href='htmllink.cgi?title=NeAT : Distribution&file=distrib/index.html' >Distribution</a>
	<a class='menu_item' href='htmllink.cgi?title=NeAT : Credits&file=neat_credits.html' >Credits</a>
	<a class='menu_item' href='htmllink.cgi?title=NeAT : Links&file=neat_links.html' >Links</a>
	</div>
</div>

<h3>Feedback <br/>
<a href='mailto:sylvain-at-bigre.ulb.ac.be'>Sylvain Broh&eacute;e</a>
</h3>
</h4>
</div>
</div>
 </div><!-- /#sidebar-wrapper -->
");?>