        <link rel='stylesheet' href='css/simple-sidebar.css'></link>
        <link rel="stylesheet" type="text/css" href="css/menu.css" media="screen,projection,print" />
        <link rel='stylesheet' href='css/colorbox.css'></link>
	<link rel='stylesheet' type='text/css' href='js/autocomplete/css/jquery-ui.css' />
         <script src="js/RSAT_menu.js" type="text/javascript"></script>
        <script src="js/jquery.colorbox-min.js"></script>

        <script>
            $(document).ready(function(){
                $('.iframe').colorbox({iframe:true, innerWidth:'70%', innerHeight:'70%'});
		$("#menu-toggle").click(function(e){
			e.preventDefault();
			$("#wrapper").toggleClass('toggled');
		});
            });

        </script>

<?php
    require_once('functions.php');
?>

<div id="wrapper">
                   <!-- Sidebar -->
            <div id="sidebar-wrapper">
             <div id="menubody" style='padding-top:10px'>

             <!--div id="tabmenu">
             <ul class="rsat">
                 <li><a class="active" href="index.php">RSAT</a></li>
                 <li><a href="NeAT_home.php" >NeAT</a></li>
                 <!--      <li><a href="http://www.rsat.eu/index_neat.html" target="_top">NeAT</a></li>>
             </ul>
             </div> < /#tabmenu -->
             <div id="content" class="rsat">

                 <div class="menu">
                     <h2>
                         <a href="index.php" >
                             <img src="images/RSAT_icon.jpg" style="max-width:150px;max-height:60px;" alt="RSAT server" border="0">
                                 </a>
                     </h2>
                 </div>

                 <!--      <p ><h2 style="border-style:solid;border-color:#cc6600;padding:0 10px;background-color:#F6E6CA;"><a href="http://www.eccb14.org/program/tutorials/cis-r" target="_blank">RSAT tutorial<br>at ECCB'14</a></h2></p>>

                  <p ><h2 style="border-style:solid;border-color:#cc6600;padding:0 10px;background-color:#F6E6CA;">
                  <a target='_blank' href="http://rsat.ulb.ac.be/eccb14/" >RSAT tutorial<br>at ECCB'14</a>
                  </h2></p-->
		<div style='padding:60px 0 0 10px;' align='center'><i class='fa fa-bar-chart  fa-lg'></i> <b>
<?php
    if (php_sapi_name() == 'cli'){
        $out = shell_exec("perl -e 'push @INC, \"$`lib/\"; require \"RSA.lib\";require \"RSA2.cgi.lib\"; @org = &RSAT::OrganismManager::get_supported_organisms_web(); print scalar @org'");
        echo $out;
    }else{
        $path = $properties['rsat_www'].'nbOrg.cgi';
        virtual($path);
        virtual('nbOrg.cgi');
    }
?>
</b> <i>organisms</i></div>
                <div>
			<input type='search' id='searchfunc' placeholder='Search' class='searchmenu' onKeyPress='searchfunc()' onKeyUp='searchfunc()'/>
		</div>
		<div align='right'>
                     <b>New&nbsp;items</b><img src="images/onebit_49.png"  class="new"/>
                 </div>

                 <div class="menu_expand" onclick="expandAll('14')" id="expand"> > view all tools</div>

                 <div class="menu">
                     <div class="menu_heading_closed"
                         onclick="toggleMenu('1')" id="heading1">Genomes and genes</div>
                     <div id="menu1" class="menu_collapsible">
                         <a class="menu_item" href="supported-organisms.cgi">supported organisms</a>
                         <a class="menu_item" href="gene-info_form.cgi"> gene information</a>


<?php

    if ($properties['phylo_tools'] == '1') {
        echo '<a class="menu_item" href="infer-operons_form.cgi"> infer operons</a>';
        echo '<a class="menu_item" href="get-orthologs_form.cgi" >get orthologs</a>';
    }
    else{
        echo '<span class="menu_item" style="cursor:default;color:lightgray"> infer operons</span>';
        echo '<span class="menu_item" style="cursor:default;color:lightgray"> get orthologs</span>';
    }
?>


                         <a class="menu_item" href="random-genes_form.cgi" >random gene selection</a>
                     </div>

                     <div class="menu_heading_closed"
                         onclick="toggleMenu('2')" id="heading2">Sequence tools</div>
                     <div id="menu2" class="menu_collapsible">
                         <a class="menu_item" href="retrieve-seq_form.cgi" >retrieve sequence</a>
<?php

    if($properties['ensembl_tools'] == '1') {
        echo '<a class="menu_item" href="retrieve-ensembl-seq_form.cgi" >retrieve EnsEMBL seq</a><a class="menu_item" href="fetch-sequences_form.php" >fetch-sequences from UCSC</a>'; }
    else{
        echo '<span class="menu_item" style="cursor:default;color:lightgray"> retrieve EnsEMBL seq</span><span class="menu_item" style="cursor:default;color:lightgray"> fetch-sequences from UCSC</span>';
    }
?>

                         <a class="menu_item" href="retrieve-seq-bed_form.cgi" >sequences from bed/gff/vcf</a>
                         <!--	  <a class="menu_item" href="http://www.rsat.eu/retrieve-ensembl-seq_form.cgi" >retrieve EnsEMBL seq</a>-->
                         <a class="menu_item" href="purge-sequence_form.cgi" >purge sequence</a>
                         <a class="menu_item" href="convert-seq_form.cgi" >convert sequences</a>
                         <a class="menu_item" href="random-seq_form.cgi" >random sequences</a>
                     </div>


                     <div class="menu_heading_closed"
                         onclick="toggleMenu('12')" id="heading12">Matrix tools </div>
                     <div id="menu12" class="menu_collapsible">
                         <a class="menu_item" href="retrieve-matrix_form.cgi">retrieve matrix </a>
                         <a class="menu_item" href="convert-matrix_form.cgi" >convert matrix</a>
                         <a class="menu_item" href="compare-matrices_form.cgi" >compare matrices</a>
                         <a class="menu_item" href="matrix-clustering_form.cgi" >matrix-clustering</a>
                         <a class="menu_item" href="matrix-distrib_form.cgi" >matrix distrib</a>
                         <a class="menu_item" href="matrix-quality_form.cgi" >matrix quality</a>
                     </div>


                     <div class="menu_heading_closed"
                         onclick="toggleMenu('11')" id="heading11">Build control sets</div>
                     <div id="menu11" class="menu_collapsible">
                         <a class="menu_item" href="random-genes_form.cgi" >random gene selection</a>
                         <a class="menu_item" href="random-seq_form.cgi" >random sequences</a>
                         <a class="menu_item" href="random-genome-fragments_form.cgi" >random genome fragments</a>
                         <!--	  <a class="menu_item" href="http://www.rsat.eu/random-genome-fragments_form.cgi" >random genome fragments</a>-->
                         <a class="menu_item" href="random-motif_form.cgi" >random-motif</a>
                         <a class="menu_item" href="permute-matrix_form.cgi" >permute-matrix</a>
                         <a class="menu_item" href="random-sites_form.cgi" >random-sites</a>
                         <a class="menu_item" href="implant-sites_form.cgi" >implant-sites</a>
                     </div>
                 </div>

                 <div class="menu">
                     <div class="menu_heading_closed"
                         onclick="toggleMenu('3')" id="heading3">Motif discovery</div>
                     <div id="menu3" class="menu_collapsible">
                         <a class="menu_separator">strings</a>
                         <a class="menu_item" href="oligo-analysis_form.cgi" >oligo-analysis (words)</a>
                         <a class="menu_item" href="oligo-diff_form.cgi" >oligo-diff (words)</a>
                         <a class="menu_item" href="dyad-analysis_form.cgi" >dyad-analysis (spaced pairs)</a>
                         <a class="menu_item_last" href="pattern-assembly_form.cgi" >pattern assembly</a>
                         <a class="menu_separator">strings with positional biais</a>
                         <a class="menu_item" href="position-analysis_form.cgi" >position-analysis (words)</a>
                         <a class="menu_item" href="local-word-analysis_form.cgi" >local-word-analysis (word and spaced pairs)</a>
                         <a class="menu_separator">matrices</a>
                         <a class="menu_item" href="info-gibbs_form.cgi" >info-gibbs</a>
                         <!--<a class="menu_item" href="consensus_form.cgi" >consensus</a>
                         	  <a class="menu_item" href="gibbs_form.cgi" >gibbs</a>-->
                     </div>


                     <div class="menu_heading_closed"
                         onclick="toggleMenu('4')" id="heading4">Pattern matching</div>
                     <div id="menu4" class="menu_collapsible">
                         <a class="menu_separator">matrices</a>
                         <a class="menu_item" href="matrix-scan_form.cgi" >matrix-scan<br>(full options)</a>
                         <a class="menu_item" href="matrix-scan-quick_form.cgi" >matrix-scan (quick)</a>
                         <a class="menu_item" href="crer-scan_form.cgi" >crer-scan</a>
                         <!--<a class="menu_item" href="matrix-enrichment_form.cgi" >matrix-enrichment <img src="images/onebit_49.png" height="30" class="new"></img></a>-->
                         <!--	  <a class="menu_item" href="patser_form.cgi" >patser [discontinued]</a>-->
                         <!--	  <a class="menu_item" href="genome-scale-patser_form.cgi" >genome-scale patser [discontinued]</a>-->
                         <a class="menu_separator">strings</a>
                         <a class="menu_item" href="dna-pattern_form.cgi" >dna-pattern</a>
                         <a class="menu_item_last" href="genome-scale-dna-pattern_form.cgi" >genome-scale dna-pattern</a>
                     </div>
                 </div>

                 <div class="menu">
                     <div class="menu_heading_closed"
                         onclick="toggleMenu('10')" id="heading10">Comparative genomics&nbsp;</img></div>
                     <div id="menu10" class="menu_collapsible">
<?php

    if($properties['phylo_tools'] == '1') {
        echo '<a class="menu_item" href="get-orthologs_form.cgi" >get orthologs</a>'; }
    else{
        echo '<span class="menu_item" style="cursor:default;color:lightgray"> get orthologs</span>';
    }
?>

                         <?php

    if($properties['compara_tools'] == '1') {
        echo '<a class="menu_item" href="get-orthologs-compara_form.cgi" >get orthologs-compara <img src="images/onebit_49.png" height="30" class="new"></img></a>'; }
    else{
        echo '<span class="menu_item" style="cursor:default;color:lightgray"> get orthologs-compara </span>';
    }
?>


<?php

    if($properties['phylo_tools'] == '1') {
        echo '                         <a class="menu_item" href="footprint-discovery_form.cgi" >footprint-discovery</a>
                         <a class="menu_item" href="footprint-scan_form.cgi" >footprint-scan <img src="images/onebit_49.png" height="30" class="new"></img>
                         </a>'; }
    else{
        echo '<span class="menu_item" style="cursor:default;color:lightgray"> footprint-discovery</span><span class="menu_item" style="cursor:default;color:lightgray"> footprint-scan</span>';
    }
?>

                     </div>


                     <div class="menu_heading_closed"
                         onclick="toggleMenu('13')" id="heading13">NGS - ChIP-seq</div>
                     <div id="menu13" class="menu_collapsible">
                         <a class="menu_item" href="peak-motifs_form.cgi" >peak-motifs</a>
                         <!--<a class="menu_item" href="peak-motifs2_form.cgi" >peak-motifs2 - beta</a>-->
                         <?php

    if($properties['ucsc_tools'] == '1') {
        echo '<a class="menu_item" href="fetch-sequences_form.php" >fetch-sequences from UCSC</a>'; }
    else{
        echo '<span class="menu_item" style="cursor:default;color:lightgray"> fetch-sequences from UCSC</span>';
    }
?>
                        <a class="menu_item" href="retrieve-seq-bed_form.cgi" >sequences from bed/gff/vcf </a>

                         <a class="menu_item" href="random-genome-fragments_form.cgi" >random genome fragments</a>
                         <!--	  <a class="menu_item" href="random-genome-fragments_form.cgi" >random genome fragments</a>-->
                     </div>


                     <div class="menu_heading_closed"
                         onclick="toggleMenu('9')" id="heading9">Genetic variations (Var-tools)</div>
                     <div id="menu9" class="menu_collapsible">

                     <?php

    if($properties['variations_tools'] == '1') {
        echo '<!-- <a class="menu_item" href="variation-info_form.cgi" >Variation info </a> -->
                         <a class="menu_item" href="retrieve-variation-seq_form.cgi" >Retrieve variation sequences </img></a>
                         <a class="menu_item" href="variation-scan_form.cgi" >Scan variations with motifs</a>
                         '; }
    else{
        echo '<span class="menu_item" style="cursor:default;color:lightgray"> Variation information </span>
        <span class="menu_item" style="cursor:default;color:lightgray"> Retrieve variation sequences</span>
        <span class="menu_item" style="cursor:default;color:lightgray"> Scan variations with motifs</span>';
    }
?>

                         <a class="menu_item" href="convert-variations_form.cgi" >Convert variation formats </a>
                     </div>
                 </div>


                 <div class="menu">
                     <div class="menu_heading_closed"
                         onclick="toggleMenu('8')" id="heading8">Conversion/Utilities
                         <?php
                             if($properties['ucsc_tools'] == '1'){
                               echo '<img src="images/onebit_49.png"  class="new"/>'; }
                         ?>
                       </div>
                     <div id="menu8" class="menu_collapsible">
                         <a class="menu_separator">Set comparisons / enrichment</a>
                         <!--<a class="menu_item_last" href="compare_classes_form.php?menu=RSAT" >compare classes/clusters</a> -->
                         <a class="menu_separator">Stats</a>
                         <a class="menu_item_last" href="classfreq_form.cgi" >Frequency distribution</a>
                         <a class="menu_separator">sequences</a>
                         <a class="menu_item_last" href="convert-seq_form.cgi" >convert sequences</a>
                         <a class="menu_separator">matrices</a>
                         <a class="menu_item_last" href="convert-matrix_form.cgi" >convert matrix / logo</a>
                         <a class="menu_separator">background models</a>
                         <a class="menu_item_last" href="create-background-model_form.cgi" >create background</a>
                         <a class="menu_item" href="convert-background-model_form.cgi" >convert background</a>
                         <a class="menu_item" href="seq-proba_form.cgi" >sequence probability</a>
                         <a class="menu_separator">features</a>
                         <a class="menu_item" href="convert-features_form.cgi" >convert features</a>
                         <a class="menu_item_last" href="compare-features_form.cgi" >compare features</a>

<?php

    if($properties['ucsc_tools'] == '1') {
      echo '                         <a class="menu_separator">networks</a>
      <a class="menu_item" href="network-interactions_form.cgi" >network interactions <img src="images/onebit_49.png"  class="new"/></a>'; }
    else{
      echo '                         <a class="menu_separator">networks</a>
      <a class="menu_item" style="cursor:default;color:lightgray" >network interactions</a>';
    }
?>
                     </div>


                     <div class="menu_heading_closed"
                         onclick="toggleMenu('5')" id="heading5">Drawing</div>
                     <div id="menu5" class="menu_collapsible">
                         <a class="menu_item" href="feature-map_form.cgi" >feature map</a>
                         <a class="menu_item" href="feature-map_form2.cgi" >feature map 2</a>
                         <a class="menu_item" href="XYgraph_form.cgi" >XY graph</a>
                     </div>


                 </div>

                 <div class="menu">
                     <div class="menu_heading_closed"
                         onclick="toggleMenu('6')" id="heading6">Web services <img src="images/onebit_49.png"  class="new"/></div>
                     <div id="menu6" class="menu_collapsible">
                      <a class="menu_item" href="htmllink.cgi?title=RSAT : REST API &file=web_services_REST.html">REST API <img src="images/onebit_49.png"  class="new"/></a>
                      <a class="menu_separator">deprecated SOAP</a>
                       <a class="menu_item" href="htmllink.cgi?title=RSAT : Web services&file=web_services.html">Programmatic interface to RSAT</a>
                      <a class="menu_item" href="htmllink.cgi?title=RSAT : Web services documentation&file=web_services/RSATWS_documentation.xml">WSDL Documentation</a>


                         <!--<a class="menu_item" href="htmllink.cgi?title=RSAT : Web services&file=web_services.html">Programmatic interface to RSAT</a>
                         <a class="menu_item" href="web_services/RSATWS.wsdl" target=tools>WSDL</a>
                          <a class="menu_item" href="ws_clients.html" target=tools>Clients</a>
                          <a class="menu_item" href="ws_workflows.html" target=tools>Taverna workflows </a-->
                     </div>

                 </div>

                 <div class="menu">


                     <div class="menu_heading_open"
                         onclick="toggleMenu('7')" id="heading7">Help & Contact</div>
                         <div id="menu7" class="menu_collapsible_display">
                         <a class="menu_item" href="people.php" >RSAT team</a>
                         <a class="menu_item" href="https://rsat-doc.github.io/teaching/" target="_blank">Training material</a>
                         <a class="menu_item" href="tutorials.php" >Tutorials</a>

                         <a class="menu_item" href="publications.cgi" >Publications</a>
                         <a class="menu_item" href="credits.php" >Credits</a>
                         <a class="menu_item" href="https://github.com/rsa-tools/" >Download</a>
                         <a class="menu_item_last" href="htmllink.cgi?title=RSAT-motif&file=motif_databases/" >Motif databases</a>
                         <a class="menu_item_last" href="htmllink.cgi?title=RSAT-Data&file=data/" >Data</a>

                        </div>

                         </div>

                        <A target=_top href="#" style='font-size:3px;visibility:hidden'>bug/dont understand</A>

                    </div>
               </div>

            </div><!-- /#sidebar-wrapper -->
