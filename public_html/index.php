<HTML>
<HEAD>
<TITLE>Regulatory Sequence Analysis Tools (RSAT)</TITLE>
<link rel='stylesheet' type='text/css' href='css/bootstrap.min.css' />
<link rel='stylesheet' type='text/css' href='css/home.css' />

<link rel='stylesheet' type='text/css' href='css/font-awesome.css' />
<script src="js/jquery.js"></script>
<script src='https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js'></script>
</head>

<BODY>
<?php
    require_once ('functions.php');
    include('menu.php');
?>

<div class='page-content-wrapper'>
<div class='container'>
	<div class='homeheader'>
		<div class='row'>
			<div class='col-sm-4' align='center'><img src="images/RSAT22.png" style="max-width:250px;max-height:100px" alt="RSAT server" border="0"><br/>
            <span style='color:#0D73A7;font-size:11px;font-family:Arial rounded MT bold'>Regulatory Sequence Analysis Tools</span>
			</div>
			<div class='col-sm-2 homemenu'>
				<a href='publications.cgi'>Publications</a><br/>
				<a href='tutorials.php'>Tutorials</a><br/>
				<a href='people.php'>About us</a><br/>
				<a href='code_of_conduct.html'>Code of Conduct</a>
			</div>
            <div class='col-sm-6 homeheaderbar'><div class='menu-toggle'><a href='#' id='menu-toggle'><i class='fa fa-bars fa-2x' style='color:#F58634; padding-top: 25px;'></i></a></div></div>
		</div>	
	</div>

	<div class='homecontent'>
		<div class='row' style='padding-bottom:20px'>
			<div class='col-sm-9'>
				<div class='row' style='padding-bottom:20px'>
                    <div class='col-sm-2'><br/></div>
                    <div class='col-sm-8'><div align='center' style='font-size:17px;color:#F58634'><b>What we do</b></div><hr style='margin:10px 0 10px 0;border:2px solid #eee'>
					We offer tools to analyse cis-regulatory elements in genome sequences:
						<ul><li>motif discovery (support genome-wide data sets like ChIP-seq)</li>
						<li>transcription factor binding motif analysis (quality assessment, comparisons and clustering)</li>
						<li>comparative genomics</li>
						<li>analysis of regulatory variations</li>
						</ul>
					</div>
                    <div class='col-sm-2'></div>
				</div>
<div align='center' style='font-size:12px;padding-bottom:10px'><i>This website is free and open to all users and there is no login requirement</i></div>
				<div class='row' align='center'>
					<div class='col-sm-6'><div class='panel panel-default' data-toggle='modal' data-target='#programModal'><div class='panel-body'><div class='box-heading'>Which program to use?</div>Guide to main tools for new users<br/><i class='fa fa-cogs fa-3x fa-style'></i></div></div></div>
					<div class='col-sm-6'><div class='panel panel-default' data-toggle='modal' data-target='#tutModal'><div class='panel-body'><div class='box-heading'>Tutorial and help</div>RSAT tutorial and all training material<br/><i class='fa fa-graduation-cap fa-3x fa-style'></i></div></div></div>				
				</div>
				<div class='row' align='center'>
					<div class='col-sm-6'><div class='panel panel-default' data-toggle='modal' data-target='#serverModal'><div class='panel-body'><div class='box-heading'>Choose your server</div><div class='panel-image'><div style='visibility:hidden'>Citing RSAT complete suite of tools<br/><i class='fa fa-graduation-cap fa-3x fa-style'></i></div></div></div></div></div>
					<div class='col-sm-6'><div class='panel panel-default' data-toggle='modal' data-target='#citeModal'><div class='panel-body'><div class='box-heading'>How to cite?</div>Citing RSAT complete suite of tools<br/><i class='fa fa-book fa-3x fa-style'></i></div></div></div>			
				</div>
			</div>
			<div class='col-sm-3'><span class='homepapertext'><span style="color:red"><i class="fa fa-book fa-lg"></span></i> Check <span style="color:red"><b>latest RSAT paper</b></span><a target='_blank' href="https://doi.org/10.1093/nar/gkac312" target="_blank"> in NAR</b></a></span><br/><br/>
			  <span style="color:red"><i class="fa fa-book
			  fa-lg"></span></i>
			  Check <span style="color:red"><b>RSAT
			  Variation-Tools</b></span> <b><a target='_blank'
			  href="https://www.sciencedirect.com/science/article/pii/S2001037019301898?via%3Dihub"
			  target="_blank"> in Computational and
			  Structural Biotechnology
			      Journal</b></a></span><br/><br/>
<a class="twitter-timeline" data-border-color="#F6E6AC" data-chrome="nofooter" data-height="500" href="https://twitter.com/RSATools">Tweets by RSATools</a> <script async src="//platform.twitter.com/widgets.js" charset="utf-8"></script>
			</div>	
		</div>
	</div>
	<div class='homefooter'>
		<div class='social-circle'><div class='row'>
			<div class='col-sm-4'><a class='fa fa-envelope' href='javascript:void(0)' onclick="window.location='mailto:rsat-contact@list01.biologie.ens.fr?subject='"></a> contact RSAT team</div>
			<div class='col-sm-4'><a class='fa fa-twitter' href='https://twitter.com/RSATools' target='_blank'></a> Twitter</div>
			<div class='col-sm-4'><a class='fa fa-code' href='#' style='cursor:default'></a> Code version: <?php echo $properties['git_date'];?></div>
		</div><br/>
		<div class='social-circle'><div class='row'>
			<div class='col-sm-8'><a class='fa fa-paint-brush' href='http://www.altamirastudio.com/' target='_blank'></a> RSAT logos designed by Mauricio Guzman (<a target="_blank" href="http://www.altamirastudio.com.mx/">http://www.altamirastudio.com.mx/</a>)</div>
			<div class='col-sm-4'><a class='fa fa-bug' href='#' style='cursor:default'></a><?php echo " Group specificity: <b>".$properties['group_specificity'] . "</b></p>"; ?></div>
		</div>
		<div align='right'><?php echo $properties['rsat_site'] . ' (' . $properties['rsat_www'] . ')'; ?></div>
	</div></div></div>


<!-- Modal Program to use-->
<div class='modal' id='programModal' role='dialog'>
	<div class='modal-dialog'>
		<div class='modal-content'>
			<div class='modal-header'>
				<button type='button' class='close' data-dismiss='modal'>&times;</button>
				<h4 class='modal-title'>Which program to use?</h4>
			</div>
			<div class='modal-body'>
				<div class="holder">
    <b>1 - Choose your type of data to analyse</b></br>
        <select id="datatypes">
            <option value="">Choose Data type</option>
            <option value="data1">ChIP-seq</option>
            <option value="data2">List of gene names</option>
            <option value="data3">Sequences</option>
            <option value="data4">Matrices (PSSM)</option>
            <option value="data5">Coordinates (BED)</option>
            <option value="data6">List of variants</option>
        </select>
    <br/>
            </div>
        
        <div class="holder">
         <b>2 - Choose your biological question / analysis to perform </b></br>
            <select id="questions" disabled="true">
                <option value="">Choose selection</option>
                <option value="peak-motifs" class="data1">Which TF motifs are over-represented in this data set ?</option>
                
                <option value="retrieve-seq-programs" class="data2">I want to extract their promoter sequences (or other sequence features) </option>
                <option value="footprint-programs" class="data2">Which regulatory elements are conserved in promoters of orthologs ? (only for prokaryotes and fungi)</option>
                
                <option value="disco-programs" class="data3">Are there over-represented motifs in these sequences ? </option>
                <option value="scan-programs" class="data3">I want to scan these sequences with a motif </option>
                
                
                <option value="matrix-scan" class="data4">I want to scan sequences with these matrices </option>
                <option value="matrix-compa-programs" class="data4">I want to compare matrices (with known collections) </option>
                <option value="matrix-compa-programs" class="data4">I want to cluster and align matrices </option>
                <option value="convert-matrix" class="data4">I want to convert the matrix format </option>
                <option value="retrieve-matrix" class="data4">I want to extract specific matrices from a collection</option>
                                                                                                                                                                  
                <option value="fetch-sequences" class="data5">I want to extract the sequences corresponding to these coordinates </option>
                
                 <option value="retrieve-variation-seq" class="data6">Obtain the variants and their flanking sequences </option>
                 <option value="scan-variations" class="data6">Which transcription factor binding sites are affected by these variants ? </option>
            </select>
        </div>
    
     <div class="holder">
         <b>3 - Relevant RSAT programs </b></br>
            <select id="tools" disabled="true">
                <option value="">Choose selection</option>
                <option value="peak-motifs_form.cgi" class="peak-motifs">peak-motifs</option>
                
                <option value="retrieve-seq_form.cgi" class="retrieve-seq-programs">retrieve sequences</option>
                <option value="http://metazoa.rsat.eu/retrieve-ensembl-seq_form.cgi" class="retrieve-seq-programs">retrieve Ensembl sequences (only for organisms supported at Ensembl.org)</option>
                
                 <option value="http://prokaryotes.rsat.eu/footprint-discovery_form.cgi" class="footprint-programs">footprint discovery (Prokaryotes+Fungi)</option>
                 <option value="http://prokaryotes.rsat.eu/footprint-scan_form.cgi" class="footprint-programs">footprint-scan (Prokaryotes+Fungi)</option>
                 
                 <option value="oligo-analysis_form.cgi" class="disco-programs">oligo-analysis (words)</option>
                 <option value="dyad-analysis_form.cgi" class="disco-programs">dyad-analysis (spaced pairs)</option>
                 
                 <option value="dna-pattern_form.cgi" class="scan-programs">dna-pattern</option>
                 <option value="matrix-scan-quick_form.cgi" class="scan-programs">matrix-scan (quick)</option>
                 <option value="matrix-scan_form.cgi" class="scan-programs">matrix-scan (full options)</option>
                 
                 <option value="matrix-scan-quick_form.cgi" class="matrix-scan">matrix-scan (quick)</option>
                 <option value="matrix-scan_form.cgi" class="matrix-scan">matrix-scan (full options)</option>
                 
                 <option value="compare-matrices_form.cgi" class="matrix-compa-programs">compare matrices</option>
                 <option value="matrix-clustering_form.cgi" class="matrix-compa-programs">matrix clustering</option>
                 
                 <option value="convert-matrix_form.cgi" class="convert-matrix">convert matrix</option>
                 
                 <option value="retrieve-seq-bed_form.cgi" class="fetch-sequences">sequences from bed/gff/vcf (locally-installed organisms)</option>
                 <option value="http://metazoa.rsat.eu/fetch-sequences_form.php" class="fetch-sequences">fetch sequences from UCSC (only for organisms supported at genome.ucsc.edu)</option>
                 
                 <option value="http://metazoa.rsat.eu/retrieve-variation-seq_form.cgi" class="retrieve-variation-seq">retrieve variation sequences (Metazoa)</option>
                 <option value="http://metazoa.rsat.eu/variation-scan_form.cgi" class="scan-variations">scan variations (Metazoa)</option>
                 
                 <option value="retrieve-matrix_form.cgi" class="retrieve-matrix">retrieve matrix</option>
            </select>
        </div>
			</div>		
		</div>	
	</div>
</div>

<!-- Modal tutorial -->
<div class='modal' id='tutModal' role='dialog'>
	<div class='modal-dialog'>
		<div class='modal-content'>
			<div class='modal-header'>
				<button type='button' class='close' data-dismiss='modal'>&times;</button>
				<h4 class='modal-title'>Tutorials and Help</h4>
			</div>
			<div class='modal-body'>
				<p><i class="fa fa-graduation-cap fa-lg"></i> <a href="https://rsa-tools.github.io/teaching/" target="tools"><b>All training material</b></a> </p>
	
	<p><i class="fa fa-graduation-cap fa-lg"></i> Learn how to use <b>Peak-motifs</b> with a <b>Nature Protocol</b> <a href='http://www.nature.com/nprot/journal/v7/n8/full/nprot.2012.088.html' target=_blank>[view article]</a></font>
	</p><p><i class="fa fa-graduation-cap fa-lg"></i> Check <b>RSAT tutorial</b> at <b><a target='_blank' href="http://rsa-tools.github.io/tutorial_eccb14/index.html" target="tools">ECCB'14</a>
	</p><p><i class="fa fa-graduation-cap fa-lg"></i><a href='tutorials.php' target='_blank'> All tutorials</a></b>
			</div>
		</div>
	</div>
</div>

<!-- Modal Server -->
<div class='modal' id='serverModal' role='dialog'>
	<div class='modal-dialog modal-lg'>
		<div class='modal-content'>
			<div class='modal-header'>
				<button type='button' class='close' data-dismiss='modal'>&times;</button>
				<h4 class='modal-title'>Choose your server</h4>
			</div>
			<div class='modal-body'>
				<div align="center">
       <table class='serverlist' cellspacing='3' style="font-size:14px;">
        
  	<tr align=center valign=bottom>
  	  
  	  <td>
  	    <a href="http://rsat.france-bioinformatique.fr/fungi/" target="_top">
  	      <img src="images/logo_fungi.jpg" height='100' border='0'
  		   alt="fungi"></a>
  	    <br>hosted by <a target='_blank' href="https://www.france-bioinformatique.fr/">Institut Fran&Ccedil;ais de Bioinformatique (IFB), France </a>
  	  </td>
  
  	  <td>
  	    <a href="http://prokaryotes.rsat.eu/" target="_top">
  	      <img src="images/logo_prokaryotes.jpg" height='100' border='0' alt="prokaryotes"></a>
  	   <br>hosted by <a target='_blank' href="https://www.ccg.unam.mx/genomica-computacional/">Centro de Ciencias Genomicas (CCG), Cuernavaca, Mexico</a>
  	  </td>
  
      <td>
  	    <a href="http://rsat.france-bioinformatique.fr/metazoa/" target="_top">
  	      <img src="images/logo_metazoa.jpg" height='100' border='0' alt="metazoa"></a>
  	   <br>hosted by <a target='_blank' href="https://www.france-bioinformatique.fr/">Institut Fran&Ccedil;ais de Bioinformatique (IFB), France </a>
  	  </td>
  
  	</tr>
  	<tr align=center valign=bottom>
  
  	  <td>
  	    <a href="http://protists.rsat.eu/" target="_top">
  	      <img src="images/logo_protists.jpg" height='100' border='0' alt="protists"></a>
  	    <br>hosted by <a target='_blank' href="https://www.france-bioinformatique.fr/">Institut Fran&Ccedil;ais de Bioinformatique (IFB), France </a>
  	  </td>
  
  	  <td>
  	    <a href="http://plants.rsat.eu/" target="_top">
  	      <img src="images/logo_plants.jpg" height='100' border='0' alt="plants"></a>
  	    <br>deployed by <a target='_blank' href="http://www.eead.csic.es/compbio/staff.html">Bruno Contreras Moreira, Zaraagoza, Spain</a>
  	  </td>
  
      <td>
  	    <a href="http://teaching.rsat.eu" target="_top" >
  	      <img src="images/logo_teaching.jpg" height='100' border='0' alt="teaching"></a>
  	    <br>hosted by <a target='_blank' href="http://france-bioinformatique.fr/">IFB - Institut Fran&ccedil;ais de Bioinformatique (IFB), France</a>
  	  </td>
  	  
  	</tr>

       </table>
</div>
			</div>
		</div>
	</div>
</div>

<!-- Modal citation -->
<div class='modal' id='citeModal' role='dialog'>
	<div class='modal-dialog modal-lg'>
		<div class='modal-content'>
			<div class='modal-header'>
				<button type='button' class='close' data-dismiss='modal'>&times;</button>
				<h4 class='modal-title'>How to cite?</h4>
			</div>
			<div class='modal-body'>
				<span style="color: #cc6600;"><i class="fa fa-pencil fa-lg"></i> <b>Citing RSAT complete suite of tools:</b></span>
	<div align="left">
      <ul>
       <li> Santana-Garcia W, Castro-Mondragon JA, Padilla-G&aacute;lvez M, Nguyen NTT, Elizondo-Salas A, Ksouri N, Gerbes F, Thieffry D, Vincens P, Contreras-Moreira B, van Helden J, Thomas-Chollier M, Medina-Rivera A (2022)
      <b>RSAT 2022: regulatory sequence analysis tools</b>. Nucleic Acids Res. 2022 (Web Server issue)
      <a href='https://doi.org/10.1093/nar/gkac312'>[Full text]</a>
    </li>
      <li> Nguyen NTT, Contreras-Moreira B, Castro-Mondragon J, Santana-Garcia W, Ossio R, Robles-Espinoza CD, Bahin M, Collombet S, Vincens P, Thieffry  D, van Helden J, Medina-Rivera A, Thomas-Chollier M (2018)
	  <b>RSAT 2018: regulatory sequence analysis tools 20th anniversary</b>. Nucleic Acids Res. 2018 (Web Server issue).
	  <a href='https://doi.org/10.1093/nar/gky317'>[Full text]</a>
	</li>
<li>  Medina-Rivera A, Defrance M, Sand O, Herrmann C, Castro-Mondragon J, Delerce J, Jaeger S, Blanchet C, Vincens P, Caron C, Staines DM, Contreras-Moreira B, Artufel M, Charbonnier-Khamvongsa L, Hernandez C, Thieffry D, Thomas-Chollier M, van Helden J (2015)
	  <b>RSAT 2015: Regulatory Sequence Analysis Tools </b>. Nucleic Acids Res. 2015 (Web Server issue) 43:W50–W56.
	  <!--a href=http://www.ncbi.nlm.nih.gov/pubmed/21715389>[Pubmed 21715389]</a-->
	  <a href='http://nar.oxfordjournals.org/content/early/2015/04/21/nar.gkv362.full'>[Full text]</a>
	</li>
	<li>Thomas-Chollier M, Defrance M, Medina-Rivera A, Sand O, Herrmann C, Thieffry D, van Helden J. (2011)
	  <b>RSAT 2011: regulatory sequence analysis tools</b>. Nucleic Acids Res. 2011 Jul;39(Web Server issue):W86-91.
	  <a href=http://www.ncbi.nlm.nih.gov/pubmed/21715389>[Pubmed 21715389]</a>
	  <a href='http://nar.oxfordjournals.org/content/39/suppl_2/W86.long'>[Full text]</a>
	</li>
	<li>Thomas-Chollier, M., Sand, O., Turatsinze, J. V., Janky, R.,
	  Defrance, M., Vervisch, E., Brohee, S. & van Helden, J. (2008). <b>RSAT:
	    regulatory sequence analysis tools</b>. Nucleic Acids
	  Res.
	  <a href=http://www.ncbi.nlm.nih.gov/pubmed/18495751>[Pubmed 18495751]</a>
	  <a href='http://nar.oxfordjournals.org/cgi/content/full/36/suppl_2/W119'>[Full text]</a>
	</li>
	<li>van Helden, J. (2003). <b>Regulatory sequence analysis
	    tools</b>. Nucleic Acids Res. 2003 Jul 1;31(13):3593-6. [<a target=_blank
									href="http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=12824373&dopt=Abstract">Pubmed
	    12824373</a>] [<a target=_blank
			      href=http://nar.oupjournals.org/cgi/content/full/31/13/3593>Full
	    text</a>] [<a target=_blank
			  href=http://nar.oupjournals.org/cgi/reprint/31/13/3593.pdf>pdf</a>]
	</li>
      </ul>
</div>

<i class="fa fa-pencil fa-lg"></i> For citing <b>individual tools</b>: the reference of each tool is indicated on top of their query form.</br>
<i class="fa fa-pencil fa-lg"></i> <a href="publications.cgi">All RSAT publications</a> and all publications for <a href="htmllink.cgi?title=NeAT : Publication&file=neat_publications.html"> the Network Analysis Tools</a>.
			</div>
		</div>
	</div>
</div>

<script type="text/javascript" charset="utf-8"> 

/*! Chained 1.0.0 - MIT license - Copyright 2010-2014 Mika Tuupola */
!function(a,b){"use strict";a.fn.chained=function(c){return this.each(function(){function d(){var d=!0,g=a("option:selected",e).val();a(e).html(f.html());var h="";a(c).each(function(){var c=a("option:selected",this).val();c&&(h.length>0&&(h+=b.Zepto?"\\\\":"\\"),h+=c)});var i;i=a.isArray(c)?a(c[0]).first():a(c).first();var j=a("option:selected",i).val();a("option",e).each(function(){a(this).hasClass(h)&&a(this).val()===g?(a(this).prop("selected",!0),d=!1):a(this).hasClass(h)||a(this).hasClass(j)||""===a(this).val()||a(this).remove()}),1===a("option",e).size()&&""===a(e).val()?a(e).prop("disabled",!0):a(e).prop("disabled",!1),d&&a(e).trigger("change")}var e=this,f=a(e).clone();a(c).each(function(){a(this).bind("change",function(){d()}),a("option:selected",this).length||a("option",this).first().attr("selected","selected"),d()})})},a.fn.chainedTo=a.fn.chained,a.fn.chained.defaults={}}(window.jQuery||window.Zepto,window,document);
	
$(function(){ 
	$("#questions").chained("#datatypes"); 
	$("#tools").chained("#questions"); 
	
	document.getElementById("tools").onchange = function() {
        if (this.selectedIndex!==0) {
            window.location.href = this.value;
        }        
    };
	
});

$("#citeModal a").attr('target', '_blank'); 
</script>
</BODY>
</HTML>
