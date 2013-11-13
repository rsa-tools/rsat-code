<?php
## Load RSAT configuration
require ('functions.php');
#print_r($properties);
UpdateLogFile("rsat","","");
?>
<script type="text/javascript">
	var multiz_genome = new Array();
	var multiz_path = new Array();
	multiz_genome[0] = new Array();
	multiz_path[0] = 'none';
<?php
//////////////////////////
//Get info on avalaible genome on ucsc
$file = fopen ($properties['RSAT'].'/public_html/data/supported_organisms_ucsc.tab', 'r');
$genome_ucsc= array();

while (!feof($file)) {
	list($genome, $organism, $species, $date) = explode("\t", substr(fgets($file),0,-1));
	if ($genome != "") {
		$genome_ucsc[$genome]['organism'] = $organism;
		$genome_ucsc[$genome]['species'] = $species;
		$genome_ucsc[$genome]['date'] = $date;
	}
}
fclose(file);

//Get multiz alignements supported on this local server (ref genomes + aligned species)
$multiz_supported = array();
$filename = $properties['RSAT'].'/public_html/data/supported_organisms_ucsc_multiz_local.tab';

if (is_file($filename)) {
  $file = fopen ($filename, 'r');
  while (!feof($file)) {
	list($genome_reference,$nb_species,$genome,$path) = explode("\t", substr(fgets($file),0,-1));
	if ($genome_reference != "") {
		$multiz_supported[$genome_reference]['aligned_species'][] = $genome;
		$multiz_supported[$genome_reference]['path']= $path;
	}
  }
  fclose(file);
}

//Write array() with key multiz available genome and values genome in multiz.
foreach($multiz_supported as $genome_reference => $value) {
	echo "multiz_genome['$genome_reference'] = new Array(";

	foreach($multiz_supported[$genome_reference]['aligned_species'] as $key => $species) {
		if ($species == end($multiz_supported[$genome_reference]['aligned_species'])) {
			echo "'$species - ".$genome_ucsc[$species]['organism']."'";
		}	else {
			echo "'$species - ".$genome_ucsc[$species]['organism']."',";
		}
	}
	echo ")\n";

	if ( is_dir( $properties['RSAT']."/data/UCSC_multiz/".$multiz_supported[$genome_reference]['path'])) {
		echo "multiz_path['$genome_reference'] = '".$multiz_supported[$genome_reference]['path']."';\n";
	} else {
		echo "multiz_path['$genome_reference'] = '';\n";
	}	
}

 ?>	
 function fill_species_ali(code,fill_area,clear_area) {
	 var genomes = multiz_genome[code];
	 if (code == 0) {
		 on_reset(fill_area);
		 on_reset(clear_area);
	 } else {
		 document.getElementById(fill_area).options.length = genomes.length;
		 on_reset(clear_area);

		 for (i=0; i<genomes.length; i++) {
			 document.getElementById(fill_area).options[i].value = genomes[i].split(" ")[0];
			 document.getElementById(fill_area).options[i].id = genomes[i].split(" ")[0];
			 document.getElementById(fill_area).options[i].text = genomes[i];
		 }		
	 }	
 }

 function fill_multiz_path(ref_species) {
	 var path = multiz_path[ref_species];
	 if (path == '') {
		 document.getElementById('maf_path').type = "text";
		 document.getElementById('maf_path').value = "Can't find maf for "+ref_species+" on the server. Please contact system administrator.";
	 } else {
		 document.getElementById('maf_path').type = "hidden";
		 document.getElementById('maf_path').value = path;

	 }
 }

 function deselect_all(id) {
	 if (document.getElementById(id).options.length != 0) {
		 document.getElementById(id).options.selectedIndex = 0;
		 document.getElementById(id).options[0].selected = false;
	 }
 }

 function on_reset(id) {
	 document.getElementById(id).options.length = 0;
 }

 function select_all(id) {
	 for (i=0; i<document.getElementById(id).options.length; i++) {
		 document.getElementById(id).options[i].selected = true;
	 }
 }

 function move(from, to) {
	 for (i=0; i<document.getElementById(from).options.length; i++) {
		 if (document.getElementById(from).options[i].selected) {
			 document.getElementById(to).options.length = document.getElementById(to).options.length +1;
			 document.getElementById(to).options[document.getElementById(to).options.length-1].value = document.getElementById(from).options[i].value;
			 document.getElementById(to).options[document.getElementById(to).options.length-1].text = document.getElementById(from).options[i].text;
			 document.getElementById(from).options[i]=null;
			 i--;
		 }
	 }					
 }	

 function add_demo(ref_species, align_species, bed_file) {

	 document.forms[0].bed_url.value = bed_file;
	 document.forms[0].bed.value="";
	 document.getElementById('bed_file').value=''; 

	 document.getElementById(ref_species).selected = true;
	 fill_species_ali(ref_species,'species_ali','species_ali_keep');
	 fill_multiz_path(ref_species)
	 for (Species in align_species) {
		 document.getElementById(align_species[Species]).selected = true;
	 }
	 move('species_ali','species_ali_keep');

	 document.forms[0].r_motif.value="MA0102.2";
	 document.forms[0].custom_motifs.value="";
	 document.getElementById('custom_motifs_file').value=''; 

	 document.forms[0].nb_peaks.value=300;
	 document.forms[0].max_chi2_pvalue.value=1;
	 document.forms[0].max_hyp_pvalue.value=0.001;	
	 document.forms[0].cons_thres.value=0.7;
	 document.forms[0].window_cons_thres.value=0.7;
	 document.forms[0].window_size.value=5;
	 document.forms[0].motif_number.value=50;
	 document.forms[0].motif_number_family.value=4;

	 document.forms[0].output[0].checked=true;
 }

 </script>
 <?php
 //Print header
 $prog_name = "peak-footprints";
 $result = false;
 require ('RSAT_header.php');

 if ($multiz_supported == array()) {
     error("<span style='color:red'>No supported multiz on this server. Please contact system administrator.</span>");
 }

 ?>
 <form id="form" method='post' action='peak-footprints.php' enctype='multipart/form-data' onreset="on_reset('species_ali');on_reset('species_ali_keep');" onsubmit="select_all('species_ali_keep')">
   <fieldset>  
   <legend><b>Genomic sequences</b></legend>
   <p>
   <b>Genomic coordinates</b> 
   <span style='color:orange'>(mandatory)</span>
   <br/>Should be provided as a bed file (<a target='_blank' href='http://genome.ucsc.edu/FAQ/FAQformat.html#format1'>bed format</a>), in any of the three following ways:				
   <ul style="list-style-type:square">
   <li>Paste coordinates<br/><textarea name='bed' rows='6' cols='45'></textarea></li>
   <li>Specify the URL of bed file on remorte server (e.g. Galaxy)<br/><input type="text" name="bed_url" size="52" /></li>
   <li>Upload a file from your computer<br/><input type='file' name='bed_file' id='bed_file' size='40' /></li>
   </ul>
   </p>
   <p>
   <b>Reference genome</b> 
   <span style='color:orange'>(mandatory)</span>
   &nbsp;&nbsp;
<select name="genome" id="genome" onChange="fill_species_ali(this.options[this.selectedIndex].value,'species_ali','species_ali_keep');fill_multiz_path(this.options[this.selectedIndex].value);">
	   <option value="0">----Choose a genome----</option>								
	   <?php
	   foreach($multiz_supported as $key => $value) {
   echo "<option value='$key' id='$key'>$key - ".$genome_ucsc[$key]['organism']."</option>", "\n";
 }
?>
</select>
<input type='hidden' name='maf_path' id='maf_path' style='background-color:#EFFBFB;border-style:none;color:red' size="100%" readonly="readonly"/>
	   </p>									

	   <p>
	   <b>Aligned species</b> <span style='color:orange'>(mandatory)</span><br/>
	   <table>
	   <tr>
	   <td><select name="species_ali[]" id="species_ali" style="width:200px" size="8" multiple="multiple"/></select></td>
	   <td><input type="button" value="Add >>" name="add" onclick="move('species_ali','species_ali_keep');" style="width:80px"/><br/>
	   <input type="button" style="width:80px" value="&lt;&lt; Remove" name="remove" onclick="move('species_ali_keep','species_ali');" size='10'/></td>
	   <td><select name="species_ali_keep[]" id="species_ali_keep" style="width:200px" size="8" multiple="multiple"/></select></td>
	   </tr>
	   <tr><td align="center" colspan='3'>
	   <input type="button" style="width:80px" value="Add All" name="addall" onclick="fill_species_ali(genome.options[genome.selectedIndex].value,'species_ali_keep','species_ali')"/>
	   <input type="button" style="width:80px" value="Deselect all" name="deselectall" onclick="deselect_all('species_ali');deselect_all('species_ali_keep');"/>
	   <input type="button" style="width:80px" value="Remove All" name="remove_all" onclick="fill_species_ali(genome.options[genome.selectedIndex].value,'species_ali','species_ali_keep')"/>
	   </td></tr>					  
	   </table>
	   </p>
	   </fieldset><br/>					

	   <fieldset>
	   <legend><b>Motifs</b></legend>
	   <p>
	   <b>ID of Reference Motif</b> 
	   (<span style='color:orange'>must belong to the motif database</span>, e.g. MA0102.2)
	   <br/><input type="text" name="r_motif">
	   </p>
	   
	   <p>
	   <b>Custom motif(s)</b>  (Optional)
	   <br/>Should be provided in transfac format, in any of the two following ways:				
	   <ul style="list-style-type:square">
	   <li>Paste matrices<br/><textarea name='custom_motifs' rows='6' cols='45'></textarea></li>
	   <li>Upload a file from your computer<br/><input type='file' name='custom_motifs_file' id='custom_motifs_file' size='40' /></li>
	   </ul>
	   </p>					
	   </fieldset><br/>
	   
	   <div class="menu_heading_closed" onclick="toggleMenu('101')" id="heading101"><b>Advanced options</b></div><br/>
     
	     <div id="menu101" class="menu_collapsible">
	       <fieldset>
					<table>
						<tr><th>Maximun numbers of peaks ]1..1000]<br/>
							<span style="font-weight:normal;">The numbers of peaks is limited to 1000 peaks.<br/>
							For biggest data use <a target='_blank' href='<?php echo $properties['rsat_www'].'distrib/' ?>'>stand-alone version</a> of peak-footprints.</span>
							</th><th><input type="text" size="3" name="nb_peaks" value="1000"/></th></tr>
						<tr><th>Maximun chi2 Pvalue (>0)</th><th><input type="text" size="3" name="max_chi2_pvalue" value="1"/></th></tr>
						<tr><th>Maximun Hypergeometric P-value (>0)</th><th><input type="text" size="3" name="max_hyp_pvalue" value="0.001"/></th></tr>		
						<tr><th>Column-wise conservation threshold ]0..1[</th><th><input type="text" size="3" value="0.7" name="cons_thres"/></th></tr>
						<tr><th>Block-wise conservation threshold ]0..1[</th><th><input type="text" size="3" value="0.7" name="window_cons_thres"/></th></tr>	
						<tr><th>Initial sliding window size (>1)</th><th><input type="text" size="3" value="5" name="window_size"/></th></tr>
						<tr><th>Maximum number of reported motifs (>1)</th><th><input type="text" size="3" value="50" name="motif_number"/></th></tr>
						<tr><th>Maximum number of reported motifs per family (>1)</th><th><input type="text" size="3" value="4" name="motif_number_family"/></th></tr>
					</table>
					
	       </fieldset>
	     </div>
        <b>Output</b>&nbsp;<input type="radio" name="output" value="display" />display <input type="radio" name="output" value="email" checked="checked" />email <input type="text" name="user_email"  size="30" />
      <ul><table class='formbutton'>
        <tr style="valign:middle">
          <td><input type="submit" name="submit" value="GO" /></td>
          <td><input type="reset"  name="reset"/></td> 
<?php
if ( array_key_exists('mm9', $multiz_supported)) {
	$demo_button = '<td><input type="button" name="demo" value="Demo (mm9)" onclick="'."add_demo('mm9', Array ('hg18', 'rn4', 'bosTau3', 'monDom4'), '";
	$demo_button .= $properties['rsat_www']."demo_files/SWEMBL_mmus_CEBPA_vs_mmus_Input_peaks_R0.12_nof.bed";
	$demo_button .= "')\"/></td>";
	echo $demo_button;
}

?>
 <!--         <td><b><a href='help.fetch-sequences.html'>[MANUAL]</a></b></td>  -->
          <td><b><a href='http://www.bigre.ulb.ac.be/forums/' target='_top'>[ASK A QUESTION]</a></b></td>
        </tr></table>
      </ul>
		</form>
	</body>
</html>
