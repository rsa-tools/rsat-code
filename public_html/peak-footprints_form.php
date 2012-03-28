<?php
## Load RSAT configuration
require ('functions.php');
#print_r($properties);
UpdateLogFile("rsat","","");
?>	
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
	<head>
		<title>RSAT - peak-footprints</title>
		<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
		<link rel="stylesheet" type="text/css" href="main.css" media="screen,projection,print"/>
		<link rel="stylesheet" type="text/css" href="tabs.css" media="screen,projection,print"/>
		<script src="RSAT_menu.js" type="text/javascript"></script>
		<script src="RSAT_tabs.js" type="text/javascript"></script>
		<script type="text/javascript">
			
		var multiz_genome = new Array();
		multiz_genome[0] = new Array();
<?php
//////////////////////////
//Get info on avalaible genome on ucsc
$file = fopen ($properties['RSAT'].'/public_html/data/available_organism_ucsc.tab', 'r');
$genome_ucsc= array();

while (!feof($file)) {
	list($genome, $name, $scientist_name, $date) = explode("\t", substr(fgets($file),0,-1));

	if ($genome != "") {
		$genome_ucsc[$genome]['name'] = $name;
		$genome_ucsc[$genome]['scientist_name'] = $scientist_name;
		$genome_ucsc[$genome]['date'] = $date;
	}
}
fclose(file);
#print_r($genome_ucsc);


//Get supported multiz alignement and genone in each multiz
$file = fopen ($properties['RSAT'].'/public_html/data/supported_organism_ucsc_multiz_local.tab', 'r');
$multiz_supported = array();

while (!feof($file)) {
	list($genome, $genome_reference) = explode("\t", substr(fgets($file),0,-1));

	if ($genome_reference != "") {
		$multiz_supported[$genome_reference][] = $genome;
	}
}
fclose(file);
#print_r($multiz_supported);

//Write array() with key multiz available genome and values genome in multiz.
foreach($multiz_supported as $key => $value) {
	echo "multiz_genome['$key'] = new Array(";
	
	foreach($value as $key1 => $value1) {
		if ($value1 == end($value)) {
			echo "'$value1 - ".$genome_ucsc[$value1]['name']."'";
		}	else {
			echo "'$value1 - ".$genome_ucsc[$value1]['name']."',";
		}
	}
	echo ")\n";
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
					document.getElementById(fill_area).options[i].text = genomes[i];
				}		
			}	
		}

		function deselect_all(id) {
			document.getElementById(id).options.selectedIndex = 0;
			document.getElementById(id).options[0].selected = false;
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
		</script>
	</head>

	<body class="form">
		<div>
			<h3 align='center'><a href="<?php echo $properties['rsat_www']?>">RSA-tools</a> - peak-footprints</h3><br/>
			<form name="form" method='post' action='test.php' enctype='multipart/form-data' onreset="on_reset('species_ali');on_reset('species_ali_keep');" onsubmit="select_all('species_ali_keep')">
				<fieldset>  
					<legend><b>??</b></legend>
					<p>
						<b>Genomic coordinates</b> <span style='color:red'>(mandatory)</span>
						<br/>Should be provided as a bed file (<a target='_blank' href='http://genome.ucsc.edu/FAQ/FAQformat.html#format1'>bed format</a>), in any of the two following ways:				
						<ul type='square'>
							<li>Paste coordinates<br/><textarea name='bed' rows='6' cols='45'></textarea></li>
							<!--  <li>Specify the URL of bed file on remorte server (e.g. Galaxy)<br/><input type="text" name="sequence_url" size="62" /><br/></li> -->
							<li>Upload a file from your computer<br/><input type='file' name='bedfile' id='bedfile' size='40' /></li>
						</ul>
					</p>
					<p>
		         <b>UCSC genome</b> <span style='color:red'>(mandatory)</span>&nbsp;&nbsp;			
						<select name="genome" onChange="fill_species_ali(this.options[this.selectedIndex].value,'species_ali','species_ali_keep');">
							<option value="0">----Choose a genome----</option>								
<?php

foreach($multiz_supported as $key => $value) {
	echo "<option value='$key'>$key - ".$genome_ucsc[$key]['name']."</option>", "\n";
}
?>
						</select>
					</p>							
					<p>
					  <b>Aligned species</b> <span style='color:red'>(mandatory)</span><br/>
					  <table>
						  <tr>
							  <td><select name="species_ali[]" id="species_ali" style="width:200px" size="8" multiple></select></td>
							  <td><input type="button" value="Add >>" name="add" onclick="move('species_ali','species_ali_keep');" style="width:80px"><br/><input type="button" style="width:80px" value="&lt;&lt; Remove" name="remove" onclick="move('species_ali_keep','species_ali');" size='10'></td>
							  <td><select name="species_ali_keep[]" id="species_ali_keep" style="width:200px" size="8" multiple></select></td>
						  </tr>
						  <tr><td align="center" colspan='3'>
							  <input type="button" style="width:80px" value="Add All" name="addall" onclick="fill_species_ali(genome.options[genome.selectedIndex].value,'species_ali_keep','species_ali')">
								<input type="button" style="width:80px" value="Deselect all" name="deselectall" onclick="deselect_all('species_ali');deselect_all('species_ali_keep');">
								<input type="button" style="width:80px" value="Remove All" name="remove_all" onclick="fill_species_ali(genome.options[genome.selectedIndex].value,'species_ali','species_ali_keep')">
							</td></tr>					  
					  </table>
					</p>
				</fieldset><br/>					
				<fieldset>
	         <legend><b>??</b></legend>
					<p>
					  <b>Reference Motif</b> <span style='color:red'>(mandatory)</span><br/>	
						<input type="text" name="r_motif">
					</p>
					<p>
						<b>Motif TF</b>
						<br/>Should be provided as a transfac file, in any of the two following ways:				
						<ul type='square'>
							<li>Paste coordinates<br/><textarea name='bed' rows='6' cols='45'></textarea></li>
							<li>Upload a file from your computer<br/><input type='file' name='bedfile' id='bedfile' size='40' /></li>
						</ul>
					</p>					
				</fieldset><br/>
				
				<div class="menu_heading_closed" onclick="toggleMenu('101')" id="heading101"><b>Advanced option.</b></div><br/>
      
	      <div id="menu101" class="menu_collapsible">
	        <fieldset>
	          <legend><b>??</b></legend>
						<table>
							<tr><th>??</th><th><input type="text" size="3" value="0,7" name="s_position"></th></tr>
							<tr><th>??</th><th><input type="text" size="3" value="0,7" name="s_window"></th></tr>
							<tr><th>??</th><th><input type="text" size="3" value="5" name="t_window"></th></tr>
							<tr><th>Max number of motif in output</th><th><input type="text" size="3" value="50" name="output_motif"></th></tr>
							<tr><th>Max number of motif per family</th><th><input type="text" size="3" value="4" name="motiF_family"></th></tr>																														
						</table>
	        </fieldset>
	      </div>
	  	      
        <b>Output</b>&nbsp;<input type="radio" name="output" value="display" />display <input type="radio" name="output" value="email" checked="checked" />email <input type="text" name="user_email"  size="30" />
         
	      <ul><table class='formbutton'>
	        <tr valign=middle>
	          <td><input type="submit" name="submit" value="GO" /></td>
	          <td><input type="reset"  name="reset" onclick="deselectall();"/></td> 
<!--          <td><input type="button" name="demo" value="Demo" onclick="add_demo()"/></td>  
	          <td><b><a href='help.fetch-sequences.html'>[MANUAL]</a></b></td>  -->
	          <td><b><a href='http://www.bigre.ulb.ac.be/forums/' target='_top'>[ASK A QUESTION]</a></b></td>
	        </tr></table>
	      </ul>
			</form>
		</div>
	</body>
</html>