<?php
    session_start();
    $_SESSION['form_uri_bed_to_matrix'] = $_SERVER['REQUEST_URI'];

?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
  <head>
    <title>RSAT - Features matrix</title>
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
    <link rel="stylesheet" type="text/css" href="css/main.css" media="screen,projection,print"/>
    <link rel="stylesheet" type="text/css" href="css/tabs.css" media="screen,projection,print"/>
    <script src="js/RSAT_menu.js" type="text/javascript"></script>
    <script src="js/RSAT_tabs.js" type="text/javascript"></script>
     <script>
     
	var selDiv = "";
		
	document.addEventListener("DOMContentLoaded", init, false);
	
	function init() {
		document.querySelector('#wig').addEventListener('change', handleFileSelect, false);
		selDiv = document.querySelector("#selectedFiles");
	}
		
	function handleFileSelect(e) {
		
		if(!e.target.files) return;
		
		selDiv.innerHTML = "";
		
		var files = e.target.files;
		for(var i=0; i<files.length; i++) {
			var f = files[i];
			
			selDiv.innerHTML += f.name + "<br/>";

		}

	}
	</script>
    <script type="text/javascript">
			
			function add_demo1() {
				
				document.getElementById('genome').value = 'dm3';
				document.getElementById('shuffling').value = '3';
				<?php
				## Load RSAT configuration
				require('functions.php');
				#print_r($properties);#
				UpdateLogFile("rsat","","");
				echo "document.forms[0].backgroundFreqs_url.value = '".$properties['rsat_www']."/demo_files/SVM_2nt_intergenic_Drosophila_melanogaster.freq'";
				echo "document.forms[0].background_seqs_url.value = '".$properties['rsat_www']."/demo_files/SVM_dmel_all_intergenic_and_introns_promoters_500bp_r5.57.bed'";
				?>

}
</script>
    <script type="text/javascript">
function add_demo2() {
				
				document.getElementById('genome').value = 'dm3';
				<?php
				echo "document.forms[0].backgroundFreqs_url.value = '".$properties['rsat_www']."demo_files/SVM_2nt_intergenic_Drosophila_melanogaster.freq'";
?>
}

function display_file(f)
{
 var httpRequest = new XMLHttpRequest();
 httpRequest.open("GET", f, true);
 httpRequest.send(null);
 httpRequest.onreadystatechange = function()
 {
  if(this.readyState == 4 && this.status == 200)
  {
   document.getElementById("matrix").value = this.responseText;
   }
  }
 }
 
 function display_filebed(f)
{
 var httpRequest = new XMLHttpRequest();
 httpRequest.open("GET", f, true);
 httpRequest.send(null);
 httpRequest.onreadystatechange = function()
 {
  if(this.readyState == 4 && this.status == 200)
  {
   document.getElementById("bed").value = this.responseText;
   }
  }
 }
 function display_filebedNeg(f)
{
 var httpRequest = new XMLHttpRequest();
 httpRequest.open("GET", f, true);
 httpRequest.send(null);
 httpRequest.onreadystatechange = function()
 {
  if(this.readyState == 4 && this.status == 200)
  {
   document.getElementById("Negbed").value = this.responseText;
   }
  }
 }
  
</script>
		
  </head>
    <body class="form">
    <div>
    <?php
    echo '
      <h3 align="center"><a href="http://pedagogix-tagc.univ-mrs.fr/rsat/">RSAT</a> - Features matrix</h3>
      '
      ?>
      <br/>
        <form method="POST" action="bed_to_matrix.php" enctype="multipart/form-data">
        <fieldset>  
         <legend><b>Genome</b></legend>    
         <b>Genome </b> <font style='color:orange'>(mandatory)</font>&nbsp;&nbsp;&nbsp;
          <select name='genome' id='genome'>
           <option value ='none'> ---UCSC genome--- </option>
           <option value = 'dm3'>dm3 - D. melanogaster Apr. 2006 (BDGP R5/dm3)</option>
			</select>
       </fieldset>  

            <fieldset>
    <legend><b>Genomic coordinates</b></legend>  
	  <b>Genomic coordinates of positive training set</b> <font style='color:orange'>(mandatory)</font>
	  <br/>should be provided as a short bed file (3 columns) (<a target='_blank'
	  href='http://genome.ucsc.edu/FAQ/FAQformat.html#format1'>bed
	  format</a>), in any of the three following ways:

	  <ul type='square'>
	    <li>Paste coordinates<br/>
	      
              <textarea name='bed' id='bed' rows='6' cols='45'></textarea></li>
	   <li>Specify the URL of bed file on remote server (e.g. Galaxy)<br/>
              <input type="text" name="sequence_url" id="sequence_url" size="62" /><br/></li>
            <li>Upload a file from your computer<br/>
	      <input type='file' name='bedfile' id='bedfile' size='40' value=''/>
	      </li>
	  </ul>
	  </p>

</fieldset>
<fieldset>
        <b>Genomic coordinates of negative training set</b> <font style='color:orange'>(mandatory)</font>
	  <br/><b>Could be provided as a short bed file (3 columns) (<a target='_blank' href='http://genome.ucsc.edu/FAQ/FAQformat.html#format1'>bed format</a>), in any of the three following ways:</b>
     
	  <ul type='square'>
	    <li>Paste coordinates<br/>
	      
              <textarea name='Negbed' id='Negbed' rows='6' cols='45'></textarea></li>
	   <li>Specify the URL of bed file on remote server (e.g. Galaxy)<br/>
              <input type="text" name="Negsequence_url" id="Negsequence_url" size="62" /><br/></li>
            <li>Upload a file from your computer<br/>
	      <input type='file' name='Negbedfile' id='Negbedfile' size='40' value='Negbedfile'/>
	      </li>
	  </ul>
	  </p>
	    
          <!-- <p><b><a href='help.fetch-sequences.html#header'>Header Format</a></b>
            <input type="radio" name="header" value="galaxy" checked="checked"/>Galaxy
	    <input type="radio" name="header" value="ucsc"/>UCSC
      <br/><br/>
<p> -->
<br/><br/>OR<br/><br/>
	<b>Generate a negative random training set</b> 
	Specify the size of the desired random training set (X times positive set) <i>(default: 2)</i><input type="text" id="shuffling" name="shuffling" size="5" value=""/><br/>
	<b>Area of sequence shuffling (e.g. intergenic, intronic, promoters..): </b><br/>
		<ul type='square'>
	     <li>Upload a bed file from your computer: <input type='file' name='background_seqs' id='background_seqs' size='40' value='background_seqs'/></li>
	     <li>Specify the URL of background sequences file on remote server (e.g. Galaxy)<input type="text" name="background_seqs_url" id="background_seqs_url" size="62"/><br/></li>
			</ul>
	     <!-- <li>
		<input type="checkbox" Name="backseqGeneric" Value ="intergenic">Intergenic<br/>
		<input type="checkbox" Name= "backseqGeneric" Value ="intronic">Intronic<br/>
		<input type="checkbox" Name= "backseqGeneric" Value ="5kbpromoter">5kb upstream genes promoter<br/>
		</li>-->


         
      </fieldset><br/>
            </div>
<p>
<fieldset>
<legend><b>Matrices</b><font style='color:orange'>(mandatory)</font></legend>
<ul type='square'>

<div><div style='float:left;'>
<li><b>Paste matrices in transfac format (<a href='help.convert-matrix.html#io_format'><b>Convert matrices</b></a>)<br/>
<textarea name="matrix" id="matrix" rows="4" cols="60" ></textarea>
</div>
<!-- <div class="menu">
<div class="menu_heading_closed"
onclick="toggleMenu('98')" id="heading98"> Where to find matrices ?</div>
<div id="menu98" class="menu_collapsible">
	<a class="menu_item" href="http://www.pazar.info/" target="_blank">PAZAR</a>
	<a class="menu_item" href="http://the_brain.bwh.harvard.edu/uniprobe/" target="_blank">UniProbe</a>
	<a class="menu_item" href="http://www.gene-regulation.com/pub/databases.html" target="_blank">Transfac</a>
	<a class="menu_item" href="http://jaspar.cgb.ki.se/" target="_blank">Jaspar</a>
	<a class="menu_item" href="http://regulondb.ccg.unam.mx/download/Data_Sets.jsp" target="_blank">RegulonDB</a>

</div>
</div>
</div>-->
</li>
<br/><br/><br/><br/><br/><br/><br/>
	      
 <li><b>Upload a matrix file (transfac format) from your computer</b><br/>
	      <input type='file' name='matrix_file' id='matrix_file' size='40' value=''/>
	      </li>
 <li>Specify the URL of matrix file on remote server (e.g. Galaxy)<br/>
              <input type="text" name="matrix_file_url" id="matrix_file_url" size="62" /><br/></li>
	  </ul>

<b>Upload a background frequencies file (<a href='help.convert-background-model.html#io_format'><b>Background frequencies file</b></a>) from your computer <font style='color:orange'>(mandatory)</font><br/>
	      <ul type='square'>
	      <li><input type='file' name='backgroundFreqs' id='backgroundFreqs' size='40' value=''/></li>
			<li>Specify the URL of background frequencies file on remote server (e.g. Galaxy)<input type="text" name="backgroundFreqs_url" id="backgroundFreqs_url" size="62" /><br/></li>
			</ul>
</fieldset><br/>
<fieldset>
       <legend> <b>NGS data</b> <font style='color:green'>(optionnal)</font></legend>
	  <br/><b>Could be provided as wig(bigwig) (<a target='_blank'
	  href='http://genome.ucsc.edu/FAQ/FAQformat.html#format1'>wig
	  format</a>) or coordinates (bed) file, in any of the following ways:</b>
     
	  <ul type='square'>
	    			<li>Specify the URL of wig(bigWig) file on remote server (e.g. Galaxy)<br/>
              <input type="text" name="wig_url" id="wig_url" size="62" /><br/></li>
            <li>Upload a file (or multiple files) from your computer<br/>
	      <input type="file" id="wig" name="wigs[]" multiple><br/>
        <div id="selectedFiles"></div>
	      </li>
	  </ul>
	  </p>
	<b>For wig(bigwig) file, select statistic to compute:</b> 
	<select name='stat' id='stat'>
           <option value ='' selected='selected'><i>Select statistic</i></option>
           <option value ='mean' >mean</option>
           <option value = 'min'>min</option>
           <option value = 'max'>max</option>
           <option value = 'total'>total</option>
           <option value = 'coverage'>coverage</option>
			</select><br/>
	<!-- <b>For bed file, specify minimum overlap:</b> <i>(default: 1bp; for 1bp write 1e-9 / for percentage, write number between 0 and 1)</i><input type="text" id="overlap" name="overlap" size="5" value=""/><br/>
	-->      
	</fieldset><br/>
            <b>Output</b>&nbsp;<input type="radio" name="output" value="display"  checked="checked"/>display <input type="radio" name="output" value="email"/>email <input type="text" name="user_email"  size="30" />

      <ul><table class='formbutton'>
        <tr valign=middle>
          <td><input type="submit" name="submit" value="GO" /></td>
          <td><input type="reset"  name="reset" /></td> 
         <?php
          
    		echo "
    		
			 <td><input type='button' name='demo' id='open' value='Demo(Random)' onclick='add_demo1();display_file('".$properties['rsat_www']."demo_files/SVM_259_matrices_small.tf);display_filebed(".$properties['rsat_www']."demo_files/SVM_heart_CC_CRM_3cols.40.bed')/></td>
			 <td><input type='button' name='demo' id='open' value='Demo(Negative bed)' onclick='add_demo2();display_file('".$properties['rsat_www']."'demo_files/SVM_259_matrices_small.tf');display_filebed('".$properties['rsat_www']."demo_files/SVM_heart_CC_CRM_3cols.40.bed');display_filebedNeg('".$properties['rsat_www']."demo_files/SVM_background_CRM_13-16.bed')/></td>
         ";
			echo $properties['rsat_www']."demo_files/SVM_259_matrices_small.tf";
         ?>
          <td><b><a href='help.fetch-sequences.html'>[MANUAL]</a></b></td>
          <td><b><a href='http://www.bigre.ulb.ac.be/forums/' target='_top'>[ASK A QUESTION]</a></b></td>
        </tr></table>
      </ul>
    </div>
   </form>
  </body>
</html>
