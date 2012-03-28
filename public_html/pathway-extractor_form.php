<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
  <head>
    <title>RSAT - fetch-sequences</title>
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
    <link rel="stylesheet" type="text/css" href="http://rsat.ulb.ac.be/rsat//main.css" media="screen,projection,print"/>
    <link rel="stylesheet" type="text/css" href="http://rsat.ulb.ac.be/rsat//tabs.css" media="screen,projection,print"/>
    <script src="RSAT_menu.js" type="text/javascript"></script>
    <script src="RSAT_tabs.js" type="text/javascript"></script>
    <script type="text/javascript">
			function add_demo() {
				document.getElementById('gnn').value = 'Escherichia_coli_strain_K12-83333-Uniprot'; 
				document.getElementById('network').value = 'MetaCyc_v141_directed'; 	
				document.getElementById('seeds').innerHTML="NP_416523.1\nNP_416524.1\nNP_416525.1\nNP_416526.4\nNP_416527.1\nNP_416528.2\nNP_416529.1\nNP_416530.1";
			}
	</script>												
<?php

## Load RSAT configuration
require ('functions.php');
#print_r($properties);#
UpdateLogFile("rsat","","");

#echo "document.forms[0].sequence_url.value = '".$properties['rsat_www']."/demo_files/fetch-sequences_Schmidt_2011_mm9_CEBPA_SWEMBL_R0.12_702peaks.bed'";

?>
		
	
  </head>

  <body class="form">
    <div>
      <h3 align='center'><a href="<?php echo $properties['rsat_www']?>">RSA-tools</a> - Pathway Extractor</h3>
      <br/>
      <form method='post' action='pathway-extractor.php' enctype='multipart/form-data'>

        <fieldset>  
         <legend><b>Select a Genome</b></legend>    
         <b>Gene Reaction mapping</b> <font color='red'>(mandatory)</font>&nbsp;&nbsp;&nbsp;
          <select name='gnn' id='gnn'>
           <option value ='none'> --- Gene Reaction --- </option>
<?php

   ## Get the list of supported organisms for metabolic analysis and display it in the pop-up menu
   $cmd = $properties['RSAT'].'/perl-scripts/supported-organisms-metab';
   exec($cmd, $metab_organisms);
   sort($metab_organisms);
   

   foreach ($metab_organisms as $ligne) {
     list($mapping) =  explode("\t", $ligne);
     echo "<option value = '$mapping'>", $mapping, "</option>\n";
   }
   
?>
        </select><br/><br/>
<legend><b>Select a Network</b></legend>    
         <b>Netwoks </b> <font color='red'>(mandatory)</font>&nbsp;&nbsp;&nbsp;
          <select name='network' id='network'>
           <option value ='none'> ---Networks--- </option>
<?php

   ## Get the list of supported organisms from METAB and display it in the pop-up menu
   $cmd = $properties['RSAT'].'/perl-scripts/supported-metabolic-networks';
   exec($cmd, $metab_networks);
   sort($metab_networks);
   foreach ($metab_networks as $ligne) {
     list($network, $source) =  explode("\t", $ligne);
     echo "<option value = '$network'>", $network, "</option>\n";
   }
   
?>
        </select><br/><br/>
    <p>
	  <b><a href='help.pathway-extractor.html#input_format'>input format</a></b> <font color='red'>(mandatory)</font>
	  <br/>should be provided as a bed file (<a target='_blank'
	  href='http://genome.metab.edu/FAQ/FAQformat.html#format1'>bed
	  format</a>), in any of the three following ways:

	  <ul type='square'>
	    <li>Genes,Compounds,ECs,Reactions<br/>
	          <textarea name='seeds' id='seeds'rows='6' cols='45'></textarea>
       	</li>

	  </ul>
	  </p>
      </fieldset><br/>
      
      
      <b>Output</b>&nbsp;<input type="radio" name="output" value="display" checked="checked" />display <input type="radio" name="output" value="email" />email <input type="text" name="user_email"  size="30" />

      <ul><table class='formbutton'>
        <tr valign=middle>
          <td><input type="submit" name="submit" value="GO" /></td>
          <td><input type="reset"  name="reset" /></td> 
          <td><input type="button" name="demo" value="Demo" onclick="add_demo()"/></td>  
          <td><b><a href='help.fetch-sequences.html'>[MANUAL]</a></b></td>
          <td><b><a href='http://www.bigre.ulb.ac.be/forums/' target='_top'>[ASK A QUESTION]</a></b></td>
        </tr></table>
      </ul>
    </div>
  </body>
</html>
