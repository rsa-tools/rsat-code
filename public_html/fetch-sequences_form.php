<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
  <head>
    <title>RSAT - fetch-sequences</title>
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
    <link rel="stylesheet" type="text/css" href="main.css" media="screen,projection,print"/>
    <link rel="stylesheet" type="text/css" href="tabs.css" media="screen,projection,print"/>
    <script src="RSAT_menu.js" type="text/javascript"></script>
    <script src="RSAT_tabs.js" type="text/javascript"></script>
    <script type="text/javascript">
			function add_demo() {
				document.getElementById('genome').value = 'mm9'; 	
				document.forms[0].bed.value="";
				document.getElementById('bedfile').value=''; 				
				document.forms[0].header[0].checked=true;			
				document.forms[0].reference[0].checked=true;
				document.forms[0].upstr_ext.value="0";	
				document.forms[0].downstr_ext.value="0";
				document.forms[0].output[0].checked=true;									
<?php

## Load RSAT configuration
require ('functions.php');
#print_r($properties);#
UpdateLogFile("rsat","","");

echo "document.forms[0].sequence_url.value = '".$properties['rsat_www']."/demo_files/fetch-sequences_Schmidt_2011_mm9_CEBPA_SWEMBL_R0.12_702peaks.bed'";

?>
			}
		</script>
  </head>

  <body class="form">
    <div>
      <h3 align='center'><a href="<?php echo $properties['rsat_www']?>">RSA-tools</a> - fetch-sequence</h3>
      <br/>
<?php
   ## Get the list of supported organisms from UCSC and display it in the pop-up menu
   $cmd = $properties['RSAT'].'/perl-scripts/supported-organisms-ucsc';
   exec($cmd, $ucsc_organisms);
   sort($ucsc_organisms);

  if ($ucsc_organisms=="")  {
    error("<span style='color:red'>Can't connect on UCSC. Please contact system administrator.</span>");
   }
?> 
      <form method='post' action='fetch-sequences.php' enctype='multipart/form-data'>

        <fieldset>  
         <legend><b>Genomic coordinates</b></legend>    
         <b>Genome </b> <font style='color:orange'>(mandatory)</font>&nbsp;&nbsp;&nbsp;
          <select name='genome' id='genome'>
           <option value ='none'> ---UCSC genome--- </option>
<?php
   foreach ($ucsc_organisms as $ligne) {
     list($genome,$description) =  explode("\t", $ligne);
     echo "<option value = '$genome'>", $genome, " - ", str_replace(" Genome at UCSC", "", $description), "</option>\n";
   }
?>
        </select><br/><br/>          
    <p>
	  <b><a href='help.fetch-sequences.html#input_format'>Genomic
	  coordinates</a></b> <font style='color:orange'>(mandatory)</font>
	  <br/>should be provided as a bed file (<a target='_blank'
	  href='http://genome.ucsc.edu/FAQ/FAQformat.html#format1'>bed
	  format</a>), in any of the three following ways:

	  <ul type='square'>
	    <li>Paste coordinates<br/>
	      
              <textarea name='bed' rows='6' cols='45'></textarea></li>

          <li>Specify the URL of bed file on remote server (e.g. Galaxy)<br/>
              <input type="text" name="sequence_url" size="62" /><br/></li>

            <li>Upload a file from your computer<br/>
	      <input type='file' name='bedfile' id='bedfile' size='40' /></li>
	  </ul>
	  </p>
	    
          <p><b><a href='help.fetch-sequences.html#header'>Header Format</a></b>
            <input type="radio" name="header" value="galaxy" checked="checked"/>Galaxy
	    <input type="radio" name="header" value="ucsc"/>UCSC
      </fieldset><br/>
      
      <div class="menu_heading_closed" onclick="toggleMenu('101')" id="heading101">
        <span><b>Reference from which the sequences should be fetched.</b></span></div><br/>
      
      <div id="menu101" class="menu_collapsible">
        <fieldset>
          <legend><b><a href='help.fetch-sequences.html#reference'>Reference</a></b></legend>
          <p>
            <b><a href='help.fetch-sequences.html#reference'>Reference</a></b><br/>   
            <ul>     
              <input type="radio" name="reference" value="segment" checked="checked"/><b>Segment</b>
              <input type="radio" name="reference" value="start"/><b>Start</b>
              <input type="radio" name="reference" value="center"/><b>Center</b></p>
              <input type="radio" name="reference" value="end"/><b>End</b></p>
            </ul>
          <p>
            <b><a href='help.fetch-sequences.html#options'>Extend</a></b>
            <ul>
              <li><b><a href='help.fetch-sequences.html#upstr_ext'>Upstream side&nbsp;</a></b><input type="input" name="upstr_ext" size ='3' value='0'/>&nbsp;bp</li>
              <li><b><a href='help.fetch-sequences.html#downstr_ext'>Downstream side&nbsp;</a></b><td style='padding-right:15px;'><input type="input" name="downstr_ext" size = '3' value='0'/>&nbsp;bp</li>
              </tr>
            </ul></p>  
        </fieldset><br/>
      </div>
      
      <b>Output</b>&nbsp;<input type="radio" name="output" value="display"  checked="checked"/>display <input type="radio" name="output" value="email"/>email <input type="text" name="user_email"  size="30" />

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
