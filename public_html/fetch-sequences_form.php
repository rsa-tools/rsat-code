<html>
  <head>
    <title>Network Analysis Tools - compare-graphs</title>
    <META HTTP-EQUIV="Content-Type" CONTENT="text/html; charset=UTF-8">
    <link rel="stylesheet" type="text/css" href="http://rsat.ulb.ac.be/rsat//main.css" media="screen,projection,print"/>
    <link rel="stylesheet" type="text/css" href="http://rsat.ulb.ac.be/rsat//tabs.css" media="screen,projection,print"/>
    <script src="RSAT_menu.js" type="text/javascript"></script>
    <script src="RSAT_tabs.js" type="text/javascript"></script>
    <script language="javascript">
			function add_demo() {
<?php
## Load RSAT configuration
require ('functions.php');
#   print_r($properties);
UpdateLogFile("rsat","","");

echo "document.forms[0].sequence_url.value = '".$properties['RSAT']."/public_html/demo_files/fetch-sequences_Schmidt_2011_mm9_CEBPA_SWEMBL_R0.12_702peaks.bed'";
?>
			}
		</script>
  </head>

  <body class="form">
    <h3 align='center'><a href='http://rsat.ulb.ac.be/rsat//RSAT_home.cgi'>RSA-tools</a> - fetch-sequence</h3>
    <br/>
    <form method='post' action='fetch-sequences.php' enctype='multipart/form-data'>

      <fieldset>  
       <legend><b>Genomic coordinates</b></legend>    
       <b>Genome </b> <font color='red'>(mandatory)</font>&nbsp;&nbsp;&nbsp;
        <select name='genome'>
         <option value = 'none'> ---UCSC genome---
<?php
   ## Get the list of supported organisms from UCSC and display it in the pop-up menu
   $cmd = $properties['RSAT'].'/perl-scripts/supported-organisms-ucsc';
   exec($cmd, $ucsc_organisms);
   sort($ucsc_organisms);
   foreach ($ucsc_organisms as $ligne) {
   list($genome,$description) =  explode("\t", $ligne);
   echo "<option value = '$genome'>", $genome,"\n";
   }
   ?>
        </select><br/><br/>

        <p>
	  <b><a href='help.fetch-sequences.html#input_format'>Genomic
	  coordinates</a></b> <font color='red'>(mandatory)</font>
	  <br/>should be provided as a bed file (<a target='_blank'
	  href='http://genome.ucsc.edu/FAQ/FAQformat.html#format1'>bed
	  format</a>), in any of the three following ways:

	<ul type='square'>
	  <li>Paste coordinates<br/>
	    
            <textarea name='bed' rows='6' cols='45'></textarea></li>

          <li>Specify the URL of bed file on remorte server (e.g. Galaxy)<br/>
            <input type="text" name="sequence_url" size="62" /><br/></li>

          <li>Upload a file from your computer<br/>
	    <input type='file' name='bedfile' size='40' />
	  </p>
	</ul>

	    
          <p><b><a href='help.fetch-sequences.html#header'>Header Format</a></b><input type="radio" name="header" value="ucsc" checked="checked"/>UCSC
            <input type="radio" name="header" value="galaxy"/>Galaxy     
      </fieldset><br/>
      
      <div class="menu_heading_closed" onclick="toggleMenu('101')" id="heading101">
        <span><b>Reference from which the sequences should be fetched.</b></span></div><br/>
      
      <div id="menu101" class="menu_collapsible">
        <fieldset>
          <legend><b><a href='help.fetch-sequences.html#reference'>Reference</a></b></legend>
          <input type="radio" name="reference" value="segment" checked="checked"/><b>Segment</b>
          <input type="radio" name="reference" value="start"/><b>Start</b>
          <input type="radio" name="reference" value="end"/><b>End</b>
        </fieldset><br/>
      </div>

      <div class="menu_heading_closed" onclick="toggleMenu('102')" id="heading102">
          <span><b>Extend sequence retrieve</b></span></div><br/>

      <div id="menu102" class="menu_collapsible">
        <fieldset>
          <table border='0'>
            <legend><b><a href='help.fetch-sequences.html#options'>Extend</a></b></legend>
            <!--  <tr align='center'><td colspan='4'>OR</td></tr>-->
            <tr><td><b><a href='help.fetch-sequences.html#upstr_ext'>Upstream side</a></b></td><td style='padding-right:15px;'><input type="input" name="upstr_ext" size ='3'/>&nbsp;bp</td>
              <!--  <td rowspan="2" align='center' style='border-left:1px solid #2D282E;' width='120px'><b><a href='help.fetch-sequences.html#extend'>Both side</a></b></td>
                <td rowspan="2" ><input type="input" name="extend" size = '3'/>&nbsp;bp </td>
            </tr>
            <tr>--><td><b><a href='help.fetch-sequences.html#downstr_ext'>Downstream side</a></b></td><td style='padding-right:15px;'><input type="input" name="downstr_ext" size = '3'/>&nbsp;bp</td></tr>
          </table>
        </fieldset><br/>
      </div>
      
      <b>Output</b>&nbsp;<input type="radio" name="output" value="display" />display <input type="radio" name="output" value="email" checked="checked" />email <input type="text" name="user_email"  size="30" />

      <ul><ul><table class='formbutton'>
        <tr valign=middle>
          <td><input type="submit" name="submit" value="GO" /></td>
          <td><input type="reset"  name="reset" /></td> 
          <td><input type="button" name="demo" value="Demo" onclick="add_demo()"/></td>  
          <td><b><a href='help.fetch-sequences.html'>[MANUAL]</a></b></td>
          <td><b><a href='tutorials/tut_peak-motifs.html'>[TUTORIAL]</a></b></td>
          <td><b><a href='http://www.bigre.ulb.ac.be/forums/' target='_top'>[ASK A QUESTION]</a></b></td>
        </tr></table>
      </ul></ul>

  </body>
</html>
