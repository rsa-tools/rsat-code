<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
   <head>
   <title>RSAT - pathway-extractor</title>
   <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
   <link rel="stylesheet" type="text/css"
   href="main.css"
   media="screen,projection,print" />
   <link rel="stylesheet" type="text/css"
   href="tabs.css"
   media="screen,projection,print" />
   <script src="RSAT_menu.js" type="text/javascript"></script>
   <script src="RSAT_tabs.js" type="text/javascript"></script>
   <script src="pathwayextractor.js" type="text/javascript"></script>
   
   <script type="text/javascript">
   function add_demo() {
     document.getElementById('gnn').value = 'Escherichia_coli_strain_K12-83333-Uniprot'; 
     document.getElementById('network').value = 'MetaCyc_v141_directed';
     document.getElementById('seeds').innerHTML=""; 	
     document.getElementById('seeds').innerHTML="NP_416523.1\nNP_416524.1\nNP_416525.1\nNP_416526.4\nNP_416527.1\nNP_416528.2\nNP_416529.1\nNP_416530.1";
     document.getElementById('directedgraph').checked=true;
   }
   </script>
   <?php
   
   ## Load RSAT configuration
   require ('functions.php');
# print_r($properties);
#UpdateLogFile("rsat","","");
#echo "document.forms[0].sequence_url.value = '".$properties['rsat_www']."/demo_files/pathway-extractor_Schmidt_2011_mm9_CEBPA_SWEMBL_R0.12_702peaks.bed'";
   ?>
   
   </head>

   <body class="form">
   <div>
   <h3 align='center'><a href="<?php echo $properties['rsat_www']?>">RSA-tools</a>
   - Pathway Extractor</h3>
   <br />
   <form method='post' action='pathway-extractor_seeds.php'
   enctype='multipart/form-data'>

   <fieldset><legend><b>Parameters</b></legend>
   <ul>
   <li><span class="fieldtitle">Select a Genome for seeds/node mapping </span> <font color='red'>(mandatory)</font>
     <select
     name='gnn' id='gnn'>
     <option value='none'>--- Gene Reaction ---</option>
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
</select>
</li>
<li>
<span class=fieldtitle >Select a Network</span>
  <font color='red'>(mandatory)</font>&nbsp;&nbsp;&nbsp; <select name='network' id='network'>
						    <option value='none'>---Networks---</option>
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
</select>
</li>
<li>
<span class=fieldtitle >Directed Network</span><input type="checkbox" name="directedgraph"
  id="directedgraph" value="directedgraph" checked />
  </li>

  <li>
  <span class=fieldtitle>Provide seeds</span> <font color='red'>(mandatory)</font>: Genes,Compounds,ECs,Reactions
  (Minimun seed size: 2 characters)
  <ul>
  <li type="none"><textarea name='seeds' id='seeds' rows='6' cols='45'></textarea></li></ul>

  </ul>
  </p>
  </fieldset>
  <fieldset><legend class=clickable onClick="javascript:collapseAll(document.getElementById('advparams'));" >
  <b><span id=advparamsSign class=sign><img src="images/arrow_box.gif" alt="+"></img></span></b> Advanced parameters</legend>


  <div  id=advparams style="display: none">
  <table>
  <tbody>
  <tr>
  <td><a href="help.pathwayinference.jsp#algorithm"><b>Algorithm</b></a>

  <select name="algorithm">
  <option value="takahashihybrid" selected="selected">kWalks +
  Takahashi-Matsuyama</option>
  <option value="takahashi">Takahashi-Matsuyama</option>

  <option value="kwalks">kWalks</option>
  <option value="rea">iterative REA</option>
  <option value="hybrid">kWalks + iterative REA</option>
  </select></td>
  </tr>
  <tr>
  <td><a href="help.pathwayinference.jsp#preproc"><b>Preprocess</b></a>
  <input name="preproc" value="on" type="checkbox"> <font
  color="#cc6600">(for RPAIR networks only)</font></td>
  </tr>
  <tr>
  <td><a href="help.pathwayinference.jsp#postproc"><b>Postprocess</b></a>
  <input name="postproc" value="on" type="checkbox" checked></td>
  </tr>
  <tr>
  <td><a href="help.pathwayinference.jsp#postproc"><b>Iteration number</b></a>
  <input name="iter" size="10" value="1" type="text"> <font
  color="#cc6600">(for kWalks or hybrid algorithms only)</font></td>
  </tr>
  <tr>
  <td>
  <table bgcolor="#f6e6ca">
  <tbody>
  <tr>
  <td valign="top" width="200">Hybrid algorithms only</td>

  <td><a href="help.pathwayinference.jsp#hybridoptions"><b>Percentage
  of input network</b></a> <input name="percentage" size="10"
  value="0.05" type="text"></td>
  </tr>
  <tr>
  <td width="100"></td>
  <td><a href="help.pathwayinference.jsp#hybridoptions"><b>Compute
  weights with kWalks</b></a> <input name="kwalksweight" value="on"
  type="checkbox"></td>
  </tr>
  </tbody>
  </table>
  </td>
  </tr>
  </tbody>
  </table>



  </div>
  </fieldset>
  <br />


  <b>Output</b>
  &nbsp;
<input type="radio" name="output" value="display" checked="checked" />
  display
  <input type="radio" name="output" value="email" />
  email
  <input type="text" name="user_email" size="30" />

  <ul>
  <table class='formbutton'>
  <tr valign=middle>
  <td><input type="submit" name="submit" value="GO" /></td>
  <td><input type="reset" name="reset" /></td>
  <td><input type="button" name="demo" value="Demo"
  onclick="add_demo()" /></td>
  <td><b><a href='help.pathway-extractor.html'>[MANUAL]</a></b></td>
  <td><b><a href='http://www.bigre.ulb.ac.be/forums/' target='_top'>[ASK
								     A QUESTION]</a></b></td>
  </tr>
  </table>
  </ul>

  </div>
  </body>
  </html>
