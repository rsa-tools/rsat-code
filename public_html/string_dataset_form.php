<html>
<head>
   <title>Network Analysis Tools - download string dataset</title>
   <link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
</head>
<body class="form">
<?php
  require ('functions.php');
  require ('string_functions.php');
  title('String dataset download');
  
  $evidence = array('combined_confidence','automated_textmining','experimental_interaction_data','gene_cooccurence', 'gene_fusion_events','gene_coexpression');
  # demo genes
  $demo = $_REQUEST['demo'];
  if ($demo == 1) {
    $demo_genes = "GAP1\nMEP2\nMEP1";
  }  
  
   echo ("<center>Download a subgraph from the <a href = 'http://string.embl.de/'>STRING database</a></center>\n");
   echo ("<br>
   <a href = 'http://string.embl.de/'>STRING</a> is a database of known and predicted protein-protein interactions.
The interactions include direct (physical) and indirect (functional) associations. <a href = 'http://string.embl.de/'>STRING</a> quantitatively integrates interaction data from various sources for a large number of organisms, and transfers information between these organisms where applicable. The database currently contains 1,513,782 proteins in 373 species.<br> This program allows you to download a subgraph composed of the neighbours of some query genes in the <a href = 'http://string.embl.de/'>STRING database</a>. 
For more information on <a href = 'http://string.embl.de/'>STRING</a>, see <a href = 'http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=17098935&dopt=Abstract'>Von Mering et al, 2007</a>          \n");
   
   echo ("<form method='post' action='string_names.php' enctype='multipart/form-data'> ");
   
  # LIST OF ORGANISM
  echo("<br>&nbsp;&nbsp;&nbsp;<a href = 'help.string_dataset.html#organism'><B>Organism</B></a><br>");
  $organisms = readStringOrganisms();
  echo ("<select name='organism'>");
  for ($i = 1; $i < count($organisms); $i++) {
    $selected = "";
    if ($organisms[$i][1] == 'Saccharomyces cerevisiae') {
      $selected = "selected";
    } 
    echo("<option $selected value = '".$organisms[$i][0]."'> ".$organisms[$i][1]);
  }
  echo ("</select>");
  
  # LIST OF PROTEINS
  echo ("<br><br>&nbsp;&nbsp;&nbsp;<b><a href = 'help.string_dataset.html#genes'>Genes</a></b><br>");
  echo ("<textarea name='genes' rows='6' cols='65'>$demo_genes</textarea><br>");
  
  # SIGNIFICANCE THRESHOLD
  echo ("<br><br>&nbsp;&nbsp;&nbsp;<b><a href = 'help.string_dataset.html#threshold'>String score threshold</a> (between 0 and 1)</b><br>");
  echo("  
    <table>
    <tr><td><B>Lower threshold</B></td><td><input type = 'text' value = 'none' name='lth' size = 4></input></td></tr>
    <tr><td><B>Upper threshold</B></td><td><input type = 'text' value = 'none' name='uth' size = 4></input></td></tr>
    <table>");
  
  # EVIDENCE CHANNELS
  echo ("<br><br>&nbsp;&nbsp;&nbsp;<b><a href = 'help.string_dataset.html#channels'>Evidence channel</a></b><br><br>");

    for ($i = 0; $i < count ($evidence); $i++) {
      $evidence_channel = $evidence[$i];
      $evidence_channel = str_replace("_", " ", $evidence_channel);
      $checked = "";
      if ($evidence[$i] == "combined_confidence") {
        $checked = "checked";
      }
      echo(" <input type='checkbox' $checked name='$evidence[$i]' value='on'/>&nbsp;<B>$evidence_channel</B><br>");
    }

  
  
echo("  
  <ul><ul><table class='formbutton'>
  <TD><input type='submit' name='.submit' value='GO' /></TD>
  <TD><B><A HREF='string_dataset_form.php?demo=0'>RESET</A></B></TD>
  <TD><B><A HREF='string_dataset_form.php?demo=1'>DEMO</A></B></TD>
  </form>

  <TD><B><A target = '_blank' HREF='".checkNeatTutorial("")."'>TUTORIAL</A></B></TD>
  <TD><B><A target = '_blank' HREF='help.string_dataset.html'>HELP</A></B></TD>
  <TD><B><A HREF='mailto:sylvain@bigre.ulb.ac.be'>MAIL</A></B></TD>
  </TR></TABLE></ul></ul>
 ");

?>
</body>
</html>
