<html>
<head>
   <title>Network Analysis Tools - download string dataset</title>
   <link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
</head>
<body class="form">
<?php
  require ('functions.php');
  title('String dataset download');
  
  $evidence = array('automated_textmining','experimental_interaction_data','gene_cooccurence', 'gene_fusion_events','gene_coexpression', 'combined_confidence');
  # demo genes
  $demo = $_REQUEST['demo'];
  if ($demo == 1) {
    $demo_genes = "GAP1\nMEP2\nMEP1";
  }  
  
   echo ("<center>Donwload a subgraph from the String database</center>\n");
   echo ("<br>
   STRING is a database of known and predicted protein-protein interactions.
The interactions include direct (physical) and indirect (functional) associations; they are derived from four sources. STRING quantitatively integrates interaction data from these sources for a large number of organisms, and transfers information between these organisms where applicable. The database currently contains 1,513,782 proteins in 373 species.<br> This program allows you to download a subgraph composed of the neighbours of the genes of the String database.
          \n");
   
   echo ("<form method='post' action='string_dataset.php' enctype='multipart/form-data'> ");
   
  # LIST OF ORGANISM
  echo("<br>&nbsp;&nbsp;&nbsp;<B>Organism</B><br>");
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
  echo ("<br><br>&nbsp;&nbsp;&nbsp;<b>Genes</b><br>");
  echo ("<textarea name='genes' rows='6' cols='65'>$demo_genes</textarea><br>");
  
  # EVIDENCE CHANNELS
  echo ("<br><br>&nbsp;&nbsp;&nbsp;<b>Evidence channel</b><br><br>");

    for ($i = 0; $i < count ($evidence); $i++) {
      $evidence_channel = $evidence[$i];
      str_replace("_", " ", $evidence_channel);
      echo(" <input type='checkbox' checked name='$evidence[$i]' value='on'/>&nbsp;<B>$evidence_channel</B><br>");
    }

  
  
echo("  
  <ul><ul><table class='formbutton'>
  <TD><input type='submit' name='.submit' value='GO' /></TD>
  <TD><B><A HREF='string_dataset_form.php?demo=0'>RESET</A></B></TD>
  <TD><B><A HREF='string_dataset_form.php?demo=1'>DEMO</A></B></TD>
  </form>

  <TD><B><A target = '_blank' HREF='".checkNeatTutorial("")."'>TUTORIAL</A></B></TD>

  <TD><B><A HREF='mailto:sylvain@scmbb.ulb.ac.be'>MAIL</A></B></TD>
  </TR></TABLE></ul></ul>
 ");

?>
</body>
</html>
