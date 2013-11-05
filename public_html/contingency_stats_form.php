<html>
<head>
   <title>Network Analysis Tools - contingency-stats</title>
   <link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
</head>
<body class="form">
<?php
  require ('functions.php');
  require ('demo_dataset.php');
  # variable definition
  $default_stats = "checked";
  $default_rowstats = "checked";
  $default_colstats = "checked";
  
  # PIPE VALUES
  $pipe = $_REQUEST['pipe'];
  $matrix_file = $_REQUEST['matrix_file'];
  
  # demo graph
  $demo = $_REQUEST['demo'];
  if ($demo == 1) {
    $demo_matrix = storeFile("demo_files/gavin_mcl_clusters_inf2.1_cc_mips_complexes_matrix.tab");
    $pipe = 0;
    $demo_comment = "This contingency table consists in the comparaison between <a href = 'http://micans.org/mcl/' target = 'top'>MCL</a> obtained clusters applied on the <a target = '_blank' href = 'http://www.ncbi.nlm.nih.gov/sites/entrez?Db=pubmed&Cmd=ShowDetailView&TermToSearch=16429126&ordinalpos=1&itool=EntrezSystem2.PEntrez.Pubmed.Pubmed_ResultsPanel.Pubmed_RVDocSum'>Gavin <i>et al</i> (2006)</a> co-immunoprecipitation dataset (inflation 1.8) and the complexes described in the <a target = '_blank' href = 'http://mips.gsf.de/'>MIPS database</a>.";
  }

   title('contingency-stats');
   echo ("<center>This programs takes as input a contingency table, and calculates various
	  matching statistics between the rows and columns. <br>The description of
	  these statistics can be found in <a href = 'http://www.biomedcentral.com/1471-2105/7/488'>Broh&eacute;e and van Helden (2006)</a>\n");
   echo ("<br>This program was developed by <a target=_blank href=http://www.bigre.ulb.ac.be/Users/jvanheld/>Jacques van Helden</a>.
	  </center>\n");


  echo ("<form method='post' action='contingency_stats.php' enctype='multipart/form-data'>");
  
  echo("<br>");
  if (!$pipe) {
    if ($demo) {
      demo($demo_comment);
    }
    echo ("<b>Contingency table</b><br>");
    echo ("<textarea name='matrix' rows='6' cols='65'>$demo_matrix</textarea>
    <br>Upload contingency table from file : <br>
    <input type='file' name='matrix_file' size='45' /><br>");
  } else {
    info_link("Contingency table loaded from the previous treatment", rsat_path_to_url($matrix_file));
    echo "<input type='hidden' NAME='pipe_matrix_file' VALUE='$matrix_file'>";
  }

  echo("<table>
  <tr><td  colspan = 2><b><a href = 'help.contingency_stats.html#return_fields'>Return fields</a></b></td></tr>
  <tr><td><input type='checkbox' name='stats' value='on' $default_stats/></td><td><B>table-wise statistics</B></td>
  <tr><td><input type='checkbox' name='rowstats' value='on' $default_rowstats/></td><td><B>row-wise statistics</B></td>
  <tr><td><input type='checkbox' name='colstats' value='on' $default_colstats/></td><td><B>column-wise statistics</B></td>
  <tr><td><input type='checkbox' name='tables' value='on' /></td><td><B>full tables for each statistics</B></td>
  <tr><td><input type='checkbox' name='margins' value='on' /></td><td><B>marginal statistics besides the tables</B></td>
  </table>
  
  
  <ul><ul><table class='formbutton'>
  <TD><input type='submit' name='.submit' value='GO' /></TD>
  <TD><B><A HREF='contingency_stats_form.php?demo=0'>RESET</A></B></TD>
  <TD><B><A HREF='contingency_stats_form.php?demo=1'>DEMO</A></B></TD>");
  echo("</form>
  <TD><B><A HREF='help.contingency_stats.html'>MANUAL</A></B></TD>
  <TD><B><A target = '_blank' HREF='".checkNeatTutorial("tutorials/neat_tutorial/Influence_graph_alteration.html")."'>TUTORIAL</A></B></TD>

  <TD><B><A HREF='mailto:sylvain@bigre.ulb.ac.be'>MAIL</A></B></TD>
  </TR></TABLE></ul></ul>
 ");

?>
</body>
</html>
