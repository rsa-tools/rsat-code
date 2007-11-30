
<html>
<head>
   <title>GrA-tools - convert-graph</title>
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
    $demo_matrix = "";
  }

  title('contingency-stats');
  echo ("<center>This programs takes as input a contingency table, and calculates various
    matching statistics between the rows and columns. The description of
    these statistics can be found in <a href = 'http://www.biomedcentral.com/1471-2105/7/488'>Broh&eacute;e and van Helden (2006)</a>.</center>\n");
  echo ("<form method='post' action='contingency_stats.php' enctype='multipart/form-data'>");
  
  echo("<br>");
  if (!$pipe) {
    if ($demo) {
      demo("");
    }
    echo ("<b>Contingency table</b><br>");
    echo ("<textarea name='matrix' rows='6' cols='65'>$demo_matrix</textarea>
    <br>Upload contingency table from file : <br>
    <input type='file' name='matrix_file' size='45' /><br>");


  } else {
    info("Contingency table uploaded from the previous treatment");
    echo "<input type='hidden' NAME='pipe_matrix_file' VALUE='$matrix_file'>";
  }

  echo("<table>
  <tr><td  colspan = 2><b><a href = 'help.contingency_stats.html#return'>Return fields</a></b></td></tr>
  <tr><td><input type='checkbox' name='stats' value='on' $default_stats/></td><td><B>table-wise statistics</B></td>
  <tr><td><input type='checkbox' name='rowstats' value='on' $default_rowstats/></td><td><B>row-wise statistics</B></td>
  <tr><td><input type='checkbox' name='colstats' value='on' $default_colstats/></td><td><B>column-wise statistics</B></td>
  <tr><td><input type='checkbox' name='tables' value='on' /></td><td><B>full tables for each statistics</B></td>
  <tr><td><input type='checkbox' name='margins' value='on' /></td><td><B>marginal statistics besides the tables</B></td>
  
  
  
  <ul><ul><table class='formbutton'>
  <TD><input type='submit' name='.submit' value='GO' /></TD>
  <TD><B><A HREF='contingency_stats_form.php?demo=0'>RESET</A></B></TD>");
  // No demo at the moment
  //<TD><B><A HREF='contingency_stats_form.php?demo=1'>DEMO</A></B></TD>
  echo("</form>
  <TD><B><A HREF='help.contingency_stats.html'>MANUAL</A></B></TD>
  <TD><B><A HREF='mailto:sylvain@scmbb.ulb.ac.be'>MAIL</A></B></TD>
  </TR></TABLE></ul></ul>
 ");

?>
</body>
</html>