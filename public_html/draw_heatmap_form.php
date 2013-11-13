<html>
<head>
   <title>Network Analysis Tools - draw-heatmap</title>
   <link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
</head>
<body class="form">
<?php
   require ('functions.php');

# variable definition
$default_rowh = 10;
$default_colw = 10;
$default_min = "";
$default_max = "";
$default_gradient = " selected ";
  
# PIPE VALUES
$pipe = $_REQUEST['pipe'];
$query_table_file = $_REQUEST['table_file'];
  
# demo classes
$demo = $_REQUEST['demo'];
$demo_output_format = "selected";
if ($demo == 1) {
  $demo_table = storeFile("demo_files/krogan_2006_mcl_vs_mips_table.tab");
 }

title('draw-heatmap');

## Description
echo ("<center>Draw a heatmap from a table <br>");

## Authors
echo ("This program was developed by
	  <a target=_blank
          href='http://www.bigre.ulb.ac.be/~sylvain/'>Sylvain
          Broh&eacute;e</a> and <a target=_blank
          href=http://www.bigre.ulb.ac.be/Users/jvanheld/>Jacques van
          Helden</a>.
          </center>\n");

## Open the form
echo ("<form method='post' action='draw_heatmap.php'
	 enctype='multipart/form-data'>");

################################################################
## Input options
echo("<hr>");
echo("<p><b>Input</b></p>");

## Input format   

if ($pipe) {
    info_link("Table uploaded from the previous treatment", rsat_path_to_url($query_table_file));
    echo "<input type='hidden' NAME='pipe_query_table_file' VALUE='$query_table_file'/>";
} else {
  if ($demo) {
    demo("This demonstration consists in a contingency table that results from the comparison of the complexes annotated in the MIPS database to clusters returned by the MCL algorithm when applied to the yeast interaction dataset produced by Krogan et al (2006). Each cell contains the Jaccard coefficient (intersection size / union size) of a cluster - complex comparison.");
  }

  ## Classes
  echo ("<b>Table</b><br>");

  echo ("<textarea name='table' rows='6' cols='65'>$demo_table</textarea>
    <br>Upload table from file : <br>
    <input type='file' name='table_file' size='45' /><br>
    &nbsp;&nbsp;&nbsp;");
  

 }




################################################################
## Output options

echo("<hr>");
echo("<p><b>Output options</b></p>");
echo("
    <br>
    <table>
    <tr><td><B><a href = 'help.draw_heatmap.html#rowh'>Row height</a></B></td><td><input type = 'text' name='rowh' value = '$default_rowh' size = 2></input></td></tr>
    <tr><td><B><a href = 'help.draw_heatmap.html#colw'>Col width</a></B></td><td><input type = 'text' name='colw' value = '$default_colw' size = 3></input></td></tr>
    <tr><td><B><a href = 'help.draw_heatmap.html#min'>Minimum value</a></B></td><td><input type = 'text' name='min' size = 2 value = '$default_min'></input></td></tr>
    <tr><td><B><a href = 'help.draw_heatmap.html#max'>Maximum value</a></B></td><td><input type = 'text' name='max' size = 2 value = '$default_max'></input></td></tr>
    </table>");


## Output format   
echo("<B><a href = 'help.draw_heatmap.html#formats'>Output format</a></B>&nbsp;<select name='out_format'>");
echo ("<option selected value = 'png'>png
<option value = 'jpg'>jpeg
</select><br>");

## Color gradient
echo("<B><a href = 'help.draw_heatmap.html#gradient'>Color gradient</a></B>&nbsp;<select name='gradient'>");
echo ("<option selected value = 'fire'>fire
<option selected value = 'grey'>grey
<option value = 'blue'>blue
<option value = 'green'>green
<option value = 'red'>red
</select><br><br><br>");

# Text in the cells
echo ("<input type='checkbox' value = '1' name='no_text'/>&nbsp;<B><a href = 'help.draw_heatmap.html#no_text'>Disable cell values.</a></B><br><br>");

# The fist column of the table contains row names
echo ("<input type='checkbox' value = '1' checked name='rownames'/>&nbsp;<B><a href = 'help.draw_heatmap.html#rownames'>The first column contains the names of the row.</a></B><br><br>");



## Action buttons
echo ("<hr>");
echo ("<br>
  <table class='formbutton'>
  <TD><input type='submit' name='.submit' value='GO' /></TD>
  <TD><B><A HREF='draw_heatmap_form.php?demo=0'>RESET</A></B></TD>
  <TD><B><A HREF='draw_heatmap_form.php?demo=1'>DEMO</A></B></TD>
  </form>
  <TD><B><A HREF='help.draw_heatmap.html'>MANUAL</A></B></TD>
  <TD><B><A target = '_blank' HREF='".checkNeatTutorial("tutorials/neat_tutorial/Network_visualization_forma.html")."'>TUTORIAL</A></B></TD>
  <TD><B><A HREF='mailto:sylvain@bigre.ulb.ac.be'>MAIL</A></B></TD>
  </TR></TABLE>
 ");

?>
</body>
</html>
