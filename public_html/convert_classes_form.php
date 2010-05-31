<html>
<head>
   <title>Network Analysis Tools - convert-classes</title>
   <link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
</head>
<body class="form">
<?php
   require ('functions.php');

# variable definition
$default_mcol = 1;
$default_ccol = 2;
$default_wcol = "";
  
# PIPE VALUES
$pipe = $_REQUEST['pipe'];

  
# demo classes
$demo = $_REQUEST['demo'];
$demo_output_format = "selected";
if ($demo == 1) {
  $demo_clusters = storeFile("demo_files/mcode.mde");
 }

title('convert-classes');

## Description
echo ("<center>Convert a classification / clustering result delivered in a given format by an algorithm into format that can be used by NeAT\n");

## Authors
echo ("<br>This program was developed by
	  <a target=_blank
          href='http://www.bigre.ulb.ac.be/people/Members/sylvain'>Sylvain
          Broh&eacute;e</a> and <a target=_blank
          href=http://www.bigre.ulb.ac.be/Users/jvanheld/>Jacques van
          Helden</a>.
          </center>\n");

## Open the form
echo ("<form method='post' action='convert_classes.php'
	 enctype='multipart/form-data'>");

################################################################
## Input options
echo("<hr>");
echo("<p><b>Input options</b></p>");

## Input format   
echo("<a href='help.convert_classes.html#formats'><B>Input format</B></a>");

if (!$pipe) {
  echo("&nbsp;<select name='in_format'>
         <option value = 'mcl'>MCL format
	 <option value = 'tab'>Neat tab delimited format
	 <option value = 'mcode' $demo_output_format>MCODE format
	 <option value = 'profiles'>profiles (table)
   </select><br>");
 } 


if ($pipe) {

} else {
  if ($demo) {
    demo("This demonstration clustering consist in the result of using Cytoscape on the <a href = 'http://www.ncbi.nlm.nih.gov/sites/entrez?db=pubmed&uid=10688190'>Uetz et al (2001)</a> yeast two-hybrid protein protein interaction network. It consists of 10 clusters");
  }

  ## Classes
  echo ("<b>Classes</b><br>");

  echo ("<textarea name='classes' rows='6' cols='65'>$demo_clusters</textarea>
    <br>Upload classes from file : <br>
    <input type='file' name='classes_file' size='45' /><br>
    &nbsp;&nbsp;&nbsp;
  
    <br><a href = 'help.convert_classes.html#columns'><B>Column specification (only relevant for tab-delimited input)</b></a><br>
    <table>
    <tr><td><B><a href = 'help.convert_classes.html#mcol'>Member</a></B></td><td><input type = 'text' name='m_col' value = '$default_mcol' size = 1></input></td></tr>
    <tr><td><B><a href = 'help.convert_classes.html#ccol'>Class</a></B></td><td><input type = 'text' name='c_col' value = '$default_ccol' size = 1></input></td></tr>
    <tr><td><B><a href = 'help.convert_classes.html#wcol'>Weight</a></B></td><td><input type = 'text' name='w_col' size = 1 value = '$default_wcol'></input></td></tr>
    </table>");
 }




################################################################
## Output options

echo("<hr>");
echo("<p><b>Output options</b></p>");

## Output format   
echo("<B><a href = 'help.convert_classes.html#formats'>Output format</a></B>&nbsp;<select name='out_format'>");
echo ("<option selected value = 'tab'>tab-delimited
<option value = 'profile'>Profile
</select><br>");




## Action buttons
echo ("<hr>");
echo ("<br>
  <table class='formbutton'>
  <TD><input type='submit' name='.submit' value='GO' /></TD>
  <TD><B><A HREF='convert_classes_form.php?demo=0'>RESET</A></B></TD>
  <TD><B><A HREF='convert_classes_form.php?demo=1'>DEMO</A></B></TD>
  </form>
  <TD><B><A HREF='help.convert_classes.html'>MANUAL</A></B></TD>
  <TD><B><A target = '_blank' HREF='".checkNeatTutorial("tutorials/neat_tutorial/Network_visualization_forma.html")."'>TUTORIAL</A></B></TD>
  <TD><B><A HREF='mailto:sylvain@bigre.ulb.ac.be'>MAIL</A></B></TD>
  </TR></TABLE>
 ");

?>
</body>
</html>
