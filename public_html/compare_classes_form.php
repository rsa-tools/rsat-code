<html>
<head>
<title>NeAT - compare_classes</title>
<link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
   </head>
   <body class="form">


   <?php
   require ('functions.php');
// require ('demo_dataset.php');
$default_min_sig = 0;
  
# PIPE VALUES
$pipe = $_REQUEST['pipe'];
$pipe_Q_class_file = $_REQUEST['class_file'];
  
# demo graph
$demo = $_REQUEST['demo'];

if ($demo == 1) {
  $demo_classesQ = storeFile("demo_files/gavin_mcl_clusters_inf2.1.tab") ;
  $demo_classesR = storeFile("demo_files/mips_complexes.tab") ;
      
  $demo_remark = "This demonstration consists in the comparaison between clusters obtained after application of the <a href = 'http://micans.org/mcl/' target = 'top'>MCL</a> clustering algorithm to the <a target = '_blank' href = 'http://www.ncbi.nlm.nih.gov/sites/entrez?Db=pubmed&Cmd=ShowDetailView&TermToSearch=16429126&ordinalpos=1&itool=EntrezSystem2.PEntrez.Pubmed.Pubmed_ResultsPanel.Pubmed_RVDocSum'>Gavin et al (2006)</a> interaction network and the complexes annotated in the <a target = '_blank' href = 'http://mips.gsf.de/'>MIPS database</a>.";
}
title('compare-classes');
?>

<center>Compare two classifications (clustering results, functional
  classes, ...), and assess the statistical significance of common
  members between each pair of classes.

<br>This program was
  developed <A HREF='mailto:jacques.van.helden@ulb.ac.be (Jacques van
  Helden)'>Jacques van Helden</a>, with a contribution of Joseph Tran
  for a first prototype.</center>

<hr>
  
<form method='post' action='compare_classes.php' enctype='multipart/form-data'>

<?php
if($demo) {
  demo($demo_remark);
 }
?>

<blockquote>
<table border=0 cellspacing=0 cellpadding=4>

<tr align = 'center'><th>Query classes</th><th>Reference classes</th></tr>

<?php
################################################################
## Query input classes textarea
echo "<tr><td>\n";
if (!$pipe) {
  echo ("<textarea name='classesQ' rows='6' cols='40'>$demo_classesQ</textarea>");
 } 
echo "</td>";

################################################################
## Reference input classes textarea
echo ("<td><textarea name='classesR' rows='6' cols='40'>$demo_classesR</textarea></td></tr>");
echo ("<tr align = 'center'><td>");

################################################################
## query input classes file
if (!$pipe) {
  echo ("Upload query classes from file : <br>
    <input type='file' name='Qclass_file' size='40' />");
 } else {
    info_link("Query classes loaded from the previous treatment", rsat_path_to_url($pipe_Q_class_file));
  echo "<input type='hidden' NAME='pipe_Q_class_file' VALUE='$pipe_Q_class_file'/>";
 }
echo "</td>";

################################################################
## reference input classes file
echo ("<td>Upload reference classes from file : <br>
<input type='file' name='Rclass_file' size='40' /></td>");
echo ("</tr>");
echo ("</table>");
echo ("</blockquote>");

################################################################
## Score column (for the dot product)
if (!$pipe) {

  echo ("<b><a href='help.compare_classes.html#score_col'>Score column</a></b>
        <input type = 'text' name='score_col' size = 2\> 
        (optional; if specified, must be valid for both query and reference classes)");
}
?>

<!-- Comparaison between query classes - options for self-comparison-->
<br><input type = 'radio' checked value='off' name='self_compa' size = 1> 
      <a href = 'help.compare_classes.html#query_vs_ref'>
      <b>Compare query classes to reference classes</a></b>

<br><input type = 'radio' value='on' name='self_compa' size = 1> 
      <a href = 'help.compare_classes.html#query_vs_query'>
      <b>Compare query classes to query classes</a></b> (do not specify reference classes)<br>

<ul><input type='checkbox' value='on' name='distinct' size = 1 checked/> 
      <a href = 'help.compare_classes.html#distinct'>
      <b>Prevent self-comparison</a></b><br>

<input type = 'checkbox' value='on' name='triangle' size = 1 checked/>
      <a href = 'help.compare_classes.html#triangle'>
      <b>Prevent reciprocal comparisons</a></b><br></ul></td>

<br><table>
   <tr><td colspan = 2 >&nbsp;&nbsp;&nbsp;<B><a href = 'help.compare_classes.html#formats'>Output format</a></B>
   <tr><td><input type = 'radio' checked name='out_format' value = 'class'/></td><td>Class pairs</td></tr>
   <tr><td><input type = 'radio' name='out_format' value = 'matrix'/></td><td>Reference/query matrix</td></tr>

<?php

echo "<hr>";
echo "<b>Output options</b><br>";
echo "<blockquote>";
echo ("<table border=1  cellspacing=0 cellpadding=4>");
echo ("<tr align=center><th><input type = 'radio' checked name='out_format' value = 'class'/>Class pairs</th>");
echo ("<th><input type = 'radio' name='out_format' value = 'matrix'/>Reference/query matrix</th></tr>");
//echo ("<td><b>Parameters for output type 'Classes file'</b></td><td><b>Parameters for output type 'Matrix file'</b></td></tr><tr align='center' valign='TOP'><br>");

# CLASS FILE OUTPUT PARAMETERS
echo("<td><table border='0' cellspacing='0' cellpadding='0'>
  <tr>  <th> <A HREF='help.compare_classes.html#return_fields'>Return fields</A> </th></tr> 
  <tr><td><label><input type='checkbox' name='occ' value='on' checked='checked' />Occurrences</label></td></tr> 
  <tr><td><label><input type='checkbox' name='freq' value='on' />Frequencies</label></td></tr> 
  <tr><td><label><input type='checkbox' name='proba' value='on' checked='checked' />Hypergeometric probability</label>
  &nbsp;&nbsp;<label><input type='checkbox' name='sort' value='on' checked='checked' />Sorting criterion</label></td></tr> 
  <tr><td><label><input type='checkbox' name='jac_sim' value='on' checked='checked' />Jaccard similarity</label></td></tr> 
  <tr><td><label><input type='checkbox' name='sor_sim' value='on' />Sorensen similarity</label></td></tr> 
  <tr><td><label><input type='checkbox' name='dotprod' value='on' />Dotproduct <i>(Only relevant if a score column is specified)</i></label></td></tr>
  <tr><td><label><input type='checkbox' name='entropy' value='on' />Entropy</label></td></tr> 
  <tr><td><label><input type='checkbox' name='members' value='on' /> Members <i> (Beware, this might generate large result files)</i></label></td></tr>
");

echo("</table><br><br>\n");     
echo("<table border='0' cellspacing='0' cellpadding='0'><tr>
     <th><A HREF='help.compare_classes.html#thresholds'>Thresholds</A></th>
     <th>Lower</th>
     <th>Upper</th></tr>\n"); 
echo("<tr align='left' valign='TOP'>
     <td> Query size </td> 
     <td><input type='text' name='lth_q' value='1' size='5' /></td> 
     <td><input type='text' name='uth_q'  size='5'  value='none'/></td></tr>\n"); 
echo("<tr align='left' valign='TOP'>
     <td> Reference size </td> 
     <td><input type='text' name='lth_r' value='1'  size='5' /></td> 
     <td><input type='text' name='uth_r'  size='5' value='none' /></td></tr>\n"); 
echo("<tr align='left' valign='TOP'>
     <td> Intersection size </td> 
     <td><input type='text' name='lth_qr' value='1'  size='5' /></td> 
     <td><input type='text' name='uth_qr'  size='5' value='none' /></td></tr>\n"); 
echo("<tr align='left' valign='TOP'>
     <td> Significance </td> 
     <td><input type='text' name='lth_sig' value='$default_min_sig' size='5' /></td> 
     <td><input type='text' name='uth_sig' value='none' size='5' value='none' /></td></tr>\n"); 
echo("<tr align='left' valign='TOP'>
     <td> P-value </td> 
     <td><input type='text' name='lth_pval' value='none' size='5' /></td> 
     <td><input type='text' name='uth_pval'  size='5' value='none' /></td></tr>\n"); 
echo("<tr align='left' valign='TOP'>
     <td> E-value </td> 
     <td><input type='text' name='lth_eval' value='none' size='5' /></td> 
     <td><input type='text' name='uth_eval'  size='5' value='none' /></td></tr>\n"); 
echo("<tr align='left' valign='TOP'>
     <td> Jaccard similarity </td> 
     <td><input type='text' name='lth_jac' value='none' size='5' /></td> 
     <td><input type='text' name='uth_jac'  size='5' value='none' /></td></tr>\n"); 
echo("<tr align='left' valign='TOP'>
     <td> Sorensen similarity </td> 
     <td><input type='text' name='lth_sor' value='none' size='5' /></td> 
     <td><input type='text' name='uth_sor'  size='5' value='none' /></td></tr>\n"); 
echo("<tr align='left' valign='TOP'>
     <td> Mutual information </td> 
     <td><input type='text' name='lth_mi' value='none' size='5' /></td> 
     <td><input type='text' name='uth_mi'  size='5' value='none' /></td></tr>\n"); 
echo("<tr align='left' valign='TOP'>
     <td> Dot product <i><br>(Only relevant if a score column is specified)</i></td> 
     <td><input type='text' value='none' name='lth_dp'  size='5' /></td> 
     <td><input type='text' name='uth_dp' value='none' size='5' /></td></tr>\n");
echo("</table><br>\n");
echo ("</td>\n");


# MATRIX FILE OUTPUT PARAMETERS
echo ("<td><table>");
echo ("<tr><td><input type='radio' name='matrix_field' checked='checked' value='QR'/></td><td>Intersection</td></tr>\n
         <tr><td><input type='radio' name='matrix_field' value='sig'/></td><td>Significance</td></tr>\n
         <tr><td><input type='radio' name='matrix_field' value='jac_sim'/></td><td>Jaccard similarity</td></tr>\n
         <tr><td><input type='radio' name='matrix_field' value='sor_sim'/></td><td>Sorensen similarity</td></tr>\n
         <tr><td><input type='radio' name='matrix_field' value='dotprod'/></td><td>Dot product</td></tr>\n
         <tr><td><input type='radio' name='matrix_field' value='E_val'/></td><td>E-value</td></tr>\n
         <tr><td><input type='radio' name='matrix_field' value='P_val'/></td><td>P-value</td></tr>\n
         <tr><td><input type='radio' name='matrix_field' value='I(Q,R)'/></td><td>Mutual information</td></tr>\n
         </table>");
echo ("</td></tr>\n</table>");       
echo "</blockquote>";

## Buttons
echo ("<ul><ul><table class='formbutton'>
  <TD><input type='submit' name='.submit' value='GO' /></TD>
  <TD><B><A HREF='compare_classes_form.php?demo=0'>RESET</A></B></TD>
  <TD><B><A HREF='compare_classes_form.php?demo=1'>DEMO</A></B></TD>
  </form>
  <TD><B><A HREF='help.compare_classes.html'>MANUAL</A></B></TD>
  <TD><B><A HREF='mailto:sylvain@bigre.ulb.ac.be'>MAIL</A></B></TD>
  </TR></TABLE></ul></ul>");
echo "<hr>";
?>
  


</body>
</html>
