<html>
<head>
<title>Network Analysis Tools - MCL</title>
<link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
   </head>
   <body class="form">
   <?php
   require ('functions.php');
// require ('demo_dataset.php');
# variable definition
$default_scol = 1;
$default_tcol = 2;
$default_wcol = "";
$default_inflation = 1.8;
$epsilon = 0.001;
  
# PIPE VALUES
$pipe = $_REQUEST['pipe'];
$graph_file = $_REQUEST['graph_file'];
$graph_format = $_REQUEST['graph_format'];
$scol = $_REQUEST['scol'];
$tcol = $_REQUEST['tcol'];
$wcol = $_REQUEST['wcol'];  
  
# demo graph
$demo = $_REQUEST['demo'];
if ($demo == 1) {
  $demo_graph = storeFile("demo_files/protein_interactions_gavin_2006_names.tab");
//   $demo_graph = storeFile("demo_files/rdm_test.tab");
;
  $demo_ecolors = "selected";
  $demo_ewidth = "checked";
  $default_ecolors = "";
//   $default_wcol = 3;
 }

title('MCL');
echo("<center>Fast and scalable unsupervised cluster algorithm for graphs based on simulation of (stochastic) flow in graphs.
	<br>The MCL program was developed by <a target=_blank href='http://micans.org/stijn/contact.html'>Stijn van Dongen</a> 
	(<a  target = _blank' href='http://igitur-archive.library.uu.nl/dissertations/1895620/inhoud.htm'>Van Dongen, 2000, PhD thesis</a>; 
	<a  target=_blank href='http://www.ncbi.nlm.nih.gov/sites/entrez?Db=pubmed&Cmd=ShowDetailView&TermToSearch=11917018&ordinalpos=6&itool=EntrezSystem2.PEntrez.Pubmed.Pubmed_ResultsPanel.Pubmed_RVDocSum'>Enright et al, 2002</a>)
	<br>The stand-alone version of MCL is available at <a target=_blank href='http://micans.org/mcl/'>http://micans.org/mcl/</a></center>");
  
echo ("<form method='post' action='mcl.php' enctype='multipart/form-data'>
  &nbsp;&nbsp;&nbsp;<a href = 'help.mcl.html#formats'><B>Input format</B></a>");
  
  
if (!$pipe) {
  echo("&nbsp;<select name='in_format'>
    <option selected value = 'tab'> tab-delimited format
    <option value = 'adj_matrix'> Adjacency matrix
    <option value = 'gml'> GML format
    </select><br>");
 } else {
  echo ": $graph_format<br>";
  echo "<input type='hidden' NAME='in_format' VALUE='$graph_format'>";
 } 
   
echo("<br>");
if (!$pipe) {
  if ($demo) {
    demo("This demonstration graph consists in the yeast co-immunopreciptation interaction dataset described in <a href = 'http://www.ncbi.nlm.nih.gov/sites/entrez?Db=pubmed&Cmd=ShowDetailView&TermToSearch=16429126&ordinalpos=1&itool=EntrezSystem2.PEntrez.Pubmed.Pubmed_ResultsPanel.Pubmed_RVDocSum'>Gavin et al (2006)</a>. It contains 1430 nodes and 6531 edges. The MCL algorithm is applied on it in order to highlight clusters of densely connected polypeptides.");
  }
  echo ("<b>Graph</b><br>");
  echo ("<textarea name='graph' rows='6' cols='65'>$demo_graph</textarea>
    <br>Upload graph from file : <br>
    <input type='file' name='graph_file' size='45' /><br>
    &nbsp;&nbsp;&nbsp;
  
    <br><a href = 'help.mcl.html#columns'>Column specification (only relevant for tab-delimited input)</a><br>
    <table>
    <tr><td><B><a href = 'help.mcl.html#scol'>Source node</a></B></td><td><input type = 'text' name='s_col' value = '$default_scol' size = 1></input></td></tr>
    <tr><td><B><a href = 'help.mcl.html#scol'>Target node</a></B></td><td><input type = 'text' name='t_col' value = '$default_tcol' size = 1></input></td></tr>
    <tr><td><B><a href = 'help.mcl.html#wcol'>Weight or label column</a></B></td><td><input type = 'text' name='w_col' size = 1 value = '$default_wcol'></input></td></tr>
    </table>");
 } else {
  info_link("Graph uploaded from the previous treatment", rsat_path_to_url($graph_file));
  echo "<input type='hidden' NAME='pipe_graph_file' VALUE='$graph_file'>";
 }

if ($graph_format == 'tab') {
  echo "<input type='hidden' NAME='s_col' VALUE='$scol'/>\n";
  echo "<input type='hidden' NAME='t_col' VALUE='$tcol'/>\n";
  echo "<input type='hidden' NAME='w_col' VALUE='$wcol'/>\n";
 }


echo("<B><a href = 'help.mcl.html#inflation'>Inflation value</a></B>&nbsp;<select name='inflation'>");
for ($i = 1.1; $i <= 10; $i = $i+0.1) {
  $select = "";
  if (abs($i - $default_inflation) < $epsilon) {
//     echo (abs($i - $default_inflation))."\n";
    $select = "selected";
  }
  echo ("<option $select value='$i'>$i\n"); 
}
echo ("</select><br>");

// echo("
//   <B><a href = 'help.mcl.html#inflation'>Inflation value</a></B></td><td><input type = 'text' name='inflation' value = '$default_inflation' size = 3></input> (a value between 1.0 and 5.0)</td></tr>");

echo ("<ul><ul><table class='formbutton'>
	 <TD><input type='submit' name='.submit' value='GO' /></TD>
	 <TD><B><A HREF='mcl_form.php?demo=0'>RESET</A></B></TD>
	 <TD><B><A HREF='mcl_form.php?demo=1'>DEMO</A></B></TD>
   </form>
  <TD><B><A HREF='help.mcl.html'>MANUAL</A></B></TD>
  <TD><B><A target = '_blank' HREF='".checkNeatTutorial("tutorials/neat_tutorial/Graph_clustering.html")."'>TUTORIAL</A></B></TD>

  <TD><B><A HREF='mailto:sylvain@bigre.ulb.ac.be'>MAIL</A></B></TD>
   </TR></TABLE></ul></ul>
 ");

?>
</body>
</html>
