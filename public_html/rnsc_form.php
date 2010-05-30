
<html>
<head>
<title>Network Analysis Tools - RNSC</title>
<link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
   </head>
   <body class="form">
   <?php
   require ('functions.php');
// require ('demo_dataset.php');
# variable definition
$default_scol = 1;
$default_tcol = 2;
$default_div_freq = 50;
$default_div_leng = 50;
$default_tabu_leng = 50;
$default_tabu_list_tol = 1;
$default_nb_exp = 3;
$default_shf_div_freq = 3;
$default_naive_stop = 15;
$default_scaled_stop = 15;
$default_maxnb_cluster = "";
$div_freq_values = array(10,20,50);
$tabu_leng_values = array(1,10,50,100);
$tabu_list_tolerance_values = array(1,3,5);
$scaled_stopping_tolerance_values = array(1,5,15);
$naive_stopping_tolerance_values = array(1,5,15);
$shf_diversification_len_values = array(1,3,9);
$number_of_experiment_values = array(1,3,5,10);
$diversification_frequency_values = array(10,20,50);



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
  $demo_ecolors = "selected";
  $demo_ewidth = "checked";
  $default_ecolors = "";
//   $default_wcol = 3;
 }

title('RNSC');
echo("<center> RNSC - Restricted Neighbourhood Search Cluster Algorithm
	<br>The RNSC program was developed by <a target=_blank href='http://cgm.cs.mcgill.ca/~aking6/'>Andrew King</a>. <br>RNSC is an efficient cost-based local search clustering algorithm that explores the solution space to minimize a cost function, calculated according to the numbers of intracluster and inter-cluster edges. 
	(<a  target = _blank' href='http://cgm.cs.mcgill.ca/~aking6/papers/thesis.pdf'>King, 2004, M.Sc. thesis</a>; 
	<a  target=_blank href='http://cgm.cs.mcgill.ca/~aking6/papers/bio_adv.pdf'>King et al, 2004</a>)
	<br>The stand-alone version of RNSC is available upon request.</center>");
  
echo ("<form method='post' action='rnsc.php' enctype='multipart/form-data'>
  &nbsp;&nbsp;&nbsp;<a href = 'help.rnsc.html#formats'><B>Input format</B></a>");
  
  
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
    demo("This demonstration graph consists in the yeast co-immunopreciptation interaction dataset described in <a href = 'http://www.ncbi.nlm.nih.gov/sites/entrez?Db=pubmed&Cmd=ShowDetailView&TermToSearch=16429126&ordinalpos=1&itool=EntrezSystem2.PEntrez.Pubmed.Pubmed_ResultsPanel.Pubmed_RVDocSum'>Gavin et al (2006)</a>. It contains 1430 nodes and 6531 edges. The rnsc algorithm is applied on it in order to highlight clusters of densely connected polypeptides.");
  }
  echo ("<b>Graph</b><br>");
  echo ("<textarea name='graph' rows='6' cols='65'>$demo_graph</textarea>
    <br>Upload graph from file : <br>
    <input type='file' name='graph_file' size='45' /><br>
    &nbsp;&nbsp;&nbsp;
  
    <br><a href = 'help.rnsc.html#columns'>Column specification (only relevant for tab-delimited input)</a><br>
    <table>
    <tr><td><B><a href = 'help.rnsc.html#scol'>Source node</a></B></td><td><input type = 'text' name='s_col' value = '$default_scol' size = 1></input></td></tr>
    <tr><td><B><a href = 'help.rnsc.html#scol'>Target node</a></B></td><td><input type = 'text' name='t_col' value = '$default_tcol' size = 1></input></td></tr>
    </table>");
 } else {
  info_link("Graph uploaded from the previous treatment", rsat_path_to_url($graph_file));
  echo "<input type='hidden' NAME='pipe_graph_file' VALUE='$graph_file'>";
 }

if ($graph_format == 'tab') {
  echo "<input type='hidden' NAME='s_col' VALUE='$scol'/>\n";
  echo "<input type='hidden' NAME='t_col' VALUE='$tcol'/>\n";
 }

# Parameter list
echo("
<table>");
# max number of cluster
echo ("<tr><td><B><a href = 'help.rnsc.html#max_clust'>Maximum number of cluster</a></B></td><td><input type = 'text' name='max_clust' value = '$default_maxnb_cluster' size = 1></input></td></tr>");
# tabu length
echo ("<tr><td><B><a href = 'help.rnsc.html#tabulength'>Tabu length</a></B></td><td><select name='tabulength'>");
for ($i = 0; $i < count($tabu_leng_values); $i++) {
  $select = "";
  if ($tabu_leng_values[$i] == $default_tabu_leng) {
    $select = "selected";
  }
  echo ("<option $select value='$tabu_leng_values[$i]'>$tabu_leng_values[$i]\n"); 
}
echo ("</select>");
echo ("</td></tr>");
# tabu list tolerance
echo ("<tr><td><B><a href = 'help.rnsc.html#tabulist'>Tabu list tolerance</a></B></td><td><select name='tabulist'>");
for ($i = 0; $i < count($tabu_list_tolerance_values); $i++) {
  $select = "";
  if ($tabu_list_tolerance_values[$i] == $default_tabu_list_tol) {
    $select = "selected";
  }
  echo ("<option $select value='$tabu_list_tolerance_values[$i]'>$tabu_list_tolerance_values[$i]\n"); 
}
echo ("</select>");
echo ("</td></tr>");
# naive stopping tolerance
echo ("<tr><td><B><a href = 'help.rnsc.html#naive_stop'>Naive stopping tolerance</a></B></td><td><select name='naive_stop'>");
for ($i = 0; $i < count($naive_stopping_tolerance_values); $i++) {
  $select = "";
  if ($naive_stopping_tolerance_values[$i] == $default_naive_stop) {
    $select = "selected";
  }
  echo ("<option $select value='$naive_stopping_tolerance_values[$i]'>$naive_stopping_tolerance_values[$i]\n"); 
}
echo ("</select>");
echo ("</td></tr>");
# scaled stopping tolerance
echo ("<tr><td><B><a href = 'help.rnsc.html#scale_stop'>Scaled stopping tolerance</a></B></td><td><select name='scale_stop'>");
for ($i = 0; $i < count($scaled_stopping_tolerance_values); $i++) {
  $select = "";
  if ($scaled_stopping_tolerance_values[$i] == $default_naive_stop) {
    $select = "selected";
  }
  echo ("<option $select value='$scaled_stopping_tolerance_values[$i]'>$scaled_stopping_tolerance_values[$i]\n"); 
}
echo ("</select>");
echo ("</td></tr>");
# Diversification frequency
echo ("<tr><td><B><a href = 'help.rnsc.html#div_freq'>Diversification frequency</a></B></td><td><select name='div_freq'>");
for ($i = 0; $i < count($diversification_frequency_values); $i++) {
  $select = "";
  if ($diversification_frequency_values[$i] == $default_div_freq) {
    $select = "selected";
  }
  echo ("<option $select value='$diversification_frequency_values[$i]'>$diversification_frequency_values[$i]\n"); 
}
echo ("</select>");
echo ("</td></tr>");

#  shuffling diversification length
echo ("<tr><td><B><a href = 'help.rnsc.html#shf_div_len'>Shuffling diversification length</a></B></td><td><select name='shf_div_len'>");
for ($i = 0; $i < count($shf_diversification_len_values); $i++) {
  $select = "";
  if ($shf_diversification_len_values[$i] == $default_shf_div_freq) {
    $select = "selected";
  }
  echo ("<option $select value='$shf_diversification_len_values[$i]'>$shf_diversification_len_values[$i]\n"); 
}
echo ("</select>");
echo ("</td></tr>");

#  number of experiment
echo ("<tr><td><B><a href = 'help.rnsc.html#exp_nb'>Number of experiments</a></B></td><td><select name='exp_nb'>");
for ($i = 0; $i < count($number_of_experiment_values); $i++) {
  $select = "";
  if ($number_of_experiment_values[$i] == $default_nb_exp) {
    $select = "selected";
  }
  echo ("<option $select value='$number_of_experiment_values[$i]'>$number_of_experiment_values[$i]\n"); 
}
echo ("</select>");
echo ("</td></tr>");


echo ("</table>");






// echo("
//   <B><a href = 'help.rnsc.html#inflation'>Inflation value</a></B></td><td><input type = 'text' name='inflation' value = '$default_inflation' size = 3></input> (a value between 1.0 and 5.0)</td></tr>");

echo ("<ul><ul><table class='formbutton'>
	 <TD><input type='submit' name='.submit' value='GO' /></TD>
	 <TD><B><A HREF='rnsc_form.php?demo=0'>RESET</A></B></TD>
	 <TD><B><A HREF='rnsc_form.php?demo=1'>DEMO</A></B></TD>
   </form>
  <TD><B><A HREF='help.rnsc.html'>MANUAL</A></B></TD>
  <TD><B><A target = '_blank' HREF='".checkNeatTutorial("tutorials/neat_tutorial/Graph_clustering.html")."'>TUTORIAL</A></B></TD>
  <TD><B><A HREF='mailto:sylvain@bigre.ulb.ac.be'>MAIL</A></B></TD>
   </TR></TABLE></ul></ul>
 ");

?>
</body>
</html>
