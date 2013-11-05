<html>
<head>
   <title>Network Analysis Tools - download string dataset</title>
   <link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
   <link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">

</head>
<body class="form">
<?php
  require ('functions.php');
  require ('string_functions.php');
  $genes = $_REQUEST['genes'];
  $genes = trim_text($genes);
  $gene_list = explode("\n", $genes); 
  $organism_id = $_REQUEST['organism'];
  $gene_description = resolveName($gene_list, $organism_id);
  $channels = array();
  $gene_list = array();
  $uth = $_REQUEST['uth'];
  $lth = $_REQUEST['lth'];
  
  
  
  if ($_REQUEST['automated_textmining'] != "") {
    array_push ($channels, 'automated_textmining');
  }
  if ($_REQUEST['experimental_interaction_data'] != "") {
    array_push ($channels, 'experimental_interaction_data');
  }
  if ($_REQUEST['gene_cooccurence'] != "") {
    array_push ($channels, 'gene_cooccurence');
  }
  if ($_REQUEST['gene_fusion_events'] != "") {
    array_push ($channels, 'gene_fusion_events');
  }
  if ($_REQUEST['gene_coexpression'] != "") {
    array_push ($channels, 'gene_coexpression');
  }
  if ($_REQUEST['combined_confidence'] != "") {
    array_push ($channels, 'combined_confidence');
  }  
  $channels = implode(",", $channels);
  
  title('STRING : genes selection');
  
   echo ("<center>Confirm the selected genes</center>\n");
   echo ("<form method='post' action='string_dataset.php' enctype='multipart/form-data'> ");
   if (count($gene_description) > 0) { 
    echo("<table border = 1>");
    echo("<tr><td></td><td>Usual name</td><td>String ID</td><td>Description</td></tr>");
    for ($i = 0; $i < count($gene_description); $i++) {
      $varname = $gene_description[$i][1];
      $varname = trim ($varname);
      $varname = str_replace(".", "_", $varname);
      $checked = "checked";
      if ($gene_description[$i][1] == "Not found in STRING") {
        $checked = "";
        $not_found = 1;
      }
      echo ("<tr>");
      if ($checked != "") {
        echo ("<td><input type='checkbox' name='".$varname."' value='1'/></td>");
      } else {
        echo ("<td></td>");
      }
      echo ("<td><b>".$gene_description[$i][0]."</b></td><td>".$gene_description[$i][1]."</td><td>".$gene_description[$i][2]."</td>");
      array_push ($gene_list, $varname);
    }
    echo "</table>";
  }
  $gene_list = implode ($gene_list, "-sep-");
  echo "<input type='hidden' NAME='channels' VALUE='$channels'>";
  echo "<input type='hidden' NAME='organism' VALUE='$organism_id'>";
  echo "<input type='hidden' NAME='genes' VALUE='$gene_list'>";
  echo "<input type='hidden' NAME='uth' VALUE='$uth'>";
  echo "<input type='hidden' NAME='lth' VALUE='$lth'>";
  if ((count ($gene_description) == 1 && $not_found == 1)) { 
    # IF UNKNOWN GENE -> ERROR
    error ("Your query did not match any record in STRING");
  } else if (count($gene_description) == 0) { 
    # IF NO GENE PROVIDED -> ERROR
    error ("You must at least provide one gene name");
  } else {
    echo("  
      <ul><ul><table class='formbutton'>
      <TD><input type='submit' name='.submit' value='GO' /></TD>
      </form>

      <TD><B><A HREF='mailto:sylvain-at-bigre.ulb.ac.be'>MAIL</A></B></TD>
      </TR></TABLE></ul></ul>
   ");
  }
?>
</body>
</html>
