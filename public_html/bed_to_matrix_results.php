<?php
    session_start();
?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
  <head>
    <title>RSAT <?php print $_SESSION['previous_bed_to_matrix']; ?></title>
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
    <link rel="stylesheet" type="text/css" href="css/main.css" media="screen,projection,print"/>
	 <link rel="stylesheet" type="text/css" href="css/tabs.css" media="screen,projection,print"/> 
    <link rel="stylesheet" type="text/css" href = "css/main_grat.css" media="screen">
    </head>
    
    <body class="results"> 
    
      <h3 align='center'><a href="http://pedagogix-tagc.univ-mrs.fr/rsat/">RSAT</a> - Features matrix - Results</h3>
      <br/> 
      <fieldset>  
         <legend><b>Result files</b></legend>    
        <?php
        $workingdir = $_SESSION['workingdir_bed_to_matrix'];
        

        $fastas = glob($workingdir.'/*.{fasta}', GLOB_BRACE);
        print "<fieldset><legend><b>Sequence files (fasta)</b></legend> ";
        foreach ($fastas as $fasta) {
            if ($fasta != "." && $fasta != "..") {
                print "<br><a href='$fasta'>".basename($fasta)."</a>";
                $i++;
            }
        }
        print "</fieldset>";
        
        $beds = glob($workingdir.'/*_userbed.{bed}', GLOB_BRACE);
        print "<fieldset><legend><b>Coordinates files (bed)</b></legend> ";
        foreach ($beds as $bed) {
            if ($bed != "." && $bed != "..") {
                print "<br><a href='$bed'>".basename($bed)."</a>";
                $j++;
            }
        }
        print "</fieldset>";
        
        $fts = glob($workingdir.'/*.{ft}', GLOB_BRACE);
        print "<fieldset><legend><b>Matrix-scan output files (ft)</b></legend> ";
        foreach ($fts as $ft) {
            if ($ft != "." && $ft != "..") {
                print "<br><a href='$ft'>".basename($ft)."</a>";
                $k++;
            }
        }
        print "</fieldset>";
        
        $tabs = glob($workingdir.'/feature_matrix*.{tab}', GLOB_BRACE);
        print "<fieldset><legend><b>Output matrices (tab)</b></legend> ";
        foreach ($tabs as $tab) {
            if ($tab != "." && $tab != "..") {
                print "<br><a href='$tab'>".basename($tab)."</a>";
                $l++;
            }
        }
        print "</fieldset>";
        $expected_workingdir_file_nb = $i+$j+$k+$l;
        if ($expected_workingdir_file_nb < 7) {
            print "Refreshing...";
            $page = $_SERVER['PHP_SELF'];
            $sec = "5";
            header("Refresh: $sec; url=$page");
        }
        	$outmatrixfile = $workingdir."/feature_matrix.tab";
        	
        
			        
        
       echo "

         </fieldset>  
        
        <table class = 'Nextstep'  align=center>
    <tr></tr>
    <tr>
    <td align=center>
	    <form method='post' action='bed_to_matrix_form.php'>
	    <input type='submit' value='Back to the form'>
	    </form>
      </td>
      <td align=center>
	    <form method='post' action='tune_svm_form.php'>
	    <input type='hidden' name='matrixfile' id='matrixfile' value='".$outmatrixfile."'>
	    <input type='submit' value='SVM tuning'>
	    </form>
      </td>
    </tr>

  </table>
</body>
</html>
"

?>
