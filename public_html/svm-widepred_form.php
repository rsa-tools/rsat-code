<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
  <head>
    <title>RSAT - SVM model prediction.</title>
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
    <link rel="stylesheet" type="text/css" href="css/main.css" media="screen,projection,print"/>
    <link rel="stylesheet" type="text/css" href="css/tabs.css" media="screen,projection,print"/>
    <script src="js/RSAT_menu.js" type="text/javascript"></script>
    <script src="js/RSAT_tabs.js" type="text/javascript"></script>
    <!-- <script type="text/javascript">
			function add_demo() {
				document.getElementById('genome').value = 'dm3';
				document.forms[0].bed.value="";
				document.getElementById('bedfile').value='/data/rsat/public_html/data/svm_rsat/data/CRM_13-16.bed'; 				
				document.forms[0].header[0].checked=true;			
				document.forms[0].output[0].checked=true;
				document.getElementById('xtime').value="5";
				document.getElementById('expName').value='HeartCRM_13-16h';									
 				document.forms[0].sequence_url.value = '/data/rsat/public_html/data/svm_rsat/data/CRM_13-16.bed';
}
		</script> -->
  </head>

  <body class="form">
    <div>
      <h3 align='center'><a href="http://pedagogix-tagc.univ-mrs.fr/rsat/">RSAT</a> - SVM model prediction.</h3>
      <br/>
      <form method='post' action='svm-widepred.php' enctype='multipart/form-data'>
 		<fieldset>
 		<legend><b>SVM model <font style='color:orange'>(mandatory)</font></b></legend> 
				Upload SVM model (.rda file), in any of the three following ways: 
			       <ul type='square'>
	    <li>Paste feature file<br/>
	      
              <textarea name='modelfile' rows='6' cols='45'></textarea></li>
	   <li>Specify the URL of ft file on remote server (e.g. Galaxy)<br/>
              <input type="text" name="sequence_url" size="62" /><br/></li>
            <li>Upload a file from your computer<br/>
	      <input type='file' name='modelfile' id='modelfile' size='40' /></li>
	  </ul>
	  </p>
	   <br/><br/>
	    </fieldset> 
	   <fieldset>
	   <legend><b>Experiment Name <font style='color:orange'>(mandatory)</font></b></legend> 
	   <b>Give experience name &nbsp;&nbsp;</b><input type="text" id="expName" name="expName" size="10" placeholder="ExpName" /><br/>      
          </fieldset>   
          <fieldset>
	   <legend><b>Genomic regions to scan <font style='color:orange'>(mandatory)</font></b></legend> 
	   Upload <a href='https://pedagogix-tagc.univ-mrs.fr/rsat/features-matrix_form.php'>csv file</a> of regions of interest (warning: you have to use strictly the same features as for building SVM model), in any of the three following ways: 
			       <ul type='square'>
	    <li>Paste feature file<br/>
	      
              <textarea name='modelfile' rows='6' cols='45'></textarea></li>
	   <li>Specify the URL of ft file on remote server (e.g. Galaxy)<br/>
              <input type="text" name="sequence_url" size="62" /><br/></li>
            <li>Upload a file from your computer<br/>
	      <input type='file' name='modelfile' id='modelfile' size='40' /></li>
	  </ul>
	  </p>
	   <br/><br/>   
	    </fieldset>     
    <p>
				
				
 		</fieldset>
      <ul><table class='formbutton'>
        <tr valign=middle>
          <td><input type="submit" name="submit" value="GO" /></td>
          <td><input type="reset"  name="reset" /></td> 
          <td><input type="button" name="demo" value="Demo" onclick="add_demo()"/></td>  
          <td><b><a href='help.fetch-sequences.html'>[MANUAL]</a></b></td>
          <td><b><a href='http://www.bigre.ulb.ac.be/forums/' target='_top'>[ASK A QUESTION]</a></b></td>
        </tr></table>
      </ul>
    </div>
  </body>
</html>

