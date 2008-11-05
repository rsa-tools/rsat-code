<html>
<head>
   <title>NAT - KEGG network provider</title>
   <link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
</head>
<body class="form">
<?
  require ('functions.php');

  // default variables
  $default_organisms = "";
  $default_organisms_file = "";
  $default_reactions = "";
  $default_reactions_file = "";
  $default_excludedcompounds = "";
  $default_excludedreactions = "";
  $default_excludedrpairclasses = "";
  $default_out_format = "";
  $default_email='';

  # demo
  $demo = $_REQUEST['demo'];
  if ($demo == 1) {
     $demo_organisms = "dsba/dsmi/dspd/spo/sce";
     $default_organisms = $demo_organisms;
  }

  title('KEGG network provider');
  echo ("<center>Obtain organism-specific KEGG networks from KEGG PATHWAY version 46.0. Click on <a href='help.keggnetworkprovider.html'><img src='images/question-mark-button.jpg' width='15' alt='help'></a> for help.<br>
  Web service and interface by Karoline Faust.</center>\n");
  echo ("<form method='post' action='keggnetworkprovider.php' enctype='multipart/form-data'>
  <hr>
  <h2>1. Input <a href='help.keggnetworkprovider.html#input'><img src='images/question-mark-button.jpg' width='15' alt='help'></a>
  </h2>
  <br>
  <br>");
  if($demo == 1){
  	demo("To demonstrate the KEGG network provider, we are constructing the metabolic network from five yeast species:<br>
  	dsba (Saccharomyces bayanus), dsmi (Saccharomyces mikatae), dspd (Saccharomyces paradoxus), spo (Schizosaccharomyces pombe)
  	 and sce (Saccharomyces cerevisiae).
  	  <input type='hidden' name='organisms' value='$default_organisms'></input>
  	");
  }else{
  echo("
  <table>
  <tr><td>
         Enter KEGG abbreviations for organisms, separated by '/' in text field below (e.g. eco/sce/reference):
  </td></tr>
  <tr><td>
        <input type='text' name='organisms' value='$default_organisms' size=100></input>
   </td></tr>
   <tr><td>Check<a href='data/KEGG/Kegg_organisms_list.txt'> here</a> for organism abbreviations.</td></tr>
   <tr><td>
        AND/OR<br>
        Upload a file which contains one organism KEGG identifier by line:
   </td></tr>
   <tr><td>
       <input type='file' name='organisms_file' size='45' />
   </td></tr>
   <tr><td>
   <br>
   <b>AND/OR</b>
   <br>
   </td></tr>
   <tr><td>
        Enter KEGG reaction identifiers, separated by '/' in text field below (e.g. R00015/R00551):
   </td></tr>
   <tr><td>
     <input type='text' name='reactions' value='$default_reactions' size=100>
   </td></tr>
    <tr><td>
    AND/OR<br>
    Upload a file which contains one reaction KEGG identifier by line:
    </td></tr>
    <tr><td>
          <input type='file' name='reactions_file' size='45' />
    </td></tr>
  </table>
  <br>");
  }
  echo("
   <br>
   <hr>
   <h2>2. Options
   <a href='help.keggnetworkprovider.html#options'><img src='images/question-mark-button.jpg' width='15' alt='help'></a>
   </h2>
   <br>
   <br>
   <br>
  <table>
     <tr><td><b>Optional:</b> KEGG identifiers of compounds to exclude (separated by '/')</td><td><input type='text' name='excludedcompounds' value='$default_excludedcompounds' size=80></input></td></tr>
     <tr><td><b>Optional:</b> KEGG identifiers of reactions to exclude (separated by '/')</td><td><input type= 'text' name='excludedreactions' value='$default_excludedreactions' size=80></input></td></tr>
     <tr><td><b>Optional:</b> RPAIR classes to exclude (separated by '/')</td><td><input type= 'text' name='excludedrpairclasses' value='$default_excludedrpairclasses' size=80></input></td></tr>
  </table>
  <br>
  <br>
    <table>
         <tr><td><B>directed network</B><input type='checkbox' name='directed' value='on'></input></td></tr>
         <tr><td><B>RPAIR network</B>&nbsp;&nbsp;&nbsp;<input type='checkbox' name='rpair' value='on'></input></td></tr>
         <tr><td><br></td></tr>
         <tr><td><b>Optional:</b> Select attributes from the menu below. Note that additional attributes to the default ones will slow down the network construction.<br>
         Note also that you can select more than one additional attribute by pushing the Shift key.<br>
         <select multiple name='attributes[]'>
             <option selected value=''>default</option>
             <option value='ECNumber' >EC numbers</option>
             <option value='Label'>Label</option>
             <option value='Equation'>Reaction equation (not for RPAIR network)</option>
             <option value='FORMULA'>Compound formula</option>
         </select>
         </td></tr>
    </table>
  <br>
  <br>
  <br>
  <hr>
  <h2>4. Output
  <a href='help.keggnetworkprovider.html#output'><img src='images/question-mark-button.jpg' width='15' alt='help'></a>
  </h2>
  <br><br>
  <b>Network format&nbsp;&nbsp;</b><br>
    <select name='outFormat'>
        <option value='tab'>tab-delimited</option>
        <option value='gml'>gml</option>
        <option value='dot'>dot</option>
        <option value='visml'>visML</option>
    </select>
  <br><br>
  <table>
  <tr><td>
  <B><a href = 'help.keggnetworkprovider.html#email'>Email (optional)</a></B>
  </td><td>
  <input type='text' name='email' value='$default_email' size='30'></input>
  </td></tr>
  </table>
  <br>
  <hr>
  <table class='formbutton'>
  <TD><input type='submit' name='.submit' value='GO' /></TD>
  <TD><B><A HREF='keggnetworkprovider_form.php?demo=0'>RESET</A></B></TD>
  <TD><B><A HREF='keggnetworkprovider_form.php?demo=1'>DEMO</A></B></TD>
  </form>
  <TD><B><A HREF='help.keggnetworkprovider.html'>MANUAL</A></B></TD>
    <TD><B><A target = '_blank' HREF='".checkNeatTutorial("tutorials/neat_tutorial/KEGG_network_provider.html")."'>TUTORIAL</A></B></TD>
  <TD><B><A HREF='mailto:kfaust@ulb.ac.be'>MAIL</A></B></TD>
  <TD><B><A HREF='help.keggnetworkprovider.html#webservice'>WSDL</A></B></TD>
  </TR></TABLE></ul></ul>");
?>
</body>
</html>