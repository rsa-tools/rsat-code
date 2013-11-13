<html>
<head>
   <title>NAT - KEGG network provider</title>
   <link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">
</head>
<body class="form">
<?php
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
  $kgmlVersionFile = $rsat_main.'/public_html/data/KEGG/kgmlVersion.txt';
  $kgmlVersion = storeFile($kgmlVersionFile);

  # demo
  $demo = $_REQUEST['demo'];
  if ($demo == 1) {
     # Comment: in higher KEGG versions, the three other yeast species are no longer available as KGML files
     if($kgmlVersion < 50.0){
        $default_organisms = "dsba/dsmi/dspd/spo/sce";
     }else{
        $default_organisms = "spo/sce";
     }
  }else if ($demo == 2){
     $default_organisms = "reference";
     $default_excludedcompounds = "C01371/C06058/C01872/C15521/C02403/C00722/C01372/C01706/C06015/C00069C00226/C00339/C00432/C02525/C02893/C03130/C03485/C03491/C01813/C01891/C00947/C00709/C02019/C02146/C00071/C00538/C00609/C01928/C00193/C00287/C01629/C01016/C02434/C01883/C03317/C02391/C00060/C00161/C00347/C00539/C00706/C00893/C01609/C01664/C02136/C02196/C02311/C02324/C02580/C05360/C00972/C03096/C03252/C03547/C00241/C00484/C01658/C01976/C02058/C01663/C02179/C02244/C16072/C00726/C01887/C16073/C01812/C01318/C03122/C00045/C00405/C00151/C01021/C01076/C03253/C05167/C06699/C06700/C02847/C02853/C03135/C03148/C04245/C03523/C04313/C05710/C06378/C02850/C16116/C00737/C00991/C01370/C01547/C01599/C01383/C00030/C00028/C00050/C01815/C02148/C02169/C02285/C00162/C00638/C01229/C02403/C03395/C03547/C04215/C04340/C05102/C05996/C00377/C06479/C15506/C15507/C00370/C00012/C00098/C00107/C00316/C01960/C02179/C01579/C02790/C00420/C01923/C00182/C00403/C02056/C00017/C00046/C00066/C00240/C00731/C00960/C02447/C03638/C03880/C01396/C01522/C02073/C02211/C02342/C02521/C03959/C05156/C05157/C00039/C00271/C00434/C00578/C00871/C02128/C02270/C02374/C02410/C02840/C03310/C03484/C03959/C04395/C04249/C00718/C01003/C02707/C01356/C00678/C00865/C05005/C01915/C07305/C02843/C03129/C01939/C03667/C06061/C00174/C01373/C00701/C06708/C05426/C10905/C11349/C03376/C03377";
  }

  title('KEGG network provider');
  echo ("<center>Obtain organism-specific KEGG networks from <a href=\"http://www.genome.ad.jp/kegg/pathway.html\" target=\"_blank\">KEGG PATHWAY</a> version $kgmlVersion. Click on <a href='help.keggnetworkprovider.html'><img src='images/question-mark-button.jpg' width='15' alt='help'></a> for help.<br>
  Web service and interface by Karoline Faust.</center>\n");
  echo ("<form method='post' action='keggnetworkprovider.php' enctype='multipart/form-data'>
  <hr>
  <h2>1. Input <a href='help.keggnetworkprovider.html#input'><img src='images/question-mark-button.jpg' width='15' alt='help'></a>
  </h2>
  <br>");
  if($demo == 1){
     # Comment: in higher KEGG versions, the three other yeast species are no longer available as KGML files
    if($kgmlVersion > 50.0){
        demo("To demonstrate the KEGG network provider, we are constructing the metabolic network from two yeast species:<br>
        spo (Schizosaccharomyces pombe) and sce (Saccharomyces cerevisiae).");
    }else{
        demo("To demonstrate the KEGG network provider, we are constructing the metabolic network from five yeast species:<br>
        dsba (Saccharomyces bayanus), dsmi (Saccharomyces mikatae), dspd (Saccharomyces paradoxus), spo (Schizosaccharomyces pombe)
        and sce (Saccharomyces cerevisiae).");
    }
   }else if($demo == 2){
   demo("We demonstrate how the KEGG network provider can be used to
   filter the KEGG network. KEGG LIGAND contains many problematic entries, such as
   generic compounds (e.g. C00030, which represents a reduced receptor)
   and unbalanced or duplicated reactions
   (see F&eacute;lix and Valiente, Biomol. Eng. 24, 327-335 (2007),
   M.G. Poolman et al., IEE Proc.-Syst. Biol. 153, 379-384 (2006) and
   Ott and Vriend, BMC Bioinformatics 7:517 (2006)).
   These entries make path finding and pathway inference more difficult. In this
   demo, we will construct the reference KEGG RPAIR network, from which we will remove a selection of
   generic compounds.");
  }
  echo("
  <table>
  <tr><td>
         Enter full species names as given in KEGG or KEGG abbreviations for organisms, separated by '/' in text field below.<br>
         Example 1: eco/ecp/ecs<br>
         Example 2: Aspergillus flavus/Aspergillus niger/afm<br>
         Example 3: reference<br>
         <a href='data/KEGG/Kegg_organisms_list.txt'><b>Click here for the full list of KEGG organism names and abbreviations</b></a>
  </td></tr>
  <tr><td>
        <input type='text' name='organisms' value='$default_organisms' size=100></input>
   </td></tr>
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
  <br>
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
         <tr><td><B>Directed network</B><input type='checkbox' name='directed' value='on' checked></input></td></tr>
         ");
         if($demo == 2){
          echo("
         <tr><td><B>RPAIR network</B>&nbsp;&nbsp;&nbsp;<input type='checkbox' name='rpair' value='on' checked></input><br>
         ");
         }else{
         echo("
         <tr><td><B>RPAIR network</B>&nbsp;&nbsp;&nbsp;<input type='checkbox' name='rpair' value='on'></input><br>
         ");
         }
         echo("
         For the conversion of reactions into reactant pairs, a file assembled from KEGG LIGAND is used, which
	     is available <a href='$host/rsat/data/KEGG/rpairs.tab'>here</a>.
         </td></tr>
         <tr><td><b>Reaction treatment:<b><br>
         <table>
         	<tr><td>&nbsp;&nbsp;&nbsp;&nbsp;</td><td>all reactions reversible&nbsp;&nbsp;&nbsp;
         	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
         	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
         	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
         	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
         	<input type='radio' name='keepIrre' value='false' checked='checked'></input></td></tr>
         	<tr><td>&nbsp;&nbsp;&nbsp;&nbsp;</td><td>keep irreversible reactions (only possible in directed networks!)
         	<input type='radio' name='keepIrre' value='true'></input></td></tr>
         </table>
         </td></tr>
         <tr><td><br></td></tr>
         <tr><td><b>Optional:</b> The network is constructed with some
         <a href='help.keggnetworkprovider.html#attributes'>default attributes</a>.
         You can add more attributes by selecting them from the menus below, but be aware that
         this will slow down the network construction.
         Note also that you can select more than one additional attribute by pushing the Shift key.<br>
         <table>
         <tr>
         <td>Compound attributes<br>
         <select multiple name='compoundattributes[]'>
             <option value='Label'>Name</option>
             <option value='FORMULA'>Formula</option>
         </select>
         </td>
         <td>
         &nbsp;&nbsp;
         </td>
         <td>Reaction/RPAIR attributes<br>
         <select multiple name='reactionattributes[]'>
             <option value='ECNumber' >EC numbers</option>
             <option value='Equation'>Equation (not for RPAIR network)</option>
         </select>
         </td>
         </tr>
         </table>
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
  <TD><B><A HREF='keggnetworkprovider_form.php?demo=1'>DEMO1</A></B></TD>
  <TD><B><A HREF='keggnetworkprovider_form.php?demo=2'>DEMO2</A></B></TD>
  </form>
  <TD><B><A HREF='help.keggnetworkprovider.html'>MANUAL</A></B></TD>
    <TD><B><A target = '_blank' HREF='".checkNeatTutorial("tutorials/neat_tutorial/KEGG_network_provider.html")."'>TUTORIAL</A></B></TD>
  <TD><B><A HREF='mailto:kfaust@ulb.ac.be'>MAIL</A></B></TD>
  <TD><B><A HREF='neat_webservices.php'>WSDL</A></B></TD>
  </TR></TABLE></ul></ul>");
?>
</body>
</html>