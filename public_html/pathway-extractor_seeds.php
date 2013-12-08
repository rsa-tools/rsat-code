<html>
<head>
<title>RSAT - Pathway-Extractor</title>
<link rel="stylesheet" type="text/css" href = "main_grat.css" media="screen">

<script language="javascript" type="text/javascript" src="TableFilter/tablefilter_all_min.js"></script>
<script src="TableFilter/sortabletable.js" language="javascript" type="text/javascript"></script>
<script src="TableFilter/tfAdapter.sortabletable.js" language="javascript" type="text/javascript"></script>
<script src="pathwayextractor.js" language="javascript" type="text/javascript"></script>
<script language="javascript" type="text/javascript">
var table3Filters = {
		sort: true,
		sort_config: { sort_types:['String', 'string', 'None', 'none', 'None'] },
		on_keyup: true,
		on_keyup_delay: 1000,
		highlight_keywords: true,
//		sort_select: true,
//		col_1: 'select',
		status_bar: true,
		display_all_text: " [ Show all ] ",
		rows_counter: true,
		btn_reset: true,
//		alternate_rows: true,
		btn_reset_text: "Clear",
		col_width: ['120px','80px','250px',null,'250px'],
		popup_filters: true
}
var seedstableFilters = {  
        sort: true,
        sort_config: { sort_types:['none', 'string', 'string', 'string', 'string', 'string'] },
		col_0: "none",  
//        grid_enable_cols_resizer: false,
        popup_filters: true,
        btn_reset: true,
		on_keyup: true,
		on_keyup_delay: 1000,
		highlight_keywords: true,        
       alternate_rows: true,
       help_instructions_html:'<ul><li>Click on column header for sorting</li><li>Click on <img alt="Column filter" src="TableFilter/TF_Themes/icn_filter.gif"> for filtering</li><li>Click on <img alt="Reset filter" src="TableFilter/TF_Themes/btn_clear_filters.png"> to reset the filters</li><li> Click on \'Select all\' to select all visible seeds</li></ul>'
//        grid_width: '900px'
    }  
function init(){
//	alert("coucou");
	//document.getElementById('hourglass').hidden="true";
	var theseedstable = new TF("seedstable",1,seedstableFilters);  
	theseedstable.AddGrid();
	var elem = document.getElementById("inf_seedstable");
//	var content="<a class=\"Button\" href=\"javascript:selectAllVisiblecb('chkID',0);\">Select All</a>&nbsp<a class=\"Button\" href=\"javascript:selectAllVisiblecb('chkID',1 );\">Inverse All</a>";
//	elem.innerHTML="<div id='myleftdiv' class='ldiv'>"+content+"</div>"+elem.innerHTML;  
}
function buildarrayforseeds(form)
{ 	var max=form.elements.length;
	
	var elements = document.getElementsByTagName("TABLE");
	//var seedsinput = document.getElementById("seeds");
	var max=form.elements.length;
	//var seeds="";
	var j = 0;
	for (var i=0;i<max;i++) {
		temp=form.elements[i];
		if (temp.name=="chkID" && temp.checked) {
			//seedsinput[j++] = temp.value;
			//alert(temp.value);
			addHiddenArrayInputElement(form,"seeds",j,temp.value);
			j++;
		}
	}
	if (j++>0) return true;
	
}
function addHiddenArrayInputElement(form,name,suffix,value) {
	  var newelem = document.createElement('input');
	  newelem.setAttribute('name',name+"["+suffix+"]");
	  newelem.setAttribute('id',name);
	  newelem.setAttribute('type',"hidden");
	  newelem.setAttribute('value',value);
	  form.appendChild(newelem);

}
function enum_formforseeds(form)
{
	var max=form.elements.length;
	var seeds="";
	for (var i=0;i<max;i++) {
		temp=form.elements[i];
		if (temp.name=="chkIDS" && temp.checked) {
			seeds+=temp.value + ";";
		}
	}
	return seeds;
}

function submitSel(){
	document.forms['actionForm'].submit();
}

function setChecked(e,b){
    for(var i=0; i<e.length; i++)
            e[i].checked = b
}

function toggleChecked(e){
    for(var i=0; i<e.length; i++)
            e[i].checked = !e[i].checked
}
</script>
</head>
   <body class="results" onload="javascript:init()"> 
   <!-- div id=hourglass name=hourglass><img src="images/animated_hourglass.gif" alt="Please Wait!!"></img></div-->
<?php
// Load RSAT configuration
   require('functions.php');

// print_r($properties);
UpdateLogFile("rsat","","");

// Import variables with prefix fs_ from form
import_request_variables('P','fs_');

//print_r($_POST);
// print_r($_FILES);


$basedir=$properties['RSAT']."/public_html/data/metabolic_networks";
#$cmd = $properties['RSAT'].'/perl-scripts/search4metapath-seeds';
$cmd = $properties['RSAT'].'/perl-scripts/network-match-nodes';

////////////////////////////////////////////////////////////////
//Print <h3>
echo "<H3><a href='".$properties['rsat_www']."'>RSAT</a> - Pathway-Extractor <BR> Step 2 : Mapping Seeds - results</H3>";

////////////////////////////////////////////////////////////////
// Check arguments
$errors = 0;

// Check that gnn has been specified
if ($fs_gnn == "none" or $fs_gnn == "" ) {
  error( "You forgot to select the genome source for the mapping.");
  $errors=1;
} else {
   $GERdir = "$basedir/GER_files/$fs_gnn";
   $argumenta .= " -p -v 0 -gnn ".$GERdir."/".$fs_gnn."-gene_ec.tab";
 }
// Check that gnn has been specified
if ($fs_network == "none" or $fs_network == "" ) {
  error( "You forgot to select network.");
  $errors=1;
} else {
  $neworkfilepattern= $fs_network;
  $neworkfilepattern = eregi_replace("_v.*","",$fs_network); 
  $neworkfilepattern = "$basedir/networks/$neworkfilepattern/$fs_network";
  $networknodenames = $neworkfilepattern."-metab-node_names.tab";
  $argumenta .= " -nnn ".$networknodenames;
}
//Check syntax of email address (ensure the texte netered in email box was an email address)
if($fs_output =="email") {
  if (!preg_match("#^[\wàáâãäåæçèéêëìíîïðñòóôõöøùúûüý._\-]+@([a-z]+.)+[a-z]{2,4}$#", $fs_user_email)) {
     error( "Email not valid");
     $errors=1;
  }
 }

////////////////////////////////////////////////////////////////
// Check that the BED has been specified once and only once

// Bed data pasted in text area
#$fs_seeds = $_POST["seeds"];

if ($fs_seeds == "") {
  error( "No seeds!");
  $errors=1;
}else {
  $seeds = preg_replace( "#\r\n|\r|\t#", "\n", $_POST["seeds"] );
  #$array_line = explode("\n",$fs_seeds);
}



///////////////////////////////////////////
// Run commands
if ($errors == 0) { 
  $cmd .= $argumenta;

// SHOW COMMAND and SEEDS can be commented
?> 
<fieldset><legend class="clickable" onclick="javascript:collapseAll(document.getElementById('showcommand'));">
      <b><span id="showcommandSign" class="sign"><img src="images/arrow_box.gif" alt="-"></span></b><u>Show Command</u></legend>
      <div id="showcommand" style="display: none;">
<?php      
  info("Command : ".str_replace($properties['RSAT'], '$RSAT', $cmd));
  info("Seeds : ".$seeds);
?>
</div>
</fieldset>
<?php
flush(); 
  //  Run the command
$descriptorspec = array(
   0 => array("pipe", "r"),
   1 => array("pipe", "w"),
   2 => array("file", "/tmp/error-output.txt", "a")
);
$process = proc_open($cmd, $descriptorspec, $pipes);
  
if (is_resource($process)) {

    //row2xfdf is made-up function that turns HTML-form data to XFDF
    fwrite($pipes[0], $seeds);
    fclose($pipes[0]);

    $ouput_content = stream_get_contents($pipes[1]);
    fclose($pipes[1]);
//	info("ouput_content : ".$ouput_content);
    $return_value = proc_close($process);

 //   header('Content-type: text/plain');
 //   header('Content-Disposition: attachment; filename="output.tab"');
//	info("Mapping : ".$return_value);
    //    echo $return_value;
}  

//if ($return_value = 0){
	
$str = rtrim($ouput_content);
//class=sortable 

?>
<form method='post' action='pathway-extractor.php' enctype='multipart/form-data' onSubmit="return buildarrayforseeds(this);">
<!--input type=hidden id="seeds" name="seeds[]" value=""-->
<?php
echo converttab2table ($str,"name=seedstable id=seedstable");
?>
<a class="Button" href="javascript:selectAllVisiblecb('chkID',0);">Select All</a>
<a class="Button" href="javascript:selectAllVisiblecb('chkID',1 );">Inverse All</a><BR><BR>
<?php 
echo "<hr>"; 
$req = array_merge($_GET, $_POST);
foreach($req as $key=>$val)
{ 
  if ($key!="submit"){
  	if ($key != "seeds"){
  		echo	"<input type=hidden name=".$key." value=".$val.">\n";
  	}
  }
}      
}	

?>
<!-- input type=button name=test id=test onClick="enum_formforseeds(this.form)"> -->
<input type="submit" name="submit" value="GO" />
 </form>
  </body>
</html>
