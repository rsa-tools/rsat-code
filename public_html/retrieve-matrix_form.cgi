#!/usr/bin/env perl
#### this cgi script fills the HTML form for the program retrieve-matrix
BEGIN {
    if ($0 =~ /([^(\/)]+)$/) {
        push (@INC, "$`lib/");
    }
    require "RSA.lib";
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

################################################################
### default values for filling the form
$default{output}="display";
$default{input}="";
$default{id}="";
$default{id_file} = "tab";
$default{table} = 1;

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
    if ($query->param($key)) {
        $default{$key} = $query->param($key);
    }
}

################################################################
### header
&RSA_header("retrieve-matrix", "form");
print "<style>
	table.result {
		border-collapse: collapse;
        table-layout: fixed;
        width: 100%;
	}
	table.result th, table.result td {
		border: 1px solid #cbcbb4;
		padding: 15px;
        word-wrap: break-word;
        font-family: monospace;
        font-size: 12px;
	}
	table.resultlink td, table.resultlink th{
		font-size: 100%;
	}
</style>";

print "<CENTER>";
print "Retrieve the matrices with identifiers.<P>\n";
print "</CENTER>";
print "<BLOCKQUOTE>\n";

################################################################
#### collections
print "<hr>";

print '<link rel="stylesheet" href="css/select2.min.css" /><script src="js/select2.full.min.js"></script>';
print '<style>.select2-results__options {font-size:10px} </style>';
print '<script src="js/retrieve_matrix.js">


 
   </script>';

print '<style> input.select2-search__field { width: 90% !important; } </style>';
print "<table><tr><td style='padding-right:20px'><b>Input</b></td>";

print "<td align='top'>1 - Select a database in list:<br/>";
print ' <select id="dbs_choice" name="dbs_choice" style="width:300px"><option></option>';
## load the various databases that can be compared against
&DisplayMatrixDBchoice("mode"=>"radio");
print "</select>";
print "<td align='top'>2 - Select a collection in list:<br/>";
print ' <select id="db_choice" name="db_choice" style="width:300px"><option></option>';
## load the various databases that can be compared against
print '</select>';


print "</td><td align='top'>";
print "3 - Select one or more matrix identifiers:<br/>";
print " <select id='db_id_retrieve' style='width:300px' multiple='multiple'></select>";
print '<div style="clear:both;"></div>';
print "</td><td><div id='wait_ids' style='display:none'><image src='images/wait.gif' style='width:50px;height:50px' /></div></td></tr></table>";

print '<br/><a class="inline" href="#matrix_descr""> View matrix descriptions & download full collections</a> <br/>';
print "<div style='display:none'><div id='matrix_descr'>";
&DisplayMatrixDBchoice("mode" => "list");
print "</div></div></div>";

print "<hr>";
####### useful link
print '<script>
function setDemo(){
	setDemo1();
	$("#wait_ids").attr("style","display:block");
	setTimeout(setDemo2, 2000);
}
function setDemo1(){
    $("#dbs_choice").val("Jaspar").change();
}
function setDemo2(){
    $("#db_choice").val("jaspar_core_nonredundant_vertebrates").change();
    setTimeout(setDemo3, 2000);
}
function setDemo3(){
    $("#db_id_retrieve").val(["MA0019.1", "MA0031.1"]).change();
    $("#wait_ids").attr("style", "display:none");
}

function reset(){
    $("#dbs_choice").val("").change();
}

function sendemail(){
    email = $("#user_email").val();
    db_name = $("#dbs_choice").val();
    collection_name = $("#db_choice").val();
    db_id = $("#db_id_retrieve").val();
    if(db_id != "" && db_id != null){
        $.ajax({
            type: "GET",
            url:"getMatrixIds.cgi?dbs_chocie=" + db_name + "&db_choice=" + collection_name + "&db_id=" + db_id + "&mode=retrieveemail&user_email=" + email,
            success: function(data){
                $("#sendemailmsg").html(data);
                document.getElementById("sendemailmsg").style.display = "block";
            }
        });
    }
}

</script>';

print '<br/><button type="reset" onclick="reset()">RESET</button>&nbsp;<button type="button" onclick="setDemo()">DEMO</button>';

print "&nbsp;<b><A class='iframe' HREF='help.retrieve-matrix.html'>MANUAL</A>&nbsp; ";


### send results by email
print "<div id='email' style='display:none'><hr>";
print " <b>Send to my email</b> <input type='text' id='user_email' size='30' /><button id='sendemail' onclick='sendemail()'>SEND</button>" unless ($ENV{mail_supported} eq "no");
print "</div>";

######### result

print "<br/><br/><hr><h2>Result</h2>The selected matrix will be displayed in the table below<br/><br/>";
print "<div id='sendemailmsg'></div>";
print "<div id='outputurl'></div>";
print "<div id='result'></div>";

### prepare data for piping
print "<div id='piping' style='display:none'></div><hr>";

print $query->end_html;

exit(0);
