#!/usr/bin/perl
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
print '<script type="text/javascript">

    function formatState(state){
        if(!state.id){return $("<div align=\'center\' style=\'text-transform:uppercase;background-color:lightgray;border:1px solid #eee;border-radius:7px;padding:8px\'><b>" + state.text + "</b></div>");}
        var $state = $("<span><i class=\'fa fa-circle fa-lg\' style=\'color:" + state.element.className + "\'></i></span>&nbsp;&nbsp;&nbsp;" + state.text + "</span>");
        return $state;
    }
    $(function(){
        $(".inline").colorbox({inline:true,width:"70%"});
        // turn the element to select2 select style
        $("#db_choice").select2({placeholder: "Select a collection",
            templateResult: formatState,
            allowClear: true,
            theme: "classic"
        });
    
        $("#db_id").select2({placeholder: "Select identifiers",
            allowClear: true,
            dropdownAutoWidth: true,
            theme: "classic"
        });
    
        $("#db_choice").change(function(){
            $("#db_id").val("").change();
            $("#matrix").val("");
            db_name = $("#db_choice").val();
            $.ajax({
                type: "GET",
                url: "getMatrixIds.cgi?db_choice=" + db_name,
                dataType: "json",
                data: {action: "request"},
                success: function(data){
                    var res = data.entries;
                    if(res.length != 0){ document.getElementById("identifier").style.display = "block"; }
                    var selectopt = "<option></option>";
                    for(i = 0; i < res.length; i++){
                        selectopt += "<option value=\'" + res[i].id + "\'>" + res[i].idac + "</option>";
                    }
                    document.getElementById("db_id").innerHTML = selectopt;
                },
                error: function(){
                    alert("Error matrix");
                }
            });
        });
    $("#db_id").change(function(){
        db_name = $("#db_choice").val();
        db_id = $("#db_id").val();
        output = $("input[name=output]:checked").val();
        if(db_id != null && db_id != ""){
            $.ajax({
                type: "GET",
                dataType:"json",
                url: "getMatrix.cgi?db_choice=" + db_name + "&db_id=" + db_id + "&mode=retrieve",
                data: {action: "request"},
                success: function(data){
                    res = data.entries;
                    outputhtml = "</br><table class=\"result\">";
                    for(i = 0; i < res.length; i++){
                        outputhtml += "<tr><td>" + res[i].info + "</td><td><pre>" + res[i].all + "</pre></td></tr>";
                    }
                    outputhtml += "</table>";
                    $("#result").html(outputhtml);
                    
                    outputfile = "<table class=\"resultlink\"><tr><th>Content</th><th>URL</th></tr>";
                    outputfile += "<tr><td>Input file</td><td><a href=\"" + data.inputfile + "\" target=\"_blank\">" + data.inputfile + "</a></td></tr>";
                    outputfile += "<tr><td>Output file</td><td><a href=\"" + data.resultfile + "\" target=\"_blank\" id=\"resultfile\">" + data.resultfile + "</a></td></tr>";
                    outputfile += "</table>";
                    $("#outputurl").html(outputfile);
                },
                error: function(){
                    alert("Error matrix");
                }
            });
        }
    });
    
    $("#email").change(function(){
        if($(this).is(":checked")){
            document.getElementById("sendemail").style.display = "inline";
            document.getElementById("sendemailmsg").style.display = "block";
        }else{
            document.getElementById("sendemail").style.display = "none";
            document.getElementById("sendemailmsg").style.display = "none";
        }
    });
    });
   </script>';


print "<table><tr><td style='padding-right:20px'><b>Input</b></td>";

print "<td align='top'>1 - Select a collection in list:<br/>";
print ' <select id="db_choice" name="db_choice" style="width:300px"><option></option>';
## load the various databases that can be compared against
&DisplayMatrixDBchoice_select2("mode"=>"radio");
print '</select>';
print '<br/><a class="inline" href="#matrix_descr""> View matrix descriptions</a> <br/>';
print "<div style='display:none'><div id='matrix_descr'>";
&DisplayMatrixDBchoice_select2("mode" => "list");
print "</div></div></div>";

print "</td><td align='top'>";

print "<div id='identifier'><div style='float:left;margin-left:27px;font-size:11px;'>2 - Select 1 or more matrix identifiers:<br/>";
print " <select id='db_id' style='width:300px' multiple='multiple'></select><br/>&nbsp;";
print "</div></div>";
print '<div style="clear:both;"></div>';

print "</td></tr></table>";


print "<br/><hr>";
### send results by email or display on the browser
print "<input type='checkbox' id='email' />Email <input type='text' id='user_email' size='30'/>" unless ($ENV{mail_supported} eq "no");

print "<button style='display:none' id='sendemail' onclick='sendemail()'>SEND EMAIL</button>&nbsp;";
####### useful link
print "<script>
function setDemo(){
    \$.ajax({
        url:setDemo1(),
        success:function(){
            setDemo2();
        }
    });
}
function setDemo1(){
    \$('#db_choice').val('jaspar_core_nonredundant_vertebrates').change();
    
}
function setDemo2(){
    \$('input[name=output][value=display]').prop('checked', true);
    \$('#db_id').val(['MA0019.1', 'MA0031.1']).change();
}

function reset(){
    \$('#db_choice').val('').change();
    \$('input[name=output][value=server]').prop('checked', true);
}

function sendemail(){
    email = \$('#user_email').val();
    db_name = \$('#db_choice').val();
    db_id = \$('#db_id').val();
    \$.ajax({
        type: 'GET',
        url:'getMatrix.cgi?db_choice=' + db_name + '&db_id=' + db_id + '&output=email&user_email=' + email,
        success: function(data){
            \$('#sendemailmsg').html(data);
            document.getElementById('sendemailmsg').style.display = 'block';
        }
    });
}
</script>";
print '<br/><button type="reset" onclick="reset()">RESET</button>&nbsp;<button type="button" onclick="setDemo()">DEMO</button>';


print "&nbsp;&nbsp;<b><A class='iframe' HREF='help.retrieve-seq.html'>MANUAL</A>&nbsp; ";
print "<A HREF='htmllink.cgi?title=RSAT : Tutorials&file=tutorials/tut_retrieve-seq.html'>TUTORIAL</A> &nbsp;";
print "<A HREF='mailto:Jacques.van-Helden\@univ-amu.fr'>MAIL</A></b>";

######### result
print "<br/><br/><hr><h2>Result</h2>The selected matrix will be displayed in the table below<br/><br/>";
print "<div id='outputurl'></div>";
print "<div id='result'></div>";
print "<div id='sendemailmsg'></div>";
print $query->end_html;

exit(0);

