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
print '<script>

    function formatState(state){
        if(!state.id){return $("<div align=\'center\' style=\'text-transform:uppercase;background-color:lightgray;border:1px solid #eee;border-radius:7px;padding:8px\'><b>" + state.text + "</b></div>");}
        var $state = $("<span><i class=\'fa fa-circle fa-lg\' style=\'color:" + state.element.className + "\'></i></span>&nbsp;&nbsp;&nbsp;" + state.text + "</span>");
        return $state;
    }
    $(function(){
        $(".inline").colorbox({inline:true,width:"70%"});
        // turn the element to select2 select style
        $("#dbs_choice").select2({placeholder: "Select a database",
            templateResult: formatState,
            allowClear: true,
            theme: "classic"
        });
        $("#db_choice").select2({placeholder: "Select a collection",
            templateResult: formatState,
            allowClear: true,
            theme: "classic"
        });
        $("#db_id").select2({placeholder: "Select identifiers",
            allowClear: true,
            dropdownAutoWidth: true,
            theme: "classic",
            closeOnSelect: false
        });
        
        $("#dbs_choice").change(function(){
            db_name = $("#dbs_choice").val();
            $("#db_id").val("").change();
            $("#matrix").val("");
            $.ajax({
                type: "GET",
                url: "getMatrixIds.cgi?dbs_choice=" + db_name + "&mode=db",
                dataType: "json",
                data: {action: "request"},
                success: function(data){
                    var res = data.entries;
                    var selectopt = "<option></option>";
                    for(i = 0; i < res.length; i++){
                        selectopt += "<optgroup label=\'" + res[i].category + "\'>";
                        for(j = 0; j < res[i].data.length; j++){
                            selectopt += "<option value=\'" + res[i].data[j].name + "\' class=\'" + res[i].color + "\'>" + res[i].data[j].descr + "</option>";
                        }
                        selectopt += "</optgroup>";
                    }
                    document.getElementById("db_choice").innerHTML = selectopt;
                },
                error: function(){
                    alert("Error matrix dbs_choice");
                }
            });
        });
        $("#db_choice").change(function(){
            $("#db_id").val("").change();
            $("#matrix").val("");
            db_name = $("#dbs_choice").val();
            collection_name = $("#db_choice").val();
            $.ajax({
                type: "GET",
                url: "getMatrixIds.cgi?dbs_choice=" + db_name + "&db_choice=" + collection_name + "&mode=id",
                dataType: "json",
                data: {action: "request"},
                success: function(data){
                    var res = data.entries;
                    var selectopt = "<option></option>";
                    for(i = 0; i < res.length; i++){
                        selectopt += "<option value=\'" + res[i].id + "\'>" + res[i].idac + "</option>";
                    }
                    document.getElementById("db_id").innerHTML = selectopt;
                },
                error: function(){
                    alert("Error matrix db_choice");
                }
            });
        });
        $("#db_id").change(function(){
            db_name = $("#dbs_choice").val();
            collection_name = $("#db_choice").val();
            db_id = $("#db_id").val();
            output = $("input[name=output]:checked").val();
            if(db_id != null && db_id != ""){
                $.ajax({
                    type: "GET",
                    dataType:"json",
                    url: "getMatrixIds.cgi?dbs_choice=" + db_name +"&db_choice=" + collection_name + "&db_id=" + db_id + "&mode=retrieveweb",
                    data: {action: "request"},
                    success: function(data){
                        res = data.entries;
                        outputhtml = "</br><div style=\"max-width:1200px\"><table class=\"result\">";
                        for(i = 0; i < res.length; i++){
                            outputhtml += "<tr><td>" + res[i].info + "</td><td>" + res[i].all + "</td></tr>";
                        }
                        outputhtml += "</table>";
                        $("#result").html(outputhtml);
                        
                        outputfile = "<table class=\"resultlink\"><tr><th>Content</th><th>URL</th></tr>";
                        outputfile += "<tr><td>Input file</td><td><a href=\"" + data.inputfile + "\" target=\"_blank\">" + data.inputfile + "</a></td></tr>";
                        outputfile += "<tr><td>Output file</td><td><a href=\"" + data.resultfile + "\" target=\"_blank\" id=\"resultfile\">" + data.resultfile + "</a></td></tr>";
                        outputfile += "</table></div>";
                        $("#outputurl").html(outputfile);
                        $("#sendemailmsg").html("");
                        
                        document.getElementById("piping").style.display = "block";
                        var piphtml = "<HR SIZE = \"3\" />\
                        <TABLE class=\'nextstep\'>\
                        <tr><td colspan = 4><h3 style=\'background-color:#0D73A7;color:#D6EEFA\'>Next step</h3></td></tr>\
                        <tr valign=\"top\" align=\"center\">\
                            <th align=center>\
                                <font size=-1>Matrix tools</font>\
                            </th>\
                            <td align=\"center\" style=\'font-size:100%\'>\
                                <input type=\"button\" onclick=\"pipto(\'convert-matrix\')\" value=\"convert-matrix\" /><br/>\
                                Convert position-specific<br/>scoring matrices (PSSM)\
                            </td>\
                            <td align=\"center\" style=\'font-size:100%\'>\
                                <input type=\"button\" onclick=\"pipto(\'compare-matrices\')\" value=\"compare-matrices\" /><br/>\
                                Compare two collections of<br/>position-specific scoring matrices\
                            </td>\
                            <td align=\"center\" style=\'font-size:100%\'>\
                                <input type=\"button\" onclick=\"pipto(\'matrix-clustering\')\" value=\"matrix-clustering\" /><br/>\
                                Identify groups (clusters) of similarities<br/>between a set of motifs and align them.\
                            </td>\
                        </tr>\
                        <tr valign=\"top\" align=\"center\">\
                            <th align=center>\
                                <font size=-1>Pattern matching</font>\
                            </th>\
                            <td align=center style=\'font-size:100%\'>\
                                <input type=\"button\" onclick=\"pipto(\'matrix-scan\')\" value=\"matrix-scan\" /><br/>\
                                Scan a DNA sequence with a profile matrix\
                            </td>\
                            <td align=center style=\'font-size:100%\'>\
                                <input type=\"button\" onclick=\"pipto(\'matrix-scan-quick\')\" value=\"matrices-scan(quick)\" /><br/>\
                                Scan a DNA sequence with a profile matrix - quick version\
                            </td>\
                        </tr></TABLE>";
                        $("#piping").html(piphtml);
                        
                        document.getElementById("email").style.display = "inline";
                    },
                    error: function(){
                        alert("Error matrix id_choice");
                    }
                });
            }
        });
        
    });
        
    function pipto(f){
        db_name = $("#dbs_choice").val();
        collection_name = $("#db_choice").val();
        db_id = $("#db_id").val();
        if(db_id != null && db_id != ""){
            $.ajax({
                type: "GET",
                url: "getMatrix.cgi?dbs_choice=" + db_name + "&db_choice=" + collection_name + "&db_id=" + db_id + "&mode=retrievepipe",
                data: {action: "request"},
                success: function(data){
                    res = data.split("</format>");
                    format = (res[0] == "tf") ? "transfac" : "tab";
                    matrix = res[1];
                    $("form").remove();
                    var html = "<form id=\"dynForm_" + f + "\" action=\"" + f + "_form.cgi\" method=\"post\"><input type=\"hidden\" name=\"matrix_format";
                    if(f == "matrix-clustering"){
                        html += "1\" value=\"" + format + "\"><input type=\"hidden\" name=\"matrix";
                    }else{
                        html += "\" value=\"" + format + "\"><input type=\"hidden\" name=\"matrix";
                    }
                    
                    if(f == "matrix-clustering"){
                        html += "1\" value=\"" + matrix + "\"><input type=\"hidden\" name=\"html_title\" value=\"from retrieve-matrix\"></form>";
                    }else{
                        html += "\" value=\"" + matrix + "\"></form>";
                    }
                    document.getElementById("piping").innerHTML += html;
                    document.getElementById("dynForm_" + f).submit();
                }
            });
        }
    }
    
    
    function setDemo(){
        $.ajax({
            url:setDemo1(),
            success:function(){
                $.ajax({
                    url:setDemo2(),
                    success:function(){
                        setDemo3();
                    }
                });
                
            }
        });
    }
    function setDemo1(){
        $("#dbs_choice").val("Jaspar").change();
    }
    function setDemo2(){
        $("#db_choice").val("jaspar_core_nonredundant_vertebrates").change();
        
    }
    function setDemo3(){
        $("input[name=output][value=display]").prop("checked", true);
        $("#db_id").val(["MA0019_1", "MA0031_1"]).change();
    }
    
    function reset(){
        $("#dbs_choice").val("").change();
    }
    
    function sendemail(){
        email = $("#user_email").val();
        db_name = $("#dbs_choice").val();
        collection_name = $("#db_choice").val();
        db_id = $("#db_id").val();
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

print '<style> input.select2-search__field { width: 90% !important; } </style>';
print "<table><tr><td style='padding-right:20px'><b>Input</b></td>";

print "<td align='top'>1 - Select a database in list:<br/>";
print ' <select id="dbs_choice" name="dbs_choice" style="width:300px"><option></option>';
## load the various databases that can be compared against
&DisplayMatrixDBChoice("mode"=>"radio");
print "</select>";
print "<td align='top'>2 - Select a collection in list:<br/>";
print ' <select id="db_choice" name="db_choice" style="width:300px"><option></option>';
## load the various databases that can be compared against
print '</select>';

print "</td><td align='top'>";
print "3 - Select one or more matrix identifiers:<br/>";
print " <select id='db_id' style='width:300px' multiple='multiple'></select>";
print '<div style="clear:both;"></div>';
print "</td></tr></table>";

print '<br/><a class="inline" href="#matrix_descr""> View matrix descriptions & download full collections</a> <br/>';
print "<div style='display:none'><div id='matrix_descr'>";
&DisplayMatrixDBChoice("mode" => "list");
print "</div></div></div>";

print "<hr>";
####### useful link
print "<script>

</script>";
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
