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
  $("#db_id_retrieve").select2({placeholder: "Select identifiers",
                      allowClear: true,
                      dropdownAutoWidth: true,
                      theme: "classic",
                      closeOnSelect: false
                      });
  
  $("#dbs_choice").change(function(demo=0){
                          var db_name = $("#dbs_choice").val();
                          var db_id = $("#db_id");
                          if(typeof db_id === 'undefined'){
                            db_id = $("#db_id_retrieve");
                          }
                          db_id.val("").change();
                          $("#db_choice").select2("data", {id:""});
                          $.ajax({
                                 type: "GET",
                                 url: "getMatrixIds.cgi?dbs_choice=" + db_name + "&mode=db",
                                 dataType: "json"
                            }).done(function(data){
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
                                });
                });
  $("#db_choice").change(function(){
                         var db_id = $("#db_id");
                         if(typeof db_id.val() === 'undefined'){
                            db_id = $("#db_id_retrieve");
                         }
                         
                         $("input[name=db_choice_radio][value=custom_motif_db]").prop("checked",true);
                         
                         var db_name = $("#dbs_choice").val();
                         var collection_name = $("#db_choice").val();
                         
                         $("#wait_ids").attr("style","display:block");
                         
                         if(typeof db_id.val() !== 'undefined'){
                            db_id.val("").change();
                            $("#matrix").val("");
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
                                db_id.html(selectopt);
                                
                                $("#wait_ids").attr("style","display:none");
                                
                                },
                                error: function(){
                                alert("Error matrix db_choice");
                                }
                                });
                         }else{
                            $("#wait_ids").attr("style","display:none");
                         }
                         });
  

  $("#dbs_choice2").select2({placeholder: "Select a database",
                            templateResult: formatState,
                            allowClear: true,
                            theme: "classic"
                            });
  $("#db_choice2").select2({placeholder: "Select a collection",
                           templateResult: formatState,
                           allowClear: true,
                           theme: "classic"
                           });
  
  $("#dbs_choice2").change(function(){
                           db_name = $("#dbs_choice2").val();
                           $("#db_choice2").val("").change();
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
                                  document.getElementById("db_choice2").innerHTML = selectopt;
                                  },
                                  error: function(){
                                  alert("Error matrix dbs_choice");
                                  }
                                  });
                           });
  
  $("#db_choice2").change(function(){
                          if($("#db_choice2").val() != ""){
                          $("#matrix").val("");
                          $("#dbs_choice").val("").change();
                          $("#db_choice").val("").change();
                          db_id = $("#db_id");
                          if(typeof db_id.val() === 'undefined'){
                            db_id = $("#db_id_retrieve");
                          }
                          db_id.val("").change();
                          
                          $("input[name=db_choice_radio][value=motif_collection_all]").prop("checked",true);
                          }
                          });
  
  
  $("#db_id_retrieve").change(function(){
                     db_name = $("#dbs_choice").val();
                     collection_name = $("#db_choice").val();
                     db_id = $("#db_id_retrieve").val();
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
  
  $("#db_id").change(function(){
                     dbs_name = $("#dbs_choice").val();
                     db_name = $("#db_choice").val();
                     db_id = $("#db_id").val();
                     if(db_id == null || db_id == ""){
                     $("#matrix").val("");
                     }
                     if(db_id != null && db_id != ""){
                     $.ajax({
                            type: "GET",
                            url: "getMatrixIds.cgi?dbs_choice=" + dbs_name +"&db_choice=" + db_name + "&db_id=" + db_id + "&mode=retrievepipe",
                            data: {action: "request"},
                            success: function(data){
                            res = data.split("</format>");
                            $("#matrix").val(res[1]);
                            format = "tab";
                            if(res[0] == "tf"){
                            format = "transfac";
                            }
                            $("#matrix_format").val(format);
                            },
                            error: function(){
                            alert("Error matrix values");
                            }
                            });
                     }
                     });
  
  });

function pipto(f){
    var db_name = $("#dbs_choice").val();
    var collection_name = $("#db_choice").val();
    var db_id = $("#db_id_retrieve").val();
    if(db_id != null && db_id != ""){
        $.ajax({
               type: "GET",
               url: "getMatrixIds.cgi?dbs_choice=" + db_name + "&db_choice=" + collection_name + "&db_id=" + db_id + "&mode=retrievepipe",
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


