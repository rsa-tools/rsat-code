/* Functions used to generate the pathway extractor form */

/*  */
function collapseAll(id)
{
    var _style = id.style.display;
    var sign = document.getElementById(id.id.toString() + "Sign");
    var s="";
    if (_style == "none") {
	_style = "block";
	s="<img src=\"images/arrow_box_down.gif\" alt=\"-\">";
    }
    else {
	_style = "none";
	s="<img src=\"images/arrow_box.gif\" alt=\"+\">";
    }
    id.style.display=_style;
    sign.innerHTML=s;
}

/*  */
function showSelected(){
    var selector = document.getElementById('gmx');
    showall(selector.options[selector.selectedIndex].value);
}

/*  */
function showall(name){
    var cont = document.getElementById('parameters');
    for (i=0; i<cont.childNodes.length; i++){
	if (cont.childNodes[i].hasAttributes()&&(cont.childNodes[i].tagName == 'SPAN'|| cont.childNodes[i].tagName == 'DIV'|| cont.childNodes[i].tagName == 'A')){
	    if (cont.childNodes[i].getAttribute('Name') == name ){
		if (cont.childNodes[i].tagName != 'DIV')cont.childNodes[i].setAttribute("style","display: ");
	    }else cont.childNodes[i].setAttribute("style","display: none");
	}
    }

}

/*  */
function deleteUploadedFiles(parent){
    var dir = document.getElementById('localcls').value.substring(0, document.getElementById('localcls').value.lastIndexOf('/'));
    //alert(parent.href);
    if(parent.href.match(".*[?].*")){
	//alert(parent.href+"with & :"+dir);
	parent.href = parent.href + "&dir2del=" + dir;
    }else{
	//alert(parent.href+"with ? :"+dir);
	parent.href = parent.href + "?dir2del=" + dir;
    }
}

/*  */
function toggleAll(id){
    var _style = id.style.display;
    var sign = document.getElementById(id.id.toString() + "Sign");
    var s="";
    //	_style = (_style == "none") ? "block" : "none" ;
    if (_style == "none") {
	_style = "";
	s="-";
    }
    else {
	_style = "none";
	s="+";
    }
    id.style.display=_style;
    sign.innerHTML=s;
    //alert(id.innerTEXT);
    //alert(id.innerHTML);
}

/*  */
function collapseallByName(name)
{
    //var elements = document.getElementsByName(name);
    var elements = getElementsByName_iefix('div', name);
    for (var i = 0; i< elements.length;i++){ 
	collapseall(elements[i]);
    }
}

/*  */
function expandallByName(name)
{
    //var elements = document.getElementsByName(name);
    var elements = getElementsByName_iefix('div', name);
    for (var i = 0; i< elements.length;i++){ 
	expandall(elements[i]);
    }
}

/*  */
function selectAll(mytable, isInverse)  {
    var oTable = document.getElementById(mytable);
    var curobj = null;
    var row_count;
    var col_count;
    row_count = oTable.rows.length;
    col_count = oTable.rows[1].cells.length;
    for (var x = 0; x < row_count; x++) { 
	for (var y = 0; y < col_count; y++) {
	    if ( oTable.rows[x].cells[y].firstChild ) {
		curobj = oTable.rows[x].cells[y].firstChild;
		//      try{
		//		curobj.focus();
	        if ( curobj.nodeType == 1 ) {
		    //alert (curobj.nodeType);
		    if (curobj.type == 'checkbox' && !curobj.disabled) { 
    	    		if (isInverse) {
			    curobj.checked = (curobj.checked) ? false : true;
			} else {
			    curobj.checked = true; 
            		}
		    }
	        }
		//    } catch(e){
		//alert('Could not focus on myinput.');
		//}
	    }
	}
    }
}

/*  */
function selectAllcb(name,isInverse)  {
    var elements = document.getElementsByName(name);
    var tests = "";
    var testVisibility="";
    for (var i = 0; i< elements.length;i++){ 
	var elem =elements[i];
	testVisibility = elements[i].style.visibility;
	if ( elem.nodeType == 1 && !(testVisibility == "none")) {
	    alert("visibility:"+testVisibility); 
	    //alert (curobj.nodeType);
	    var myrow = document.getElementById(elem.value);
	    if (elem.type == 'checkbox' && !elem.disabled) { 
		if (isInverse) {
		    elem.checked = (elem.checked) ? false : true;
		} else {
		    elem.checked = true; 
		}
	    }
	}
    }
}

/*  */
function selectAllVisiblecb(name,isInverse)  {
    var elements = document.getElementsByName(name);
    for (var i = 0; i< elements.length;i++){ 
	var elem =elements[i];
	if ( elem.nodeType == 1 ) {
	    //alert (curobj.nodeType);
	    if (elem.type == 'checkbox' && !elem.disabled) { 
		var myrow = elem.parentNode.parentNode;
		//alert(myrow.type +"   id:"+myrow.id+"   n:"+myrow.name+ "  v:"+myrow.style.visibility+ "  d:"+myrow.style.display);
		if (myrow.style.display != 'none'||myrow.style.visibility == 'visible' || myrow.style.visibility == 'show'){
		    if (isInverse) {
			elem.checked = (elem.checked) ? false : true;
		    } else {
			elem.checked = true; 
		    }
		}
	    }
	}
    }
}


/*  */
function getElementsByName_iefix(tag, name) {
    var elem = document.getElementsByTagName(tag);
    var arr = new Array();
    for(i = 0,iarr = 0; i < elem.length; i++) {
	att = elem[i].getAttribute("name");
	if(att == name) {
	    arr[iarr] = elem[i];
	    iarr++;
	}
    }
    return arr;
}
		
/*function incdecimage(name,pixels,isincrease){
  var  max = 600;
  var  min = 100;
  if (isincrease){
  var newsize = imagesize + pixels;
  if (newsize < max)imagesize= newsize;
  else return; //false;
  }else{
  var newsize = imagesize - pixels;
  if (newsize > min) imagesize= newsize;
  else return;// false;
  }
  var elements = getElementsByName_iefix('IMG',name);
  for (var i = 0; i< elements.length;i++){ 
  var elem =elements[i];
  elem.height =imagesize;
  elem.width = imagesize;
  }
  //	var  max = 600;
  //	var  min = 100;
  //for (var i = 0; i< elements.length;i++){ 
  //		var elem =elements[i];
  //		if (elem == undefined || elem == null||elem.name != "profiles" ) return false;
  //	alert (elem.name );
  //		if (isincrease){
  ////	    	var newwidth = elem.width + pixels;
  //	    	var newheigth = elem.height + pixels; 
  //	    	if (newwidth < max) elem.width = newwidth;
  //	    	if (newheigth < max) elem.height = newheigth;
  //		}else{
  //	    	var newwidth = elem.width - pixels;
  //	    	var newheigth = elem.height - pixels; 
  //	    	if (newwidth > min) elem.width = newwidth;
  //	    	if (newheigth > min) elem.height = newheigth;
  //		}
  //	}
  }
*/

/*
function validation(){
    //Analysis Name
    if(document.getElementById("rpt_label").value==""){alert("Analysis name is required !");document.getElementById("rpt_label").focus();return;}
    //Mail
    else if(document.getElementById("mail").value==""){alert("Mail address is required !");document.getElementById("mail").focus();return;}
    else if(!document.getElementById("mail").value.match(/^[\S]+@[\S]+\.[\S]+$/)){alert("Invalid mail address, please correct it !");document.getElementById("mail").focus();return;}
    //Number of permutations
    else if(document.getElementById("nperm").value==""){alert("Number of permutation is required !");document.getElementById("nperm").focus();return;}
    else if(!document.getElementById("nperm").value.match(/^\d+$/)){alert("Number of permutation must be a positive integer !");document.getElementById("nperm").focus();return;}
    //phenotypes labels selection (at least 1)
    //max size
    else if(document.getElementById("set_max").value==""){alert("Max size is required !");document.getElementById("set_max").focus();return;}
    else if(!document.getElementById("set_max").value.match(/^\d+$/)){alert("Max size must be a positive integer !");document.getElementById("set_max").focus();return;}
    //min size
    else if(document.getElementById("set_min").value==""){alert("Min size is required !");document.getElementById("set_min").focus();return;}
    else if(!document.getElementById("set_min").value.match(/^\d+$/)){alert("Min size must be a positive integer !");document.getElementById("set_min").focus();return;}
    //max vs min
    else if(parseInt(document.getElementById("set_max").value)<parseInt(document.getElementById("set_min").value)){alert("Max size must be greather than min size !");return;}
    //number of markers
    else if(document.getElementById("num").value==""){alert("Number of markers is required !");document.getElementById("num").focus();return;}
    else if(!document.getElementById("num").value.match(/^\d+$/)){alert("Number of markers must be a positive integer !");document.getElementById("num").focus();return;}
    //number of graph sites of each phenotypes
    else if(document.getElementById("plot_top_x").value==""){alert("Number of graph sets of each phenotypes is required !");document.getElementById("plot_top_x").focus();return;}
    else if(!document.getElementById("plot_top_x").value.match(/^\d+$/)){alert("Number of graph sets of each phenotypes must be a positive integer !");document.getElementById("plot_top_x").focus();return;}
    //genrate all.
    //??
    else{
	var selectedphenotypes = document.actionForm.cls;//getElementsByName_iefix("input", "cls");
	var isOneSelected = false;
	for (var i = 0; i< selectedphenotypes.length;i++){ 
	    if(selectedphenotypes[i].checked == true) {
		isOneSelected = true;
		break;
	    }
	}
	if (isOneSelected == true)document.forms['actionForm'].submit();
	else {alert("You need to select at least one phenotypes comparison");document.getElementById("Pheno").focus();return;}
    }
}
*/