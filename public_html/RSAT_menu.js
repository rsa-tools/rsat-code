function toggleMenu(menu_number) {
if (!document.getElementById) return;

var objID="menu"+menu_number;
var curr = document.getElementById(objID);
var ob = document.getElementById(objID).style;

var parentID="heading"+menu_number;
var parent = document.getElementById(parentID);


if (ob.display == 'block') { // mode "open", want it closed
	ob.display = 'none';
	parent.className ='menu_heading_closed';
}
else {
	ob.display = 'block';
	parent.className ='menu_heading_open';
}

}



function closeMenu(menu_number) {
	if (!document.getElementById) return;

	var objID="menu"+menu_number;
	var curr = document.getElementById(objID);
	var ob = document.getElementById(objID).style;

	var parentID="heading"+menu_number;
	var parent = document.getElementById(parentID);

	ob.display = 'none';
	parent.className ='menu_heading_closed';
}

function openMenu(menu_number) {
	if (!document.getElementById) return;

	var objID="menu"+menu_number;
	var curr = document.getElementById(objID);
	var ob = document.getElementById(objID).style;

	var parentID="heading"+menu_number;
	var parent = document.getElementById(parentID);

	ob.display = 'block';
	parent.className ='menu_heading_open';
}




/*Expand all the menu*/
function expandAll(last_menu_number) {
	var i=1;
	
	var objID="menu"+i;
	var curr = document.getElementById(objID);
	var ob = document.getElementById(objID).style;

	if (ob.display == 'block') { // mode "open", want it closed
		
		while (i<=last_menu_number)	{
			closeMenu(i);
			i=i+1;
		}
	}
	else {
		while (i<=last_menu_number)	{
			openMenu(i);
			i=i+1;
		}
	}
} 
