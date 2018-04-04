menu_number = 0;
    heading1 = "Supported organisms;gene-info;gene-info result;infer operons;infer-operons result;get-orthologs;random gene selection;Random gene selection result";
    heading2 = "retrieve sequence;retrieve-seq result;retrieve EnsEMBL sequence;Retrieve sequences from genomic coordinates;retrieve-seq-bed result;fetch-sequences;purge-sequence;purge-sequence result;convert-seq;convert-seq result;random sequence;Random sequence result";
    heading3 = "oligo-analysis;oligo-analysis result;oligo-diff;oligo-diff result;dyad-analysis;dyad-analysis result;pattern-assembly;pattern-assembly result;position-analysis;position-analysis result;local-word-analysis;local-word-analysis result;info-gibbs;info-gibbs result;consensus;Consensus result file";
    heading4 = "matrix-scan;matrix-scan result;matrix-scan QUICK and SIMPLE;crer-scan;crer-scan result;dna-pattern;dna-pattern result;genome-scale dna-pattern";
    heading5 = "feature map;feature-map result;XYgraph;XYgraph result";
    heading6 = "Web services;Web services documentation";
    heading7 = "People;RSAT training material;RSA-tools - tutorials;RSA - Publications;";
    heading8 = "Network Analysis Tools - convert-graph;NeAT - compare_classes;Frequency distribution (<i>classfreq</i>);classfreq result;convert-seq;convert-seq result;convert-matrix;convert-matrix result;create-background-model;create-background-model result;convert-background-model;convert-background-model result;seq-proba;seq-proba result;convert-features;convert-features result;compare-features;compare-features result";
    heading9 = "Variation information;variation-info result;retrieve variation sequence;retrieve-variation-seq result;variation-scan;variation-scan result;Convert variation formats;convert-variations result";
    heading10 = "get-orthologs;get-orthologs result;get-orthologs-compara;get-orthologs-compara result;footprint-discovery; footprint-discovery result;footprint-scan;footprint-scan result;";
    heading11 = "random gene selection;Random gene selection result;random sequence;Random sequence result;random genome fragments;Random genome fragments result;Random motif;random-motif result;permute-matrix;permute-matrix result;Random sites;random-sites result;Implant sites;implant-sites result";
    heading12 = "retrieve-matrix;convert-matrix;convert-matrix result;compare-matrices;compare-matrices result;matrix-clustering;matrix-clustering result;matrix-distrib;matrix-distrib result;matrix-quality;matrix-quality result"; 
    heading13 = "peak-motifs;peak-motifs result;fetch-sequences;random genome fragments;Random genome fragments result";
    
    menu_numbers = [heading1, heading2, heading3, heading4, heading5,heading6,heading7,heading8,heading9, heading10, heading11, heading12, heading13];

$(function(){
    title = document.title.replace("RSAT : ", "");
    title = title.replace("RSAT - ", "");
    title = title.replace("NeAT - ", "");   
    $.each(menu_numbers, function(key, value){
        if(value.indexOf(title) != -1){
            menu_number = key + 1;
            var objID="menu"+menu_number;
            var parentID="heading"+menu_number;
            var parent = document.getElementById(parentID);
    
            document.getElementById(objID).style.display = 'block';
            parent.className ='menu_heading_open';
        }
    });
    
    // for NeAT
    title = title.replace("Network Analysis Tools - ", "");
    title = title.replace("Network Analysis Tools : ", "");
    title = title.replace("NeAT - ", "");
    menu_graph_number = 0;
    network = "compare-graphs;graph-topology;graph-neighbours;random-graph;alter-graph;convert-graph;display-graph";
    pathfinding = "Pathfinder;pathway-extractor;";
    cluster = "compare_classes;contingency-stats;convert-classes;";
    cluster_network = "MCL;RNSC;graph-get-clusters;graph-get-clusters;";
    graphic = "draw-heatmap";
    roc = "roc-stats";
    data = "download string dataset;KEGG network provider;";
    neat_web_services = "Web Services;Web Services - Taverna workflows";
    menu_graph_numbers = [network, pathfinding, cluster, cluster_network, graphic, roc, data, neat_web_services];
    
    $.each(menu_graph_numbers, function(key, value){
        if(value.indexOf(title) != -1){
            menu_graph_number = key + 1;
            var objID="menu"+menu_graph_number;
            var parentID="heading"+menu_graph_number;
            var parent = document.getElementById(parentID);
    
            document.getElementById(objID).style.display = 'block';
            parent.className ='menu_heading_open';
        }
    });
 
});

////// For search tool bar
function searchfunc(){
	var term = $("#searchfunc").val();
	closeAll();
	var menuitems = document.querySelectorAll(".menu_item, .menu_item_last ");
	for(var x = 0; x < menuitems.length; x++){
		t = (menuitems[x].textContent === undefined) ? menuitems[x].innerText : menuitems[x].textContent;
		html = menuitems[x].innerHTML;	
		menuitems[x].innerHTML = t;
	
		index = html.indexOf('<img');
		if(index != -1){
			menuitems[x].innerHTML += html.substring(index);			
		}			
	}
	var open = false;
	if(term.length > 0){
		var allmenu = document.getElementsByClassName('menu_heading_closed');
		for(var x = 1; x <= menu_numbers.length; x++){
			var obj = document.getElementById("menu"+x);
			k = obj.childNodes;
			for(i = 0; i < k.length; i++){
				if(k[i].className == 'menu_item' || k[i].className == 'menu_item_last'){
					html = k[i].innerHTML;
					found = (k[i].textContent === undefined) ? k[i].innerText : k[i].textContent;				
					index = found.toLowerCase().indexOf(term.toLowerCase());
                    
					if(index != -1){
						k[i].innerHTML = html.substring(0,index) + "<span style='background-color:yellow;color:black'>" + html.substring(index,index+term.length) + "</span>" + html.substring(index+term.length);
						open = true;
					}
				}
			}
			if(open){ openMenu(x); open=false; }							
		}
	}
}


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
/*Close all menu*/
function closeAll(){
	for(var i = 1; i< 14; i++){
		closeMenu(i);	
	}	
}

