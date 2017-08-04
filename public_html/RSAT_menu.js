$(function(){
    title = document.title.replace("RSAT : ", "");
    title = title.replace("RSAT - ", "");
    
    menu_number = 0;
    genomes_genes = "Supported organisms;gene-info;gene-info result;infer operons;infer-operons result;get-orthologs;random gene selection;Random gene selection result";
    sequence_tools = "retrieve sequence;retrieve-seq result;retrieve EnsEMBL sequence;Retrieve sequences from genomic coordinates;retrieve-seq-bed result;fetch-sequences;purge-sequence;purge-sequence result;convert-seq;convert-seq result;random sequence;Random sequence result";
    motif_disco = "oligo-analysis;oligo-analysis result;oligo-diff;oligo-diff result;dyad-analysis;dyad-analysis result;pattern-assembly;pattern-assembly result;position-analysis;position-analysis result;local-word-analysis;local-word-analysis result;info-gibbs;info-gibbs result;consensus;Consensus result file";
    pattern_matching = "matrix-scan;matrix-scan result;matrix-scan QUICK and SIMPLE;crer-scan;crer-scan result;dna-pattern;dna-pattern result;genome-scale dna-pattern";
    drawing = "feature map;feature-map result;XYgraph;XYgraph result";
    web_services = "Web services;Web services documentation";
    help_contacts = "People;RSAT training material;RSA-tools - tutorials;RSA - Publications;";
    conversion = "Network Analysis Tools - convert-graph;NeAT - compare_classes;Frequency distribution (<i>classfreq</i>);classfreq result;convert-seq;convert-seq result;convert-matrix;convert-matrix result;create-background-model;create-background-model result;convert-background-model;convert-background-model result;seq-proba;seq-proba result;convert-features;convert-features result;compare-features;compare-features result";
    genetic_variations = "Variation information;retrieve variation sequence;variation-scan;Convert variation formats;";
    comparative_genomics = "get-orthologs;get-orthologs-compara;footprint-discovery;footprint-scan;";
    build_control_sets = "random gene selection;Random gene selection result;random sequence;Random sequence result;random genome fragments;Random genome fragments result;Random motif;random-motif result;permute-matrix;permute-matrix result;Random sites;random-sites result;Implant sites;implant-sites result";
    matrix_tools = "convert-matrix;convert-matrix result;compare-matrices;compare-matrices result;matrix-clustering;matrix-clustering result;matrix-distrib;matrix-distrib result;matrix-quality;matrix-quality result"; 
    ngs_chipseq = "peak-motifs;peak-motifs result;fetch-sequences;random genome fragments;Random genome fragments result";
    
    menu_numbers = [genomes_genes, sequence_tools, motif_disco, pattern_matching, drawing,web_services,help_contacts,conversion,genetic_variations, comparative_genomics, build_control_sets, matrix_tools, ngs_chipseq];
    
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


