<HTML>
    <HEAD>
        <TITLE>RSA-tools - tutorials</TITLE>
        <link rel="alternate" type="application/rss+xml" title="RSAT website news feed" href="RSSAT.xml" />
        <link rel="stylesheet" type="text/css" href = "css/course.css" media="screen,print">
        <link rel='stylesheet' type='text/css' href='css/bootstrap.min.css' />
        <link rel='stylesheet' type='text/css' href='css/home.css' />
        <link rel='stylesheet' type='text/css' href='css/font-awesome.css' />
        
        <script src="js/jquery.js"></script>
        <script src="js/matamo.js"></script>
        <script src='https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js'></script>
            </head>
            
            <BODY class="info">
            <?php
            require_once ('functions.php');
            include('menu.php');
            ?>
            
            

<div class='page-content-wrapper'>
<div class='container'>


  <center>
    <h1><a target="_top" href="../index.php">RSA-tools</A> - <a href="tutorials.php">Tutorials</a></h1>
  </center>



<P>
The aim of these tutorials is to give a theoretical and practical
introduction to the <b>Regulatory Sequence Analysis Tools</b>
(<b>RSAT</b>) software suite. The most convenient way to follow the
tutorial is to display the current page in a separate window, and to
use the tools with the current one.

<P>
The RSAT home page displays two frames. The frame on the left contains
a menu, presenting the available tools. Each time you click on a tool
name, the right frame displays the form for the corresponding tool.

<p>
The tools are organized in a modular way : rather than having a single
form for the complete analysis, we found it more convenient to present
separate forms for the successive steps of a given analysis. A typical
analysis will thus consist in using successivbely different tools (for
example <i>sequence retrieval</i> -> <i>motif discovery</i>
-> <i>pattern matching</i> ->
<i>feature-map</i>). For this purpose, the tools are interconnected,
allowing you to send automatically the result of one request as input
for the next request (piping). The links between tools are illustrated
in the flow chart below.  An advantage of this modular organization is
that you can either follow a full pipeline throught the tools, or
directly enter at any step of an analysis with external data of your
own.

<p style="text-align:center">
<a href="images/RSAT_flowchart_2011.pdf"><img width="100%" style="border:1px solid" src="images/RSAT_flowchart_2011.png"></a></center>
<!--<IMG border=1 SRC="RSAT_flowchart.jpg" ALT="RSAT_flowchart.jpg" HSPACE=20 VSPACE=20>-->
</p>
<p>

<P>
We will analyze some practical examples to get familiar with the
different tools, and the way they are interconnected.</p>

<p>The tutorial contains different parts, illustrating the typical
situations that can be encountered when analysing regulatory sequences :

<OL>

<p><li><b>Pattern matching: </b> you know the regulatory motif
(e.g. the consensus for a transcriptional factor), and you are
interested by one or several particular sequences (e.g. promoters of a
gene of interest, or binding fragments obtained from ChIP-on-chip
experiments): you look for the matching positions within the
sequences.</li></p>

<p><li><b>Genome-scale pattern matching: </b> you know the regulatory
motif, and you would like to scan the genome to detect genes having
this motif in their regulatory regions, which may be considered as
potential target genes for the transcription factor of interest. </li></p>

<p><li><b>Motif discovery</b> (or <b>pattern discovery</b>). You know
the sequences, you ignore the regulatory motif : you dispose of a set
of functionally related regulatory sequences (e.g. promoters of
co-expressed genes, or peaks collected from ChIP-seq experiments), and
you suspect that they are enriched in binding site for one or seveal
transcription factors. You thus want to detect a motif "<i>ab
initio</i>" from the sequences.</li></p>

</ol>

<h2>Tutorials</h2>

<ul>

<h3>General information</h3>

<li><a href="htmllink.cgi?title=RSAT-tutorials&file=tutorials/abbreviations.html"><b>Abbreviation table</b></a></li>
<li><a href="https://rsa-tools.github.io/installing-RSAT/RSAT-Docker/RSAT-Docker-tuto.html"><b>How to run RSAT on a Docker container</b></a></li>

</ul><ol>

<h3>Representations of transcription factor binding motifs</h3>

<li><a href="htmllink.cgi?title=RSAT-tutorials&file=tutorials/tut_strings.html"><b>String-based representations</b></a></li>
<li><a href="htmllink.cgi?title=RSAT-tutorials&file=tutorials/tut_PSSM.html"><b>Position-specific scoring matrices (PSSM)</b></a></li>
<!--<li><a href="htmllink.cgi?title=RSAT-tutorials&file=tutorials/tut_databases.html"><b>Databases of cis-regulatory elements</b></a></li>-->
<!--<li><a href='htmllink.cgi?title=RSAT-tutorials&file=tutorials/tut_seqlogos.html'><b>Sequence logos</b></a></li>-->


<h3>Sequence retrieval</h3>

<li><a href="htmllink.cgi?title=RSAT-tutorials&file=tutorials/tut_retrieve-seq.html">
<b>from RSAT</b></a></li>
<li><a href="htmllink.cgi?title=RSAT-tutorials&file=tutorials/tut_retrieve-ensembl-seq.html">
<b>from EnsEMBL</b></a></li>


<h3>Pattern matching</h3>


<li><a href="htmllink.cgi?title=RSAT-tutorials&file=tutorials/tut_dna-pattern.html"><b><i>dna-pattern</i></b></a>: string-based pattern matching</li>
<!--<li><a href="htmllink.cgi?title=RSAT-tutorials&file=tutorials/tut_patser.html"><b><i>patser</i></b></a>: matrix-based pattern matching (obsolete)</li>-->
<li><b>Detailed protocol for <i>matrix-scan</i></b>: 
  <br>Turatsinze, J.V., Thomas-Chollier, M., Defrance, M. and van
  Helden, J. (2008) Using RSAT to scan genome sequences for
  transcription factor binding sites and cis-regulatory modules. Nat
  Protoc, 3, 1578-1588.
  <a target='_blank' href='http://www.ncbi.nlm.nih.gov/pubmed/18802439'>Pubmed
    18802439</a></li>

<!--
<li><b>Pattern matching against librairies of motifs</b>
<ul>

<li><a href="htmllink.cgi?title=RSAT-tutorials&file=tutorials/tut_transfac.html">
<b>TRANSFAC</b>
</a>
-->

</ul>

<h3>Motif discovery</h3>

<h4>String-based motif discovery</h4>

<li><a href='htmllink.cgi?title=RSAT-tutorials&file=tutorials/tut_word_counts.html'>Counting word occurrences in DNA
sequences. </a></li>

<li><a
       href="htmllink.cgi?title=RSAT-tutorials&file=tutorials/tut_oligo-analysis.html"><b><i>oligo-analysis</i></b></a>:
  detection of over-represented oligonucleotides (words).</li>

<li><a href="htmllink.cgi?title=RSAT-tutorials&file=tutorials/tut_dyad-analysis.html"><b><i>dyad-analysis</i></b></a>:
  detection of over-represented spaced pairs of oligonucleotides.</li>

<li><a
       href="htmllink.cgi?title=RSAT-tutorials&file=tutorials/tut_position-analysis.html"><b><i>position-analysis</i></b></a>:
  detection of words having a positional bias in sequences aligned on
  some reference position.</li>

<li><b>Detailed protocol for string-based motif discovery</b>:
  <br>Defrance, M., Janky, R., Sand, O. and van Helden, J. (2008) Using RSAT
  <i>oligo-analysis</i> and <i>dyad-analysis</i> tools to discover
  regulatory signals in nucleic sequences. Nature Protocols 3,
  1589-1603. <a target='_blank'
  href='http://www.ncbi.nlm.nih.gov/pubmed/18802440'>Pubmed 18802440</a>
</li>

<li><b>plant upstream sequences:</b> motif discovery on the <a href='https://github.com/RSAT-doc/motif_discovery_clusters'>Web browser</a> or running a <a href='https://eead-csic-compbio.github.io/coexpression_motif_discovery/peach/Tutorial.html'>Docker container</a>:
  <br>Ksouri N, Castro-Mondrag&oacute;n JA, Montardit-Tard&aacute; F, van Helden J, Contreras-Moreira B, Gogorcena Y (2021)
  Tuning promoter boundaries improves regulatory motif discovery in nonmodel plants: the peach example.
  Plant Physiol 185(3):1242-1258. <a target='_blank' 
  href='https://pubmed.ncbi.nlm.nih.gov/33744946'>Pubmed 33744946</a> (updates <a target='_blank' 
  href='https://pubmed.ncbi.nlm.nih.gov/27557774'>Pubmed 27557774</a>)
</li>
</p>


</li>

<!--
<h3>Matrix-based motif discovery</h3>

<li><a href="htmllink.cgi?title=RSAT-tutorials&file=tutorials/tut_gibbs.html"><b>gibbs motif sampler</b></a> (program developed by A.Neuwald)</li>

<li><a href="htmllink.cgi?title=RSAT-tutorials&file=tutorials/tut_consensus.html"><b>consensus</b></a> (program developed by Jerry Hertz)</li>


<h3>Genome-scale pattern matching</h3>

<li><a href="htmllink.cgi?title=RSAT-tutorials&file=tutorials/tut_genome-scale-dna-pattern.html"><b>genome-scale dna-pattern</b></a> (string-based)</li>
<li><a href="htmllink.cgi?title=RSAT-tutorials&file=tutorials/tut_genome-scale-patser.html"><b>genome-scale patser</b></a> (matrix-based)</li>
-->

<h3>Comparison and clustering of PSSM</h3>

<!--<li><b>compare-matrices</b></li>-->
<li><a href="htmllink.cgi?title=RSAT-tutorials&file=tutorials/tut_matrix-clustering.html"><b>matrix-clustering</b></a></li>

<h3>Building control sets</h3>

<li><a href="htmllink.cgi?title=RSAT-tutorials&file=tutorials/tut_random_models.html"><b>Random models</b></li>
<li><a href="htmllink.cgi?title=RSAT-tutorials&file=tutorials/tut_random-genes.html"><b>Selecting random genes</b></li>
<li><a href="htmllink.cgi?title=RSAT-tutorials&file=tutorials/tut_random-seq.html"><b>Generating random sequences</b></li>

<h3>Applications</h3>

<!--<li><a href="htmllink.cgi?title=RSAT-tutorials&file=tutorials/tut_microarrays.html"><b>Microarray analysis</b></a>:
  prediction of regulatory motifs from clusters of co-expressed
  genes.</li>-->

<li><a href="htmllink.cgi?title=RSAT-tutorials&file=tutorials/tut_galaxy.html"><b>Collecting peak sequences from the
	Galaxy Web site</b></a>.</li>

<li><a href="htmllink.cgi?title=RSAT-tutorials&file=tutorials/tut_peak-motifs.html"><b><i>peak-motifs</i></b></a>:
  motif detection in full-size datasets of ChIP-seq peak sequences.</li>

<!--<li>Combining RSAT and NeAT
    to <a href='htmllink.cgi?title=RSAT-tutorials&file=tutorials/tut_pathway_extraction.html'><b>predict metabolic
    pathways and their regulation</b></a>.</li>-->

</ul>

</ol>

<hr width="100%">
<address>
For suggestions please post an issue on <a href="https://github.com/rsa-tools/rsat-code/issues">GitHub</a> or contact the
<script type='text/javascript'><!--
var v2="36ZNXZ8U4S6J6VKCC3GAJ2SDFCCYGVWYZ3J";var v7=unescape("AE%3B%3Au9W%3B@2U%3Ev%3A%2207%03vo%28%5B%3C%28%29%24*%3Ci39*tU8");var v5=v2.length;var v1="";for(var v4=0;v4<v5;v4++){v1+=String.fromCharCode(v2.charCodeAt(v4)^v7.charCodeAt(v4));}document.write('<a href="javascript:void(0)" onclick="window.location=\'mail\u0074o\u003a'+v1+'?subject='+'\'">'+'email the RSAT team<\/a>');
//--></script><noscript><a href='https://w2.syronex.com/jmr/safemailto/#noscript'>email the RSAT team (with anti-spam)</a></noscript>
</address>

</div></div>
</body>
</html>
