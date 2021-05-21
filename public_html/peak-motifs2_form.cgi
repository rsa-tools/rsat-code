#!/usr/bin/env perl
#### this cgi script fills the HTML form for the program matrix-clustering
BEGIN {
    if ($0 =~ /([^(\/)]+)$/) {
	push (@INC, "$`lib/");
    }
    require "RSA.lib";
}
#if ($0 =~ /([^(\/)]+)$/) {
#    push (@INC, "$`lib/");
#}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";
use RSAT::matrix;
use RSAT::MatrixReader;

$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

################################################################
## Default values for filling the form
$default{demo_descr} = "";
$default{lth_occ_sig}=0;
$default{uth_pval} = "1e-4";
$default{assembly} = "";
$default{oligo_length6}="checked='checked'";
$default{oligo_length7}="checked='checked'";
$default{oligo_length8}="";
$default{merge_lengths}="";
$default{"oligo-analysis"}="checked";
$default{"dyad-analysis"}="";
$default{"position-analysis"}="checked";
$default{'local-word-analysis'}="";
$default{'local-word-analysis_dyads'} ="";
$default{"position-analysis_dyads"} ="";
$default{'matrix-scan-quick'}="checked";
$default{compare_motif_db}="checked";
$default{title}="";
$default{"seq_format1"}="none";
$default{"seq_format2"}="none";
$default{"dust_seq1"}="";
$default{"dust_seq2"}="";
$default{max_seq_len}="500";
$default{top_sequences}="";
$default{nmotifs} = 5;
$default{disc_markov}="auto";
$default{scan_markov}="1";
$default{origin} = "center";
$default{offset} = "0";
## Plot XY graph with R rather than GD
if ($ENV{R_supported}) {
    $default{r_plot} = "checked";
} else {
    $default{r_plot} = "";
}
$default{visualize}="none";


## Motif database
#$default{compare_motif_database}="jaspar_core_nonredundant_vertebrates";
#$default{custom_motif_db_name}="custom_motif_collection";
$default{compare_motif_database} = "Jaspar";
$default{compare_motif_collection}="jaspar_core_nonredundant_vertebrates"; ## I (JvH) SHOULD ADAPT THIS IN ORDER TO PRESENT DIFFERENT DATABASES DEPENDING ON TAXON SPECIFICITY OF THE SERVER (2015-03-13)

### Replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
   if ($query->param($key) =~ /checked/i) {
    $checked{$key} = "CHECKED";
  }
  if ($key eq "visualize"){
  	$checked{$query->param($key)} = "CHECKED";
  }
}

$checked{$default{visualize}} = "checked";
&ListParameters() if ($ENV{rsat_echo} >= 2);

################################################################
### print the form ###

################################################################
### header
&RSA_header_bootstrap("peak-motifs2", "form");

print $query->start_multipart_form(-action=>"peak-motifs2.cgi");

print '
 <!-- Form with bootstrap -->
<div class="container">
	<div class="row">
        <div class="col-lg-9 col-md-5 col-sm-8 col-xs-9 bhoechie-tab-container">
            <div class="col-lg-2 col-md-3 col-sm-3 col-xs-3 bhoechie-tab-menu">
              <div class="list-group">
                <a href="#" class="list-group-item active text-center">
                  <h4 class="glyphicon"><i class="fa fa-info-circle fa-2x"></i></h4><br/>peak-motifs2
                </a>
                <a href="#" class="list-group-item text-center">
                  <h4 class="glyphicon"><i class="fa fa-tag fa-2x"></i></h4><br/>Mandatory inputs
                </a>
                <a href="#" class="list-group-item text-center">
                <h4 class="glyphicon"><i class="fa fa-tasks fa-2x"></i></h4><br/>Mandatory options
                </a>
                <a href="#" class="list-group-item text-center">
                  <h4 class="glyphicon"><i class="fa fa-tags fa-2x"></i></h4><br/>Optional inputs
                </a>
                <a href="#" class="list-group-item text-center">
                  <h4 class="glyphicon"><i class="fa fa-tasks fa-2x"></i></h4><br/>Advanced options
                </a>
                <a href="#" class="list-group-item text-center">
                  <h4 class="glyphicon"><i class="fa fa-play-circle fa-2x"></i></h4><br/>Run analysis
                </a>
              </div>
            </div>
            <div class="col-lg-9 col-md-9 col-sm-9 col-xs-9 bhoechie-tab">


 <!-- ################################################################ -->
 <!-- ### info tab ### -->
            <div class="bhoechie-tab-content active">
              <h2>
                <img src="images/RSAT_logo.jpg" style="max-width:150px;max-height:60px;padding-bottom:10px" alt="RSAT server" border="0"></img>
                peak-motifs2
              </h2>
              <span class="fa-stack fa-lg" style="color: red;">
      							<i class="fa fa-exclamation-triangle fa-stack-1x"></i>
    					</span>
              This is the peak-motifs2 development version.
              <br>&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;
              <a href="peak-motifs_form.cgi" target="_blank"> You can use the previous version of <i>peak-motifs</i> here </a><br>
              <span class="fa-stack fa-lg">
  						<i class="fa fa-info-circle fa-stack-1x"></i>
					    </span>
                Discover exceptional motifs (over-represented, positionally biased) in a collection<br>
                &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; of ChIP-seq peaks. <br>
              <span class="fa-stack fa-lg">
 							  <i class="fa fa-user fa-stack-1x"></i>
					    </span>

					    <a target="_blank" href="http://morgane.bardiaux.fr/">Morgane Thomas-Chollier</a>,
              <a target="_blank" href="http://jacques.van-helden.perso.luminy.univ-amu.fr/ ">Jacques van Helden</a>,
              <a target="_blank" href="https://www.researchgate.net/profile/Matthieu_Defrance">Matthieu Defrance</a>,
              <a target="_blank" href="http://www-good.ibl.fr/en/bio-informatique/">Olivier Sand</a>,
              <br> &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;
              <a target="_blank" href="https://www.hdsu.org/">Carl Herrmann</a>,
              <a target="_blank" href="https://www.ibens.ens.fr/spip.php?rubrique27&lang=en">Denis Thieffry</a>
              <br>

					<span class="fa-stack fa-lg">
  							<i class="fa fa-folder-open fa-stack-1x"></i>
					</span>

          <a target="_blank" href="sample_outputs/peak-motifs_demo_single_20210216/peak-motifs_synthesis.html"> Sample output single dataset </a> &nbsp;
          <span class="fa-stack fa-lg">
                <i class="fa fa-folder-open fa-stack-1x"></i>
          </span>

          <a target="_blank" href="sample_outputs/peak-motifs_demo_test_ctrl_20210216/peak-motifs_synthesis.html"> Sample output test and ctrl dataset </a> <br>
					<span class="fa-stack fa-lg">
  							<i class="fa fa-book fa-stack-1x"></i>
					</span>

					<a class="iframe" href="help.peak-motifs.html">User Manual</a> &nbsp;
					<span class="fa-stack fa-lg">
  							<i class="fa fa-graduation-cap fa-stack-1x"></i>
					</span>
					<a class="iframe" href="tutorials/tut_peak-motifs/tut_peak-motifs_new.html">Tutorial</a><br>

					<span class="fa-stack fa-lg">
  							<i class="fa fa-twitter fa-stack-1x"></i>
					</span>
					<a href="https://twitter.com/rsatools" target="_blank">Ask a question to the RSAT team</a><br>

					<span class="fa-stack fa-lg">
  							<i class="fa fa-pencil fa-stack-1x"></i>
					</span>
					Cite the publication: <a href="https://twitter.com/rsatools" target="_blank"></a><br>
					<div class="panel panel-default">
  						<div class="panel-body">
               Thomas-Chollier, M., Herrmann, C., Defrance, M., Sand, O., Thieffry, D. and van Helden, J. (2011). <i>"RSAT peak-motifs: motif analysis in full-size ChIP-seq datasets"</i>. Nucleic Acids Research, 40(4): e31. <a href="https://pubmed.ncbi.nlm.nih.gov/22156162/" target="_blank">[Pubmed]</a><a href="https://academic.oup.com/nar/article/40/4/e31/2411061?keytype=ref&ijkey=zOvloLjtKzL73F8" target="_blank">[Full text]</a>
               <br><br>
               Thomas-Chollier M, Darbo E, Herrmann C, Defrance M, Thieffry D, van Helden J. (2012). <i>"A complete workflow for the analysis of full-size ChIP-seq (and similar) data sets using peak-motifs"</i>. Nat Protoc, 7(8): 1551-1568. <a href="https://pubmed.ncbi.nlm.nih.gov/22836136/" target="_blank"> [Pubmed]</a> <a href="https://www.nature.com/articles/nprot.2012.088" target="_blank">[Full text]</a>
              </div>
					</div>
          &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;

          <br>
        </div>

 <!-- ################################################################ -->
 <!-- ### mandatory inputs ### -->
       <div class="bhoechie-tab-content">
       <!-- title -->
        <div class="panel panel-danger">
          <div class="panel-heading">
            Analysis Title
            <i class="fa fa-info-circle" data-container="body" data-toggle="tooltip" data-placement="top" title="Title that will be displayed at the top of the report page." data-original-title=""></i>
          </div>
          <div class="panel-body">
            <div class="form-group">';
print $query->textfield(-id=>'title',
                       -name=>'title',
                       -class=>'form-control',
                       -placeholder=>'Provide a Title for this analysis ',
                       -size=>40,
                       -required=>'true',
			                 -default=>$default{title});
print '
      		  </div>
			    </div>
			 </div>
       <!-- Peak sequences -->
       <div class="panel panel-danger">
        <div class="panel-heading">Peak sequences
          <i class="fa fa-info-circle" data-container="body" data-toggle="tooltip" data-placement="top" title="Input here the set of peak sequences of interest where motifs will be discovered, you can either paste the peak sequences in the text box or upload it from your computer" data-original-title=""></i>
        </div>
			  <div class="panel-body">';
        &MultiSeqFormatInput_bootstrap('seq_title'   =>'Peak sequences',
                                        'seq_num'    =>1,
                                        'seq_type'   =>'fasta sequences',
                                        'show_genome'=>1);
print ' </div>
       </div>
      </div>';
print  '<!-- ################################################################ -->
      <!-- ### mandatory options ### -->
          <div class="bhoechie-tab-content">

          <!-- title -->
           <div class="panel panel-danger">
             <div class="panel-heading">
               Motif Discovery
               <i class="fa fa-info-circle" data-container="body" data-toggle="tooltip" data-placement="top" title="Set the options that will be used for motif discovery in your set of sequences." data-original-title=""></i>
             </div>
             <div class="panel-body">';
print "<div class='form-group'>";
print "<p> <b>Oligonucleotides (<i>k</i>-mers) and spaced word pairs (dyads):</b>";
# oligo-analysis
print         '<div class="form-check">';
                print $query->checkbox(-name=>'oligo-analysis',
                                       -id=>'oligo-analysis',
                                       -class=>'form-check-input',
                		                   -checked=>$default{"oligo-analysis"},
                		                   -label=>'');
print '         <label class="form-check-label" for="oligo-analysis">
                    Discover over-represented words <a class="iframe" href="help.oligo-analysis.html"> [oligo-analysis]</a>
                </label>
                </div>';
# position-analysis
print         '<div class="form-check">';
             print $query->checkbox(-name=>'position-analysis',
                                    -id=>'position-analysis',
                                    -class=>'form-check-input',
             		                    -checked=>$default{"position-analysis"},
             		                    -label=>'');
print '         <label class="form-check-label" for="position-analysis">
                    Discover words with a positional bias <a class="iframe" href="help.position-analysis.html">[position-analysis]</a>
                </label>
                </div>';
# local-word-analysis
print         '<div class="form-check">';
              print $query->checkbox(-name=>'local-word-analysis',
                                     -id=>'local-word-analysis',
                                     -class=>'form-check-input',
                                		 -checked=>$default{"local-word-analysis"},
                                		 -label=>'');
print '         <label class="form-check-label" for="local-word-analysis">
                    Discover words with local over-representation <a class="iframe" href="help.local-word-analysis.html">[local-word-analysis]</a>
                </label>
                </div>';
# local-word-analysis
print         '<div class="form-check">';
              print $query->checkbox(-name=>'dyad-analysis',
                                     -id=>'dyad-analysis',
                                     -class=>'form-check-input',
                                		 -checked=>$default{"dyad-analysis"},
                                		 -label=>'');
print '         <label class="form-check-label" for="dyad-analysis">
                    Discover over-represented spaced word pairs </b><a class="iframe" href="help.dyad-analysis.html">[dyad-analysis] </a>
                </label>
                </div>';
print "<i>Note: position-analysis and local-word-analysis will not run if a control set is provided.</i>
       </div>";

## Word size

print   "<div class='form-group'>
            <p> <b>Oligonucleotides lengths:</b></p>
          <div class='form-group'>
            <label>
              <input type='checkbox' name='oligo_length6' id='oligo_length6' value='on' ".$default{oligo_length6}."'>6</label>
            <label>
              <input type='checkbox' name='oligo_length7' id='oligo_length7' value='on' ".$default{oligo_length7}."'>7</label>
            <label>
              <input type='checkbox' name='oligo_length8' id='oligo_length8' value='on' ".$default{oligo_length8}."'>8</label>
              &nbsp;&nbsp;&nbsp;<a class='badge badge-primary iframe' href='help.peak-motifs.html#oligo_len'>Info</a>
          </div>
          <div class='form-group'>
            <label>
              <input type='checkbox' name='merge_lengths' id='merge_lengths' value='on' ".$default{merge_lengths}."'>
              merge lengths for assembly
              &nbsp;&nbsp;&nbsp;<a class='badge badge-primary iframe' href='help.peak-motifs.html#oligo_len'>Info</a>
            </label>
            <br>
            <i>Note: motifs can be larger than word sizes (words are used as seed for building matrices).</i>
          </div>
        </div>";

## Markov  order (for oligo-analysis in single strand mode)
print "<div class='form-group'>
        <p> <b>Markov order (m) of the background model:</b>";
        print '<div class="form-check">';
        ## Markov order for scanning (site prediction + motif enrichment)
        print $query->popup_menu(-id=>'markov',
                                -name=>'markov',
                                -Values=>["auto","0", "1", "-3", "-2"],
                                -labels=>{"auto"=> "automatic (adapted to sequence length)",
                                          "0"=>"m=0 (generally not ideal)",
                                          "1"=>"m=1 (more sensitive for small data sets, e.g. 100kb)",
                                          "-3"=>"m=k-3 (intermediate size data sets)",
                                          "-2"=>"m=k-2 (more stringent for large data sets e.g. > 1Mb)"},
                                -class=>'form-control',
                                -default=>$default{disc_markov});
        print '</div>';
        print "<i>Note: only for oligo-analysis in single datasets (it will be ignored if control set is provided).</i>";
print "</div>";
## Number of motifs per algorithm
print "<div class='form-group'>";
print ' <label class="form-check-label" for="nmotifs"><b>Number of motifs per algorithm</b>:</label>';
print $query->popup_menu(-name=>'nmotifs',
                         -id=>'nmotifs',
                         -class=>'form-control',
			                   -Values=>[1,2,3,4,5,6,7,8,9,10],
			                   -default=>$default{nmotifs});
print "</div>";

### 2str or 1str
print "<div class='form-group'>";
print ' <label class="form-check-label" for="strand"><b>Search on</b>:</label>';
print $query->popup_menu(-name=>'strand',
                         -id=>'strand',
			                   -Values=>["-2str","-1str"],
                         -labels=>{"-2str"=>"both strands","-1str"=>"single strand"},
                         -class=>'form-control',
			                   -default=>'-2str');
print "</div>";

   print '
              </div>
            </div>
         </div>';

 print '<!-- ################################################################ -->
 <!-- ### optional inputs ### -->
     <div class="bhoechie-tab-content">

      <div id="accordion" role="tablist">
        <!-- Control set of sequences-->
        <div class="card">
        <div class="card-header" role="tab" id="headingOne">
        <h5> <i class="fa fa-tags"></i>
        <a data-toggle="collapse" href="#collapseOne" aria-expanded="true" aria-controls="collapseOne">
          Add a set of peak sequences to contrast during motif discovery
        </a>
      </h5>
    </div>

    <div id="collapseOne" class="collapse" role="tabpanel" aria-labelledby="headingOne" data-parent="#accordion">
      <div class="card-body">
        <div class="panel panel-warning">
          <div class="panel-heading">Contrasting sequences
            <i class="fa fa-info-circle" data-container="body" data-toggle="tooltip" data-placement="right" title="Input here an additional set of peak sequences to be used as contrast for motif discovery, you can paste the peak sequences in the text box or upload them from your computer" data-original-title=""></i>
          </div>
			<div class="panel-body">';
      &MultiSeqFormatInput_bootstrap('seq_title'   =>'Peak sequences',
                                      'seq_num'    =>2,
                                      'seq_type'   =>'fasta sequences',
                                      'show_genome'=>0);

print '</div>
      </div>
     </div>
    </div>
  </div>
  <!-- Add reference motifs-->
  <div class="card">
    <div class="card-header" role="tab" id="headingTwo">
      <h5 class="mb-0"> <i class="fa fa-tags"></i>
      <a class="collapsed" data-toggle="collapse" href="#collapseTwo" aria-expanded="false" aria-controls="collapseTwo">
       Add reference motif collections to compare to the newly discovered motifs
      </a>
      </h5>
    </div>
    <div id="collapseTwo" class="collapse" role="tabpanel" aria-labelledby="headingTwo" data-parent="#accordion">
      <div class="card-body">
        <div class="panel panel-warning">
          <div class="panel-heading">Reference motif collections
            <i class="fa fa-info-circle" data-container="body" data-toggle="tooltip" data-placement="top" title="Select a motif collection from a database or provide a personal  reference motif collection to compare discovered motifs." data-original-title=""></i>
          </div>
			    <div class="panel-body">';
print "     <b>Choose an available collection in public databases:</b><br><br>";
&MotifSelection_bootstrap("mode" => "checkbox"); #NOTE WSG. This is an extremely slow function
print "     <br><br>";

## Personal motif collection
print '<div class="form-inline">
        <p><b>Add your own private motif collection:</b></p>
        <div class="form-group">';
print  $query->textfield(-name=>'custom_motif_db_name',
                         -class=>'form-control',
                         -default=>$default{custom_motif_db_name},
                         -placeholder=>'Provide a name for your Motif Collection',
                         -size=>40);
print '</div>
       <div class="form-group">';
print $query->filefield(-name=>'custom_motif_db',
                        -class=>"form-control-file",
                        -size=>10);
print "</div>";
print '</div>';
#print '</div>';
#print "<br>Matrices should be in <b>Transfac format</b> (other formats can be converted with <a href='convert-matrix_form.cgi'><i>convert-matrix</i></a>).";

## Reference motifs
#print"</p>";
#print "<br><br>";
#print '<div class="form-group">';

#print "<label><b>Add known reference motifs for this experiment:</b><label>";
#print $query->filefield(-name=>'ref_motif',
#                        -size=>10);
#print "</div>";
print "<br><i>Note: your personal motif collection should be in <b>Transfac format</b>";
print " (other formats can be converted with <a href='convert-matrix_form.cgi'>convert-matrix</i></a>).";

print '
      </div>
    </div>
  </div>
</div>
                </div></div>
                </div>

 <!-- ################################################################-->
 <!-- ### advanced options ###-->
 <!-- ADVANCED OPTIONS -->

                <div class="bhoechie-tab-content">

                  <div id="accordion" role="tablist">


 <!-- Input squences -->
  <div class="card">
    <div class="card-header" role="tab" id="headingThree">
      <h5 class="mb-0">  <i class="fa fa-tasks"></i>
        <a class="collapsed" data-toggle="collapse" href="#collapseThree" aria-expanded="false" aria-controls="collapseThree">
          Options for reducing peak sequences
        </a>
      </h5>
    </div>
    <div id="collapseThree" class="collapse" role="tabpanel" aria-labelledby="headingThree" data-parent="#accordion">
      <div class="card-body">
        <div class="panel panel-warning">
          <div class="panel-heading">Reduce peak sequences
          <i class="fa fa-info-circle" data-container="body" data-toggle="tooltip" data-placement="top" title="Options to filter or trim all the set of input sequences."  data-original-title=""></i>
          </div>
			       <div class="panel-body">';
 # Number of top sequences to keep
 print "      <div class='form-row'>
                <label for='metric' class='col-sm-9 control-label'>Number of top sequences to retain <a class='badge badge-primary iframe' href='help.peak-motifs.html#thresholds'>Info</a></label>\n";
 print "        <div class='col-sm-3'>";

 print  $query->textfield(-name=>'top_sequences',
                          -id=>'top_sequences',
                          -class=>'form-control',
                          -default=>$default{top_sequences},
                          -size=>3);
print "         </div>
              </div>\n";


# Trim peak sequences
print "       <div class='form-row'>
                <label for='merge_stat' class='col-sm-9 control-label'>  Cut peak sequences to retain a given number of base-pairs on each side of the peak centers <a class='badge badge-primary iframe' HREF='help.peak-motifs.html#thresholds'>Info</a>  </label>\n";
print "         <div class='col-sm-3'>";

print  $query->textfield(-name=>'max_seq_len',
                        -id=>'max_seq_len',
                        -class=>'form-control',
			                  -default=>$default{max_seq_len},
			                  -size=>3);
print '
      </div>
    </div>
  </div>
    </div>
      </div>
        </div>
</div>';

print '
<!-- Locating motif in peaks-->
 <div class="card">
   <div class="card-header" role="tab" id="headingFour">
     <h5 class="mb-0">  <i class="fa fa-tasks"></i>
       <a class="collapsed" data-toggle="collapse" href="#collapseFour" aria-expanded="false" aria-controls="collapseFour">
         Options for locating the discovered motifs in the sequences
       </a>
     </h5>
   </div>
   <div id="collapseFour" class="collapse" role="tabpanel" aria-labelledby="headingFour" data-parent="#accordion">
     <div class="card-body">
       <div class="panel panel-warning">
         <div class="panel-heading">Locate motifs
          <i class="fa fa-info-circle" data-container="body" data-toggle="tooltip" data-placement="top" title="Options to locate sites in the peaks of the discovered motifs."  data-original-title=""></i>
         </div>
            <div class="panel-body">';

            ### matrix-scan
            print '<div class="form-inline">';
            print '<div class="form-group">';
            #print "<fieldset><legend><b><a class='iframe' href='help.peak-motifs.html#tasks'>Locate motifs </a></b></legend>";
            print $query->checkbox(-id=>'matrix-scan-quick',
                                   -name=>'matrix-scan-quick',
                                   -class=>'form-check-label',
          			                   -checked=>$default{'matrix-scan-quick'},
          			                   -label=>'');
                                   print "</div>";

            #print "&nbsp;\n";
            print "<label for='matrix-scan-quick' >Search putative binding sites in the peak sequences <a class='badge badge-primary iframe' href='help.matrix-scan.html'>Info</a></label>\n";
            print "</div>";

            #print "<br>";

            print '<div class="form-group">';
            print "<label for='scan_markov'>
                      <b>Markov order</b> (m) of the background model for sequence scanning <br>(site prediction + motif enrichment)
                      <a class='badge badge-primary iframe' href='help.peak-motifs.html#tasks'>Info</a>
                   </label>";
            ## Markov order for scanning (site prediction + motif enrichment)
            print $query->popup_menu(-id=>'scan_markov',
                                    -name=>'scan_markov',
             			                  -Values=>["0", "1", "2", "3"],
                                    -labels=>{"0"=>"m=0 (generally not ideal)",
                                              "1"=>"m=1",
                                              "2"=>"m=2 (slower)",
                                              "3"=>"m=3 (may be very slow)"},
             			                  -class=>'form-control',
             			                  -default=>$default{scan_markov});
            print '
                  </div>';
            print '<!--## Plotting options-->
                    <div class="form-group">
                      <b>Plotting options for motif sites:</b><br>
                      <div class="form-inline">';
            print "     <label class='form-check-label'>Origin</label>\n";
            print "&nbsp;";
            print $query->popup_menu( -name=>'origin',
                                      -class=>'form-control',
                                      -Values=>['start',
                                                'center',
                                                'end'],
                                      -default=>$default{origin});
            print "&nbsp;";
            print '<a class="badge badge-primary iframe" href="help.peak-motifs.html#thresholds">Info</a>';

            print "&nbsp;&nbsp;&nbsp;&nbsp;<label class='form-check-label'>Offset</label>\n";
            print "&nbsp;";

                print $query->textfield(
                      -name=>'offset',
                      -class=>'form-control',
               			  -default=>$default{offset},
              			  -size=>8);
            print "&nbsp;";

            print '<a class="badge badge-primary iframe" href="help.peak-motifs.html#thresholds">Info</a>';
            print '</div>';
            ## Use R to generate XY plots.
            if ($ENV{R_supported}) {
                ## This option is displayed only if R is supported on the server.
                print "<br>";
                print $query->checkbox(-name=>'r_plot',
                                       -class=>'form-check-input',
                                       -checked=>$default{"r_plot"},
                                       -label=>'');
                print "&nbsp;<b>Use R to generate plots</b> (only works for servers with R installed).\n";

            }

            print "<hr>\n";
print '
    </div>
   </div>
   </div>
 </div>
   </div>
</div>';

#print'
#<!-- Exporting motif sites-->
#<div class="card" id="adv_opt_sites">
#  <div class="card-header" role="tab" id="headingFive">
#    <h5> <i class="fa fa-tasks"></i>
#      <a data-toggle="collapse" href="#collapseFive" aria-expanded="true" aria-controls="collapseFive">
#       Options for exporting predicted sites in the peaks
#      </a>
#    </h5>
#  </div>

#  <div id="collapseFive" class="collapse" role="tabpanel" aria-labelledby="headingFive" data-parent="#accordion">
#    <div class="card-body">
#    <div class="panel panel-warning">
#    <div class="panel-heading">Export peaks and sites
#      <i class="fa fa-info-circle" data-container="body" data-toggle="tooltip" data-placement="top" title="Options to export found sites in peaks for visualization in Genome Browsers."  data-original-title=""></i>
#    </div>
#    <div class="panel-body">';
#    ## Visualize UCSC custom track
#    print '<div class="form-group">';
#    #print '<div class="radio">';
#    print '<label>';
#    print ("<INPUT TYPE='radio' NAME='visualize' id='visualize_none' value='none' $checked{'none'}>","<b>No</b>");
#    print '</label>';
#    #print '</div>';
#    #print "<br/>";

#    #print '<div class="radio">';
#    print '<label>';
#    print ("<INPUT TYPE='radio' NAME='visualize' id='visualize_getfasta' value='getfasta' $checked{'getfasta'}>",
#  	 "Peak coordinates specified in fasta headers in <a target='_blank' href='http://bedtools.readthedocs.org/en/latest/content/tools/getfasta.html'><b>bedtools getfasta</b></a> format (also for","<br>","&nbsp;"x7,"<b>retrieve-seq-bed</b> output).",
#  	 "&nbsp;Fasta headers should be in the form: <tt>>3:81458-81806(.)</tt>");
#    print '</label>';
#    #print '</div>';
#    #print "<br/>";
#
#    #print '<div class="radio">';
#    print '<label>';
#    print ("<INPUT TYPE='radio' NAME='visualize' id='visualize_galaxy' value='galaxy' $checked{'galaxy'}>",
#  	 "Peak coordinates specified in fasta headers of the test sequence file (<a target='_blank' href='https://usegalaxy.org/'><b>Galaxy</b></a> format).",
#  	 "<br>","&nbsp;"x7,"Fasta headers should be in the form: <tt>>mm9_chr1_3473041_3473370_+</tt>");
#    print '</label>';
#    #print '</div>';
#    #print "<br/>";
#
#    #print '<div class="radio">';
#    print '<label>';
#    print ("<INPUT TYPE='radio' NAME='visualize' id='visualize_bed_coord' value='bed_coord' $checked{'bed_coord'}>","Peak coordinates provided as a <b>custom BED file.</b>");
#    print "&nbsp;The 4th column of the BED file","<br>","&nbsp;"x7,"(feature name) must correspond to the fasta headers of sequences.";
#    print '</label>';
#
#    ### assembly
#    print '<div class="form-inline">';
#
#    print "<label for='asssemmbly'>","&nbsp;"x7,"<b>Assembly version (UCSC)&nbsp;</b><label>\n";
#    print '<div class="form-group">';
#    print  $query->textfield(-name=>'assembly',
#                            -class=>'form-control',
#  							            -default=>$default{assembly},
#  							            -size=>10);
#    print '</div>&nbsp;';
#    print '<div class="form-group">';
#    print  $query->filefield(-name=>'bed_file',
#                             -class=>'form-control-file',
#                             -size=>10);
#    print '</div>';
#    print '</div>';
#
#    print '</div>';
#
#

# Selection of output fields and thresholds

#print '</div></div>
#    </div>
#  </div>
#  </div>';

print '                </div></div>';


#<!--close panel-->
#</div></div>
# </div>
print '
 <!--################################################################-->
 <!--### output & run ###-->

 <div class="bhoechie-tab-content">

  <!-- ## Specific options for output files-->
  <div class="panel panel-info">
    <div class="panel-heading">Output options</div>

      <div class="panel-body">
        <div class="form-group">

        <!--A class="badge badge-primary iframe" HREF="help.matrix-clustering.html#merge_operator">Info</a></label></h5-->';



#print "<hr>\n";

################################################################
## Send results by email only
print "<p>\n";
&SelectOutput();
print " </div>
  </div> </div>";


################################################################
## Action buttons
print '<script> function formreset(){
demo_descr.innerHTML = "";
} </script>';
print $query->submit(-label=>"GO", -class=>"btn btn-success", -type=>"button");
print " ";
print "<input type='reset' id='reset' class='btn btn-warning' onclick='formreset()' value='RESET'>";
print $query->end_form;
print ' </div>
</div>
</div>
';

# NOTE WSG Might be missing one-closing "</div>"?
################################################################
## Demo area
print "<textarea id='demo' style='display:none'></textarea>";
print "<div id='demo_descr' class='col-lg-9 col-md-5 col-sm-8 col-xs-9 demo-buttons-container'></div>";



#
################################################################
## Data for the demo single-set analysis
$demo_url = $ENV{rsat_www}."/demo_files/ChIP-seq_peaks/Oct4_peaks_top1000.fa";

print '<script>
function setDemo1(demo_url){
    $("#reset").trigger("click");
    $("#dbs_choice").val("'.$default{compare_motif_database}.'").trigger("change");
    setTimeout(setDemo_motif, 1000);
    descr = "<H4>Demonstration: input = one dataset</H4>\n";
    descr = descr + "<blockquote class =\'blockquote text-justify small\'>";
    descr = descr + "In this demonstration, we apply time- and memory-efficient \
    motif discovery algorithms to discover over-represented motifs in a \
    set of 1000 peak regions bound by the mouse transcription factor Oct4 \
    (Chen et al., 2008).Check the panel <b>Mandatory inputs</b> and then <b>Run analysis</b>.</p>\n</blockquote>";

    demo_descr.innerHTML = descr;
    demo.value = descr;
    fastasequence_url1.value = demo_url;
    $("#title").val("Oct4 Chen2008 sites from Jaspar");
    max_seq_len.value = "";
    top_sequences.value = "";
    $("#visualize_galaxy").prop("checked", true);

}

function setDemo_motif(){
    $("#db_choice").val("'.$default{compare_motif_collection}.'").trigger("change");
}
</script>';

print ' <div class="col-lg-9 col-md-5 col-sm-8 col-xs-9 demo-buttons-container">

<button type="button" class="btn btn-info" onclick="setDemo1('. "'$demo_url'".')">DEMO single</button> ';



################################################################
## Data for the demo of differential analysis (test vs control)
$demo_url = $ENV{rsat_www}."/demo_files/peak-motifs_GSM559652_heart_p300_1000peaks.fa";
$ctrl_url = $ENV{rsat_www}."/demo_files/peak-motifs_GSM348066_limb_p300_1000peaks.fa";


print '<script>
function setDemo2(demo_url,ctrl_url){
    $("#reset").trigger("click");
    $("#dbs_choice").val("").trigger("change");
    $("#db_choice").val("").trigger("change");
    descr = "<H4>Demonstration: input = two datasets (control sequence and test sequences)</H4>\n";
    descr = descr + "<blockquote class =\'blockquote text-justify small\'>";
    descr = descr + "In this demonstration, we run a differential analysis \
    (test vs control) to discover the motifs that are over-represented in \
    one tissue (heart) compared to another tissue (limb), for a same \
    transcription factor (p300) (Blow et al, 2010).Check the panel <b>Mandatory \
    inputs</b> and then <b>Run analysis</b></p>\n</blockquote>";

    demo_descr.innerHTML = descr;
    demo.value = descr;
    fastasequence_url1.value = demo_url;
    fastasequence_url2.value = ctrl_url;
    $("#title").val("p300 heart versus limb Blow2010");
    max_seq_len.value = "";
    $("#position-analysis").prop("checked",false);
    $("#oligo-analysis").prop("checked",true);
    $("#oligo_length6").prop("checked",true);
    $("#oligo_length7").prop("checked",true);
    $("#oligo_length8").prop("checked",false);
    $("#nmotifs").val(5);
    top_sequences.value = "";
    $("#visualize_galaxy").prop("checked", true);

};
</script>';
print '<button type="button"class="btn btn-info" onclick="setDemo2('. "'$demo_url'".','."'$ctrl_url'".')">DEMO peak sequences vs contrasting dataset</button> ';

#print $query->end_form;


################################################################
## JS action functions
print '
<script>
$("#seq_format1").on("change", function(){
    if ($(this).val() == "bed" ) {
      document.getElementById("collapseFasta1")
              .setAttribute("style","display:none");
      document.getElementById("collapseBed1")
              .setAttribute("style","display:block");
    } else if($(this).val() == "fasta" ) {
      document.getElementById("collapseBed1")
              .setAttribute("style","display:none");
      document.getElementById("collapseFasta1")
              .setAttribute("style","display:block");
    } else {
      document.getElementById("collapseFasta1")
              .setAttribute("style","display:none");
      document.getElementById("collapseBed1")
              .setAttribute("style","display:none");
    }
});
$("#UCSC_opt").on("click", function(){
      document.getElementById("local_seq1")
              .setAttribute("style","display:none");
      document.getElementById("ucsc_seq1")
              .setAttribute("style","display:block");
      document.getElementById("local_seq2")
              .setAttribute("style","display:none");
      document.getElementById("ucsc_seq2")
              .setAttribute("style","display:block");
});
$("#local_opt").on("click", function(){
      document.getElementById("ucsc_seq1")
              .setAttribute("style","display:none");
      document.getElementById("local_seq1")
              .setAttribute("style","display:block");
      document.getElementById("ucsc_seq2")
              .setAttribute("style","display:none");
      document.getElementById("local_seq2")
              .setAttribute("style","display:block");
});
$("#seq_format2").on("change", function(){
    if ($(this).val() == "bed" ) {
      document.getElementById("collapseFasta2")
              .setAttribute("style","display:none");
      document.getElementById("collapseBed2")
              .setAttribute("style","display:block");

    } else if($(this).val() == "fasta" ) {
      document.getElementById("collapseBed2")
              .setAttribute("style","display:none");
      document.getElementById("collapseFasta2")
              .setAttribute("style","display:block");
    } else {
      document.getElementById("collapseFasta2")
              .setAttribute("style","display:none");
      document.getElementById("collapseBed2")
              .setAttribute("style","display:none");
    }
});
</script>
';




print "</div> ";

print $query->end_html;

exit(0);



################################################################
#################### SUBROUTINE DEFINITIONS  ###################
################################################################

###############################################
## Multiple Coordinates and Sequence choice
sub MultiSeqFormatInput_bootstrap {
      my (%args) = @_;
      my $seq_title = $args{seq_title};
      my $seq_num = $args{seq_num};
      my $seq_type = $args{seq_type};
      my $show_genome= $args{show_genome};
      my $sequence_file = "";
      my @supported_input_formats =  ("- select -","fasta","bed");
      my %genome_param = ('in_title'     => 'Genomic coordinates' ,
                          'in_num'       => $seq_num ,
                          'in_type'      => 'genomic coordinates' ,
                          'in_format'    => 'bed' );

      if (($sequence_file = $query->param("sequence_file".$seq_num)) ||
          ($sequence_file = $query->param("sequence_file"))) {

        ## Sequence file is already on the server machine
        ## (piped from a previous script)
        $sequence_url = $sequence_file;
        $sequence_url =~ s|$ENV{RSAT}/public_html|$ENV{rsat_www}|;
        print "<a href=$sequence_url> transferred from previous query<br></a>";
        $sequence_format = $query->param(sequence_format);
        print "<input type='hidden' id='sequence_format' name='sequence_format' value='".$sequence_format."'>\n";
        print "<input type='hidden' name='sequence_file' value='".$sequence_file."'>\n";

      } else {
        print '
        <div class="form-inline">
          <div class="form-group">
            <label for="seq_format'.$seq_num.'">Sequence Format</label>
            <a class="badge badge-primary iframe" href="help.peak-motifs.html#io_format" >Info</a>';

        #### sequence format
        print $query->popup_menu(-id=>'seq_format'.$seq_num,
                               -name=>'seq_format'.$seq_num,
    			                     -Values=>[@supported_input_formats],
                               -class=>'form-control',
    			                     -default=>$default{"seq_format".$seq_num});


        print '
              </div>
            </div>
            <br>';

        print '
            <div id="collapseInput'.$seq_num.'"  style="display:block">';
        #### Fasta sequence option
        print '
              <div id="collapseFasta'.$seq_num.'" style="display:none">';
        &MultiSequenceChoiceDust_bootstrap('seq_title'  => 'Fasta sequences',
                                           'seq_num'    => $seq_num,
                                           'seq_type'   => 'fasta sequences',
                                           'seq_format' => 'fasta',
                                           'id_prefix' => 'fasta');
        print '</div>';

        #### Bed sequence option
        print '
              <div id="collapseBed'.$seq_num.'" style="display:none">';

        #### Hide choice menu if no genome will be specified
        if ($show_genome){
        print '<div class="form-group">
                <p><b>Retrieve genomic sequences from:</b></p>
                <label>
                  <input type="radio" class="seq_source1" name="seq_source1" id="UCSC_opt" value="UCSC_genome">
                  UCSC organisms
                </label>
                <label>
                  <input type="radio" class="seq_source1" name="seq_source1" id="local_opt" value="local_genome">
                  RSAT locally-installed organisms
                </label>
               </div>';
        }

        #### Fetch sequences from UCSC
        print '<div id="ucsc_seq'.$seq_num.'" style="display:none">';
        $genome_param{'id_prefix'}  = 'UCSC';
        if ($show_genome) {
          $genome_param{'genome_menu'} = 'UCSC';
        }


        &MultiInputChoice_from_Organisms_bootstrap(%genome_param);

        #### Mask optios
        print '<div class="form-inline">';
        print ' <p><b>Mask sequences:</b></p>';
        #### Mask with dust
        print ' <div class="form-group">
                 <label for="mask'.$seq_num.'">Mask low-complexity regions (run Dust program)</label>&nbsp';
        print  $query->textfield(-name=>'dust_seq'.$seq_num,
                                 -class=>'form-control '.'dust_seq'.$seq_num,
                                 -default=>$default{'dust_seq'.$seq_num},
                                 -size=>3);
        print  ' <a class="badge badge-primary iframe" href="help.formats.html" >Info</a>';
        print ' </div>
               </div>
               </div>';

        #### Fetch sequences from locally installed organisms
        print '<div id="local_seq'.$seq_num.'" style="display:none">';
        $genome_param{'id_prefix'}  = 'local';
        if ($show_genome) {
          $genome_param{'genome_menu'}  = 'local';
        }


        &MultiInputChoice_from_Organisms_bootstrap(%genome_param);

        #### Mask optios
        print '<div class="form-inline">';
        print '<p><b>Mask sequences:</b></p>';
        #### Mask with dust
        print '<div class="form-group">
               <label for="mask'.$seq_num.'">Mask low-complexity regions (run Dust program)</label>&nbsp';
        print  $query->textfield(-name=>'dust_seq'.$seq_num,
                                -class=>'form-control '.'dust_seq'.$seq_num,
                                -default=>$default{'dust_seq'.$seq_num},
                                -size=>3);
       print  '<a class="badge badge-primary iframe" href="help.formats.html" >Info</a>';
       print '</div>';
       #### Repeat masking
       print '&nbsp&nbsp&nbsp&nbsp&nbsp';
       print '<div class="form-group">';
       print $query->checkbox(-name=>'rm',
                              -id=>  'rm'.$seq_num,
                              -class=> 'form-check',
                              -checked=>$default{'rm'},
                              -label=>'');
       print  '<label for="rm'.$seq_num.'">Mask repeats</label>&nbsp
               <a class="badge badge-primary iframe" href="help.retrieve-seq-bed.html#rm" >Info</a>';
       # Close repeat masking
       print '</div>';
       # Close  masking options
       print '</div>';
       # Close locally installed genomes
        print  '</div>';
       # Close  bed input
        print '</div>';
        # Close  input
        print '</div>';
      }
}
