#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push @INC, "$`lib/";
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";

$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

### Print the header
&NeAT_header("roc-stats result", "results");

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

#$ENV{rsat_echo}= 2; # TMP
&ListParameters() if ($ENV{rsat_echo} >= 2);

@result_files = ();
$command = "$SCRIPTS/roc-stats2";
#my $tmp_file_prefix = sprintf "roc-stats.%s", &AlphaDate();
$prefix = "roc-stats";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); ($tmp_file_dir, $tmp_file_name) = &SplitFileName($tmp_file_path);

#### read parameters ####
my $parameters = " -v 1 ";
my $img_format = $query->param('img_format')||'png';

################################################################
#### Get input

my $score_file = $tmp_file_path."_input.tab";
push @result_files, "Input score file", $score_file;
my $data = $query->param('data');
if ($data){
  $data =~ s/\r//g;
  open DATA, "> $score_file";
  print DATA $data;
  close DATA;
}elsif($query->param('uploaded_file')){
  $score_file =  $query->param('uploaded_file');
}elsif($query->param('roc-stats_graph_file')){
  $score_file =  $query->param('uploaded_file');
}else{
  &cgiError("You should specify input data");
}
$parameters .= " -i $score_file";

################################################################
#### Get input parameters

## Score Column
if (&IsInteger($query->param('sc_col'))) {
  $parameters .= " -scol ".$query->param('sc_col');
}else{
  &cgiError("You should specify a column for the scores");
}

## Status Column
if (&IsInteger($query->param('status_col'))) {
  $parameters .= " -lcol ".$query->param('status_col');
}else{
  &cgiError("You should specify a column for the status");
}

## Status label
(my $query_pos = $query->param('pos')) =~ s/\r/\n/g;
foreach my $pos_status (split /\n/,$query_pos){
  if (($pos_status =~ /\S/)&&($pos_status ne "pos")) {
    $parameters .= " -status ".$pos_status." pos";
  }
}

(my $query_neg = $query->param('neg')) =~ s/\r/\n/g;
foreach my $neg_status (split /\n/,$query_neg){
  if (($neg_status =~ /\S/)&&($neg_status ne "neg")) {
    $parameters .= " -status ".$neg_status." neg";
  }
}


## total numbers
if (&IsInteger($query->param('total'))) {
  $parameters .= " -total ".$query->param('total');
}

################################################################
## Return fields
#&CGI_return_fields();
my $result_file = $tmp_file_path.".tab";
push @result_files, "ROC statistics (tab)", $result_file;

## graphs 
if ($query->param('graphs')) {
  $parameters .= " -graphs "; # TROUBLESHOOTING ! files are not created on rsat webserver.
  $parameters .= " -img_format $img_format";
  $parameters .= " -o ".$result_file;
}

################################################################
#### Execute
################################################################

## Report the command
&ReportWebCommand($command." ".$parameters);

## execute the command

if ($query->param('graphs')){
  if ($query->param('output') =~ /display/i){
    @data_report = `$command $parameters`;

    print '<H4>Graphs</H4>';
    print "<UL>\n";
    print "<A HREF=\"#scores\">Scores</A><BR>";
    print "<UL>";
    print "<LI><A HREF=\"#scores\">normal scale</A><BR>";
    print "<LI><A HREF=\"#scores_xlog2\">xlog scale</A><BR>";
    print "</UL>";
    print "<A HREF=\"#FP_TP\">FP vs TP</A><BR>";
    print "<A HREF=\"#roc\">ROC (Receiver Operating Characteristic) curve</A><BR>";
    print "<A HREF=\"#precision_recall\">Precision-Recall curve</A><BR>";
    print "<UL>";
    print "<LI><A HREF=\"#precision_recall\">normal scale</A><BR>";
    print "<LI><A HREF=\"#precision_recall_xlog\">xlog scale</A><BR>";
    print "<LI><A HREF=\"#precision_recall_log\">xlog ylog scale</A><BR>";
    print "</UL>";
    print "</UL>\n";
    print "<HR>\n";

    if($query->param('stats')
#        ||
#        $query->param('occ')||
#        $query->param('TP')||
#        $query->param('FP')||
#        $query->param('FN')||
#        $query->param('Sn')||
#        $query->param('PPV')||
#        $query->param('FPR')||
#        $query->param('Acc_g')||
#        $query->param('Acc_a')
      ){
      print '<H4>Table</H4>';
      print @data_report;
      open RESULTS, "<$result_file";
      &PrintHtmlTable(RESULTS, $result_file.".html", true);
      close(RESULTS);
    }
    print "<H3><CENTER>Graphs</CENTER></H3>";

   ## Draw stats as a function of score
    my $cmd = "$SCRIPTS/XYgraph -i ".$result_file;
    $score_distrib_plot = $tmp_file_path."_scores.".$img_format;
    push @result_files, "Score distributions", $score_distrib_plot;
    $cmd .= " -title1 'Score distributions'";
    $cmd .= " -xcol 1 -ycol 7,8,9,10,11 -xleg1 'score' -lines -pointsize 0 -ymin 0 -ymax 1 -legend";
    $cmd .= " -ygstep1 0.1 -ygstep2 0.05 ";
    $cmd .= " -format ".$img_format;
    $cmd .= " -o ".$score_distrib_plot;
    &doit($cmd);
    print "<CENTER><B><A NAME=\"scores\"></A>";
    $img_URL = $ENV{rsat_www}."/tmp/"; $img_URL .= &RSAT::util::RelativePath($TMP, $score_distrib_plot);
    print "<a href='".$img_URL."'><IMG SRC='".$img_URL."'></a><BR>";

    ## Draw stats as a function of score, with log scale on X axis
    $score_distrib_plot_xlog2 = $tmp_file_path."_scores_xlog2.".$img_format;
    push @result_files, "Score distributions (log2 X axis)", $score_distrib_plot_xlog2;
    $cmd = "$SCRIPTS/XYgraph -i ".$result_file;
    $cmd .= " -title1 'Score distributions (xlog)'";
    $cmd .= " -xcol 1 -ycol 7,8,9,10,11 -xleg1 'score' -lines -pointsize 0 -ymin 0 -ymax 1 -legend";
    $cmd .= "  -ygstep1 0.1 -ygstep2 0.05";
    $cmd .= " -format ".$img_format;
    $cmd .= " -xlog 2 -o ".$score_distrib_plot_xlog2;
    &doit($cmd);
    print "<CENTER><B><A NAME=\"scores_xlog2\"></A>";
    $img_URL = $ENV{rsat_www}."/tmp/"; $img_URL .= &RSAT::util::RelativePath($TMP, $score_distrib_plot_xlog2);
    print "<a href='".$img_URL."'><IMG SRC='".$img_URL."'></a><BR>";

    ## Draw a graph with TP=f(FP)
    $FP_TP_plot = $tmp_file_path."_FP_TP.".$img_format;
    push @result_files, "FP versus TP", $FP_TP_plot;
    $cmd = "$SCRIPTS/XYgraph -i ".$result_file;
    $cmd .= " -title1 'True versus false positives'";
    $cmd .= " -xcol 5 -ycol 4 -xleg1 FP -yleg1 TP -lines -pointsize 0";
    $cmd .= " -format ".$img_format;
    $cmd .= " -xlog 2 -o ".$FP_TP_plot;
#    $cmd .= " -o $tmp_file_path."._FP_TP.".$img_format;
    &doit($cmd);
    print "<CENTER><B><A NAME=\"FP_TP\"></A>";
    $img_URL = $ENV{rsat_www}."/tmp/"; $img_URL .= &RSAT::util::RelativePath($TMP, $FP_TP_plot);
    print "<a href='".$img_URL."'><IMG SRC='".$img_URL."'></a><BR>";

    ## Draw a ROC curve
    $ROC_plot = $tmp_file_path."_ROC.".$img_format;
    push @result_files, "ROC curve", $ROC_plot;
    $cmd = "$SCRIPTS/XYgraph -i ".$result_file;
    $cmd .= " -title1 'ROC curve'";
    $cmd .= " -xcol 9 -ycol 7 -xleg1 'FPR' -yleg1 'Sn (=TPR)' -lines -pointsize 0 -min 0 -max 1";
    $cmd .= " -format ".$img_format;
#    $cmd .= " -o $tmp_file_path."_roc.".$img_format;
    $cmd .= " -o ".$ROC_plot;
    &doit($cmd);
    print "<CENTER><B><A NAME=\"roc\"></A>";
    $img_URL = $ENV{rsat_www}."/tmp/"; $img_URL .= &RSAT::util::RelativePath($TMP, $ROC_plot);
    print "<a href='".$img_URL."'><IMG SRC='".$img_URL."'></a><BR>";

    ## Draw a Precision-recall curve
    $precision_recall = $tmp_file_path."_precision_recall.".$img_format;
    push @result_files, "Precision/Recall curve", $precision_recall;
    $cmd = "$SCRIPTS/XYgraph -i ".$result_file;
    $cmd .= " -title1 'Precision-recall curve'";
    $cmd .= " -xcol 7 -ycol 8 -xleg1 'Sn (Recall)' -yleg1 'PPV (Precision)' -lines -pointsize 0 -min 0 -max 1";
    $cmd .= " -format ".$img_format;
    $cmd .= " -o ".$precision_recall;
    &doit($cmd);
    print "<CENTER><B><A NAME=\"precision_recall\"></A>";
    $img_URL = $ENV{rsat_www}."/tmp/"; $img_URL .= &RSAT::util::RelativePath($TMP, $precision_recall);
    print "<a href='".$img_URL."'><IMG SRC='".$img_URL."'></a><BR>";

    ## Draw a Precision-recall curve with logarithmic axis X
    ## (like in von Mering, 2002, but beware: this is over-emphasizing the poor results)
    $precision_recall_xlog = $tmp_file_path."_precision_recall_xlog.".$img_format;
    push @result_files, "Precision/Recall curve (log X axis)", $precision_recall_xlog;
    $cmd = "$SCRIPTS/XYgraph -i ".$result_file;
    $cmd .= " -title1 'Precision-recall curve (xlog)'";
    $cmd .= " -xcol 7 -ycol 8 -xleg1 'Sn (Recall)' -yleg1 'PPV (Precision)' -lines -pointsize 0 -min 0 -max 1 -xlog";
    $cmd .= " -format ".$img_format;
    $cmd .= " -o  ".$precision_recall_xlog;
    &doit($cmd);
    print "<CENTER><B><A NAME=\"precision_recall_xlog\"></A>";
    $img_URL = $ENV{rsat_www}."/tmp/"; $img_URL .= &RSAT::util::RelativePath($TMP, $precision_recall_xlog);
    print "<a href='".$img_URL."'><IMG SRC='".$img_URL."'></a><BR>";

    ## Draw a Precision-recall curve with logarithmic axes
    ## (like in von Mering, 2002, but beware: this is over-emphasizing the poor results)
    $precision_recall_log = $tmp_file_path."_precision_recall_log.".$img_format;
    push @result_files, "Precision/Recall curve (log axes)", $precision_recall_log;
    $cmd = "$SCRIPTS/XYgraph -i ".$result_file;
    $cmd .= " -title1 'Precision-recall curve (log-log)'";
    $cmd .= " -xcol 7 -ycol 8 -xleg1 'Sn (Recall)' -yleg1 'PPV (Precision)' -lines -pointsize 0 -min 0 -max 1 -xlog -ylog";
    $cmd .= " -format ".$img_format;
    $cmd .= " -o ".$precision_recall_log;
    &doit($cmd);
    print "<CENTER><B><A NAME=\"precision_recall_log\"></A>";
    $img_URL = $ENV{rsat_www}."/tmp/"; $img_URL .= &RSAT::util::RelativePath($TMP, $precision_recall_log);
    print "<a href='".$img_URL."'><IMG SRC='".$img_URL."'></a><BR>";

    &PrintURLTable(@result_files);

  }else{
    ## TO BE IMPLEMENTED
    &cgiError("Graph option is not yet supported by email output. Please choose display output.");
  }

} else {
  if ($query->param('output') =~ /display/i){
    open RESULT, "$command $parameters | ";
    print '<H4>Table</H4>';
    &PrintHtmlTable(RESULT, $result_file, true);
    close(RESULT);

  }else{
    my $mail_title = join (";", "NeAT", "roc-stats", &AlphaDate());
    &EmailTheResult("$command $parameters", $query->param('user_email'), $result_file, title=>$mail_title);
  }
}

print $query->end_html();

exit(0);

################################################################
## Concatenate return and threshold options from the CGI form
sub CGI_return_fields {
  my %field_group = (
		     occ=>"N_icum,F_icum,TP_icum,FP_icum,FN_icum",
		     TP=>"TP",
		     FP=>"FP",
		     FN=>"FN",
		     Sn=>"Sn",
		     PPV=>"PPV",
		     FPR=>"FPR",
		     Acc_g=>"Acc_g",
		     Acc_a=>"Acc_a",
		     AUC=>"AUC"
		    );
  my %return_fields = ();
  foreach my $field (sort keys %field_group) {
    my $field_group = $field_group{$field};
    if ($query->param($field_group)) {
      $return_fields{$field_group} = 1;
    }
  }
  my $return_fields = join ",", sort(keys(%return_fields));

  unless ($return_fields) {
    &cgiError("You should select at least one option in the \"Return\" box.");
  } else {
    $parameters .= " -return ".$return_fields;
  }
}


sub NeAT_header {
  my $css_body_class = "form";
  my ($title) = shift;
  $title =~ s/\"//g;
  $title =~ s/\'//g;
  if (scalar @_ > 0) {
    $css_body_class = shift;
  }


#  print &html_header();
  print $query->header();
  print sorttable_script();
  ### print the header of the result page
  print $query->start_html(-title=>"Network Analysis Tools : $title",
			   -class => "$css_body_class",
			   -author=>'jacques.van.helden@ulb.ac.be',
			   -style => { 	-src => "$ENV{rsat_www}/main.css",
                             	       	-type => 'text/css',
                             		-media => 'screen' });
  print "<H3 ALIGN='center'><A HREF='$ENV{rsat_www}/NeAT_home.html'>Network Analysis Tools</A> - $title</H3>";
}
