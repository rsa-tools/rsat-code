#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push @INC, "$`lib/";
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";

$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
$command = "$SCRIPTS/roc-stats2";

### Read the CGI query
$query = new CGI;

### Print the header
&RSA_header("roc-stats result", "results");

#### update log file ####
&UpdateLogFile();

#$ENV{rsat_echo}= 2; # TMP
&ListParameters() if ($ENV{rsat_echo} >= 2);

#### read parameters ####
my $parameters = " -v 1 ";

################################################################
#### Get input
my $tmp_file_prefix = sprintf "roc-stats.%s", &AlphaDate();
my $result_file = "$TMP/$tmp_file_prefix.res";
my $score_file = "$TMP/$tmp_file_prefix.input";
&cgiError("You should specify scores data") unless($scores = $query->param('scores'));
my $scores = $query->param('scores');
$scores =~ s/\r//g;
open SCORES, "> $score_file";
print SCORES $scores;
close SCORES;
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
if ($query->param('pos') ne "pos") {
  $parameters .= " -status ".$query->param('pos')." pos";
}
if ($query->param('neg') ne "neg") {
  $parameters .= " -status ".$query->param('neg')." neg";
}

## total numbers
if (&IsInteger($query->param('total'))) {
  $parameters .= " -total ".$query->param('total');
}

################################################################
## Return fields
#&CGI_return_fields();

## graphs 
if ($query->param('graphs')) {
  $parameters .= " -graphs "; # TROUBLESHOOTING ! files are not created on rsat webserver.
  $parameters .= " -img_format ".$query->param('img_format');
  $parameters .= " -o ".$result_file;
}

################################################################
#### Execute
################################################################

## Report the command
print "<PRE>$command $parameters </PRE>" if ($ENV{rsat_echo} >= 1);

## execute the command
if ($query->param('output') =~ /display/i){
  if ($query->param('graphs')){
    @data_report = `$command $parameters`;
    print '<H4>Result</H4>';
    print @data_report;
    open RESULTS, "<$result_file";
    &PrintHtmlTable(RESULTS, $result_file.".html", true);
    close(RESULTS);
    my $result_prefix = $tmp_file_prefix.".res";
    my $img_format = $query->param('img_format');

    print "<H3><CENTER>Graphs</CENTER></H3>";
    print "<A HREF=\"#scores\">scores<BR>";
    print "<A HREF=\"#scores_xlog2\">scores (xlog)<BR>";
    print "<A HREF=\"#FP_TP\">FP vs TP<BR>";
    print "<A HREF=\"#roc\">ROC (Receiving Operator Characteristic) curve<BR>";
    print "<A HREF=\"#precision_recall\">PR (Precision-Recall) curve<BR>";
    print "<A HREF=\"#precision_recall_xlog\">PR (xlog)<BR>";
    print "<A HREF=\"#precision_recall_log\">PR (log log)<BR>";
    print "<HR>\n";

   ## Draw stats as a function of score
    my $cmd = "$SCRIPTS/XYgraph -i ".$result_file;
    $cmd .= " -title1 'Score distributions'";
    $cmd .= " -xcol 1 -ycol 7,8,9,10,11 -xleg1 'score' -lines -pointsize 0 -ymin 0 -ymax 1 -legend";
    $cmd .= " -ygstep1 0.1 -ygstep2 0.05 ";
    $cmd .= " -format ".$img_format;
    $cmd .= " -o $TMP/".$result_prefix."_scores.".$img_format;
    &doit($cmd);
    print "<CENTER><B><A NAME=\"scores\"></A>";
    print "<IMG SRC=\"$WWW_TMP/".$result_prefix."_scores.".$img_format."\"><BR>";

    ## Draw stats as a function of score, with log scale on X axis
    $cmd = "$SCRIPTS/XYgraph -i ".$result_file;
    $cmd .= " -title1 'Score distributions (xlog)'";
    $cmd .= " -xcol 1 -ycol 7,8,9,10,11 -xleg1 'score' -lines -pointsize 0 -ymin 0 -ymax 1 -legend";
    $cmd .= "  -ygstep1 0.1 -ygstep2 0.05";
    $cmd .= " -format ".$img_format;
    $cmd .= " -xlog 2 -o $TMP/".$result_prefix."_scores_xlog2.".$img_format;
    &doit($cmd);
    print "<CENTER><B><A NAME=\"scores_xlog2\"></A>";
    print "<IMG SRC=\"$WWW_TMP/".$result_prefix."_scores_xlog2.".$img_format."\"><BR>";

    ## Draw a graph with TP=f(FP)
    $cmd = "$SCRIPTS/XYgraph -i ".$result_file;
    $cmd .= " -title1 'True versus false positives'";
    $cmd .= " -xcol 5 -ycol 4 -xleg1 FP -yleg1 TP -lines -pointsize 0";
    $cmd .= " -format ".$img_format;
    $cmd .= " -o $TMP/".$result_prefix."_FP_TP.".$img_format;
    &doit($cmd);
    print "<CENTER><B><A NAME=\"FP_TP\"></A>";
    print "<IMG SRC=\"$WWW_TMP/".$result_prefix."_FP_TP.".$img_format."\"><BR>";

    ## Draw a ROC curve
    $cmd = "$SCRIPTS/XYgraph -i ".$result_file;
    $cmd .= " -title1 'ROC curve'";
    $cmd .= " -xcol 9 -ycol 7 -xleg1 'FPR' -yleg1 'Sn (=TPR)' -lines -pointsize 0 -min 0 -max 1";
    $cmd .= " -format ".$img_format;
    $cmd .= " -o $TMP/".$result_prefix."_roc.".$img_format;
    &doit($cmd);
    print "<CENTER><B><A NAME=\"roc\"></A>";
    print "<IMG SRC=\"$WWW_TMP/".$result_prefix."_roc.".$img_format."\"><BR>";

    ## Draw a Precision-recall curve
    $cmd = "$SCRIPTS/XYgraph -i ".$result_file;
    $cmd .= " -title1 'Precision-recall curve'";
    $cmd .= " -xcol 7 -ycol 8 -xleg1 'Sn (Recall)' -yleg1 'PPV (Precision)' -lines -pointsize 0 -min 0 -max 1";
    $cmd .= " -format ".$img_format;
    $cmd .= " -o $TMP/".$result_prefix."_precision_recall.".$img_format;
    &doit($cmd);
    print "<CENTER><B><A NAME=\"precision_recall\"></A>";
    print "<IMG SRC=\"$WWW_TMP/".$result_prefix."_precision_recall.".$img_format."\"><BR>";

    ## Draw a Precision-recall curve with logarithmic axis X
    ## (like in von Mering, 2002, but beware: this is over-emphasizing the poor results)
    $cmd = "$SCRIPTS/XYgraph -i ".$result_file;
    $cmd .= " -title1 'Precision-recall curve (xlog)'";
    $cmd .= " -xcol 7 -ycol 8 -xleg1 'Sn (Recall)' -yleg1 'PPV (Precision)' -lines -pointsize 0 -min 0 -max 1 -xlog";
    $cmd .= " -format ".$img_format;
    $cmd .= " -o $TMP/".$result_prefix."_precision_recall_xlog.".$img_format;
    &doit($cmd);
    print "<CENTER><B><A NAME=\"precision_recall_xlog\"></A>";
    print "<IMG SRC=\"$WWW_TMP/".$result_prefix."_precision_recall_xlog.".$img_format."\"><BR></CENTER>";

    ## Draw a Precision-recall curve with logarithmic axes
    ## (like in von Mering, 2002, but beware: this is over-emphasizing the poor results)
    $cmd = "$SCRIPTS/XYgraph -i ".$result_file;
    $cmd .= " -title1 'Precision-recall curve (log-log)'";
    $cmd .= " -xcol 7 -ycol 8 -xleg1 'Sn (Recall)' -yleg1 'PPV (Precision)' -lines -pointsize 0 -min 0 -max 1 -xlog -ylog";
    $cmd .= " -format ".$img_format;
    $cmd .= " -o $TMP/".$result_prefix."_precision_recall_log.".$img_format;
    &doit($cmd);
    print "<CENTER><B><A NAME=\"precision_recall_log\"></A>";
    print "<IMG SRC=\"$WWW_TMP/".$result_prefix."_precision_recall_log.".$img_format."\"><BR>";

  }else{
    open RESULT, "$command $parameters | ";
    print '<H4>Result</H4>';
    &PrintHtmlTable(RESULT, $result_file, true);
    close(RESULT);
  }
} else {
  my $mail_title = join (";", "NeAT", "roc-stats", &AlphaDate(), $tmp_file_prefix);
  &EmailTheResult("$command $parameters", $query->param('user_email'), $result_file,
		  title=>$mail_title,
		 );
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
