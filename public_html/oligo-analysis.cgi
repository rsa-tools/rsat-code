#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
  push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib.pl";
require "RSA.cgi.lib.pl";

$oligo_analysis_command = "$SCRIPTS/oligo-analysis";
$tmp_file_name = sprintf "oligo-analysis.%s", &AlphaDate;

### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("oligo-analysis result");
#&ListParameters;

#### update log file ####
&UpdateLogFile;

#### read parameters ####
$parameters = "";
$parameters .= " -sort";

### sequence file
$sequence_file = &GetSequenceFile;
$parameters .= " -i $sequence_file";


### fields to return
$return_fields = "";
if ($query->param('occ')) {
  $return_fields .= "occ,";
} 
if ($query->param('freq')) {
  $return_fields .= "freq,";
} 
if ($query->param('mseq')) {
  $return_fields .= "mseq,";
} 
if ($query->param('ratio')) {
  $return_fields .= "ratio,";
} 
if ($query->param('zscore')) {
  $return_fields .= "zscore,";
} 
if ($query->param('proba')) {
  $return_fields .= "proba,";
} 
if ($query->param('pos')) {
  $return_fields .= "pos,";
} 
$return_fields =~ s/,$//;
if ($return_fields eq "") {
  &cgiError("Error: you should select at least one option in the \"Return\" box.");
} else {
  $parameters .= " -return $return_fields";
}

### single strand search
if ($query->param('strand') =~ /single/) {
  $parameters .= " -1str";
} else {
  $parameters .= " -2str";
}

### prevent overlapping matches of the same pattern
if ($query->param('noov')) {
  $parameters .= " -noov";
} 

### thresholds ###
if ($query->param('ms_threshold') =~ /^\d+$/) {
  $parameters .= " -thms ".$query->param('ms_threshold');
}  
if ($query->param('occurence_threshold') =~ /^\d+$/) {
  $parameters .= " -tho ".$query->param('occurence_threshold');
}
if ($query->param('proba_occ_threshold') =~ /^[\d\.\-+e]+$/i) {
  $parameters .= " -thpo ".$query->param('proba_occ_threshold');
}
if ($query->param('occ_significance_threshold') =~ /^-{0,1}[\d\.\-+e]+$/i) {
  $parameters .= " -thosig ".$query->param('occ_significance_threshold');
}
if ($query->param('proba_ms_threshold') =~ /^[\d\.\-+e]+$/i) {
  $parameters .= " -thpms ".$query->param('proba_ms_threshold');
}

### graphical output ###
$parameters .= " -v";

#### oligo size ####
if ($query->param('oligo_size') =~ /\d/) {
  $oligo_length = $query->param('oligo_size') ;
} 
$parameters .= " -l $oligo_length";

#### expected frequency estimation ####
if ($query->param('freq_estimate') =~ /oligo freq.* in non-coding regions/i) {
  ### check organism
  unless ($organism = $query->param('organism')) {
    &cgiError("Error : you should specify an organism to use non-coding frequency calibration");
  }
  unless (defined(%{$supported_organism{$organism}})) {
    &cgiError("Error: organism $org is not supported on this site");
  }
  ### select expected frequency file for that organism
  $freq_file = $supported_organism{$organism}->{'data'};
  $freq_file .= "/oligo-frequencies";
  $freq_file .= "/${oligo_length}nt_non-coding_${organism}.freq";
  unless (-r $freq_file) {
    &cgiError("Error: cannot read expected frequency file $freq_file");
  }
  $freq_option = " -expfreq $freq_file";
} elsif ($query->param('freq_estimate') eq "alphabet from input sequence") {
  $freq_option = " -a input";
} elsif ($query->param('freq_estimate') =~ /markov/i) {
  $freq_option = " -markov";
} else {
  $freq_option = "";
} 
$parameters .= "$freq_option";

#### neighborhood ####
if ($query->param('neighborhood') =~ /N at one position/i) {
  $parameters .= " -oneN";
} elsif ($query->param('neighborhood') =~ /one degenerated position/i) {
  $parameters .= " -onedeg ";
}




if ($query->param('output') =~ /display/i) {
  ### execute the command ###
  $result_file = "$TMP/$tmp_file_name.res";
  open RESULT, "$oligo_analysis_command $parameters |";
  
  #print "<PRE>command: oligo-analysis $parameters<P>\n</PRE>";
  
  ### prepare data for piping
  $title = $query->param('title');
  $title =~ s/\"/'/g;
  print <<End_of_form;
<TABLE>
<TR>
<TD>
<H4>Next step</H4>
</TD>
<TD>
<FORM METHOD="POST" ACTION="dna-pattern_form.cgi">
<INPUT type="hidden" NAME="title" VALUE="$title">
<INPUT type="hidden" NAME="pattern_file" VALUE="$result_file">
<INPUT type="hidden" NAME="sequence_file" VALUE="$sequence_file">
<INPUT type="hidden" NAME="sequence_format" VALUE="$sequence_format">
<INPUT type="submit" value="pattern search">
</FORM>
</TD>
</TR>
</TABLE>
End_of_form
  
  ### Print result on the web page
  print '<H2>Result</H2>';
  PrintHtmlTable(RESULT, $result_file);
  close(RESULT);
  
  #### oligonucleotide assembly ####
  if ((&IsReal($query->param('occ_significance_threshold'))) && ($query->param('occ_significance_threshold')>= -1)) {
    $fragment_assembly_command = "$SCRIPTS/pattern-assembly -v";
    if ($query->param('strand') =~ /single/) {
      $fragment_assembly_command .= " -1str";
    } else {
      $fragment_assembly_command .= " -2str";
    }
    
    print "<H2>Fragment assembly</H2>\n";
    open CLUSTERS, "$fragment_assembly_command -i $result_file |";
    print "<PRE>\n";
    while (<CLUSTERS>) {
      print;
    }
    print "</PRE>\n";
    close(CLUSTERS);
  }
  
} else {
  #### send e-mail with the result
  if ($query->param('user_email') =~ /(\S+\@\S+)/) {
    $address = $1;
    print "<B>Result will be sent to your account: <P>";
    print "$address</B><P>";
    system "$oligo_analysis_command $parameters | $mail_command $address &"; 
  } else {
    if ($query->param('user_email') eq "") {
      &cgiError("ERROR: you did not enter your e-mail address");
    } else {
      &cgiError("ERROR: the e-mail address you entered is not valid");
      print $query->param('user_email')."</B><P>";      
    }
  }
  print '<HR SIZE=3>';
}

print $query->end_html;

exit(0);

