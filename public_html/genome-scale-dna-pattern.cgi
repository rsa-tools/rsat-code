#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib.pl";
require "RSA.cgi.lib.pl";

$retrieve_seq_command = "$SCRIPTS/retrieve-seq";
$dna_pattern_command = "$SCRIPTS/dna-pattern";
$add_orf_function_command = "$SCRIPTS/add-orf-function";
$link_command = "$SCRIPTS/add-yeast-link";
$tmp_file_name = sprintf "genome-scale-dna-pattern.%s", &AlphaDate;

### Read the CGI query
$query = new CGI;

### print the header of the result page
&RSA_header("dna-pattern result ".$query->param("title"));

#&ListParameters;


#### update log file ####
&UpdateLogFile;

################################################################
#
# retrieve-seq parameters
#
$retrieve_seq_parameters = " -all ";

#### organism
if (defined($supported_organism{$query->param('organism')})) {
    $org = $query->param('organism');
} else {
    $org = "yeast";
}
$retrieve_seq_parameters .= " -org ".$org;

### sequence type
if ($query->param('sequence_type')) {
  $retrieve_seq_parameters .= " -type ".$query->param('sequence_type');
}

### output format ###
if ($accepted_output_seq{$out_format}) {
  $retrieve_seq_parameters .= " -format ".$query->param('sequence_format');;
}

### sequence label
$seq_label = lc($query->param('seq_label'));
if (($seq_label =~ /gene/) && 
    ($seq_label =~ /orf/)) {
  $retrieve_seq_parameters .= " -label orf_gene";
} elsif ($seq_label =~ /gene/) {
  $retrieve_seq_parameters .= " -label gene";
} elsif ($seq_label =~ /orf/) {
  $retrieve_seq_parameters .= " -label orf";
} elsif ($seq_label =~ /full/) {
  $retrieve_seq_parameters .= " -label full";
} else {
  &cgiError("Invalid option for sequence label '$seq_label'");
}

### limits ###
if (&IsInteger($query->param('from'))) {
  $retrieve_seq_parameters .= " -from ".$query->param('from');
}  
if (&IsInteger($query->param('to'))) {
  $retrieve_seq_parameters .= " -to ".$query->param('to');
}

### orf overlap ###
unless (lc($query->param('orf_overlap')) eq "on") {
  $retrieve_seq_parameters .= " -noorf ";
}


################################################################
#
# dna-pattern parameters
#

#### read parameters ####
$parameters_dna_pattern = "";
$parameters_dna_pattern .= " -v";
$parameters_dna_pattern .= " -format ".$query->param('sequence_format');

### pattern file ####
unless ($query->param('patterns') =~ /\S/) {
  &cgiError("The pattern box should not be empty.<P>Read on-line manual for more information.");
}
$pattern_file = "$TMP/$tmp_file_name.pat";
if (open PAT, ">$pattern_file") {
  print PAT $query->param('patterns');
  close PAT;
  &DelayedRemoval($pattern_file);
}
$parameters_dna_pattern .= " -pl $pattern_file";


### return match count ###
if ($query->param('return') =~ /count/i) {
  $parameters_dna_pattern .= " -c";
  if (($query->param('threshold') =~ /^\d+$/) && ($query->param('threshold') > 0)) {
    $parameters_dna_pattern .= " -th ".$query->param('threshold');
  }


### return match count table
} elsif ($query->param('return') =~ /table/i) {
  $parameters_dna_pattern .= " -table";
  ### add a rwo and a column with the totals
  if (lc($query->param('total')) eq "on") {
    $parameters_dna_pattern .= " -total";
  }
  
### return matching positions
} elsif ($query->param('return') =~ /positions/) { 
  ### origin ###
  if ($query->param('origin') =~ /end/i) {
    $parameters_dna_pattern .= " -origin -0";
  }
  
  if ($query->param('flanking') =~ /^\d+$/) {
    $parameters_dna_pattern .= " -N ".$query->param('flanking');
  }
}

### prevent overlapping matches
if (lc($query->param('noov')) eq "on") {
  $parameters_dna_pattern .= " -noov";
}


### strands
if ($query->param('strands') =~ /direct/i) {
  $parameters_dna_pattern .= " -D";
} elsif  ($query->param('strands') =~ /reverse/i) {
  $parameters_dna_pattern .= " -R";
}

### substitutions
if ($query->param('subst') =~ /^\d+$/) {
  $parameters_dna_pattern .= " -subst ".$query->param('subst');
}


$command = "$retrieve_seq_command $retrieve_seq_parameters ";
$command .= "| $dna_pattern_command $parameters_dna_pattern ";
if ($org eq "yeast") { #### not yet supported for other organisms
    $command .= "| $add_orf_function_command  ";
    $command .= "| $link_command  ";
}

### execute the command ###
if ($query->param("output") =~ /display/i) {

    ### execute the command ###
    $result_file = "$TMP/$tmp_file_name.res";
#  print "<PRE>$command</PRE>";
    open RESULT, "$command & |";
    
    &PipingForm if ($query->param('return') =~ /positions/);
    
    ### Print the result on Web page
    &PrintHtmlTable(RESULT, $result_file);
    close RESULT;
    
    
} else {
    ### send an e-mail with the result ###
    if ($query->param('user_email') =~ /(\S+\@\S+)/) {
	$address = $1;
	print "<B>Result will be sent to your e-mail address: <P>";
	print "$address</B><P>";
	system "$command | $mail_command $address &";
    } else {
	if ($query->param('user_email') eq "") {
	    print "<B>ERROR: you did not enter your e-mail address<P>";
	} else {
	    print "<B>ERROR: the e-mail address you entered is not valid<P>";
	    print "$query->param('user_email')</B><P>";      
	}
    } 
    print "<HR SIZE = 3>";
}

print $query->end_html;


exit(0);

sub PipingForm {
  ### prepare data for piping
  $title = $query->param("title");
  $title =~ s/\"/'/g;
  print <<End_of_form;
<CENTER>
<TABLE>
<TR>
<TD>
<H4>Next step</H4>
</TD>
<TD>
<FORM METHOD="POST" ACTION="feature-map_form.cgi">
<INPUT type="hidden" NAME="title" VALUE="$title">
<INPUT type="hidden" NAME="feature_file" VALUE="$result_file">
<INPUT type="hidden" NAME="format" VALUE="dna-pattern">
<INPUT type="hidden" NAME="fill_form" VALUE="on">
<INPUT type="submit" VALUE="feature map">
</FORM>
</TD>
</TR>
</TABLE>
</CENTER>
End_of_form
}
