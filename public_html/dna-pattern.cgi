#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

$dna_pattern_command = "$SCRIPTS/dna-pattern";
$tmp_file_name = sprintf "dna-pattern.%s", &AlphaDate;

### Read the CGI query
$query = new CGI;

### print the header of the result page
&RSA_header("dna-pattern result ".$query->param("title"));

&ListParameters if ($ECHO >= 2);


#### update log file ####
&UpdateLogFile;

#### read parameters ####
$parameters = "";
$parameters .= " -v";

### pattern file ####
unless ($query->param('patterns') =~ /\S/) {
  &cgiError("The pattern box should not be empty.");
}
$pattern_file = "$TMP/$tmp_file_name.pat";
if (open PAT, ">$pattern_file") {
  print PAT $query->param('patterns');
  close PAT;
  &DelayedRemoval($pattern_file);
}
$parameters .= " -pl $pattern_file";


### sequence file
($sequence_file,$sequence_format) = &GetSequenceFile();
$parameters .= " -i $sequence_file -format $sequence_format";


### return match positions ###
if ($query->param('match_positions')) {
    ### origin ###
    if ($query->param('origin') =~ /end/i) {
	$parameters .= " -origin -0";
    }
    
    if ($query->param('flanking') =~ /^\d+$/) {
	$parameters .= " -N ".$query->param('flanking');
    }
    $parameters .= " -pos";
} 

### return match count ###
if ($query->param('match_counts')) {
    $parameters .= " -c";
    if (($query->param('threshold') =~ /^\d+$/) && ($query->param('threshold') > 0)) {
	$parameters .= " -th ".$query->param('threshold');
    }
} 

### return match count table
if ($query->param('table')) {
    $parameters .= " -table";
    ### add a rwo and a column with the totals
    if ($query->param('total')) {
	$parameters .= " -total";
    }
}

### return match statistics
if ($query->param('stats')) {
  $parameters .= " -stats";
}
  
### prevent overlapping matches
if (lc($query->param('noov')) eq "on") {
  $parameters .= " -noov";
}


### strands
if ($query->param('strands') =~ /direct/i) {
  $parameters .= " -D";
} elsif  ($query->param('strands') =~ /reverse/i) {
  $parameters .= " -R";
}

### substitutions
if ($query->param('subst') =~ /^\d+$/) {
  $parameters .= " -subst ".$query->param('subst');
}


### execute the command ###
if ($query->param("output") =~ /display/i) {
    
    ### execute the command ###
    &PipingWarning() if ($query->param('match_positions'));
    
    $result_file = "$TMP/$tmp_file_name.res";
    open RESULT, "$dna_pattern_command $parameters |";
    print "<PRE>$dna_pattern_command $parameters </b>" if ($ECHO);
  
    ### Print the result on Web page
    print "<H3>Result</H3>";
    PrintHtmlTable(RESULT, $result_file);
    close RESULT;
    &PipingForm if ($query->param('match_positions'));
    print "<HR SIZE = 3>";

} else {
    &EmailTheResult("$dna_pattern_command $parameters", $query->param('user_email'));
#   ### send an e-mail with the result ###
#     if ($query->param('user_email') =~ /(\S+\@\S+)/) {
# 	$address = $1;
# 	print "<B>Result will be sent to your e-mail address: <P>";
# 	print "$address</B><P>";
# 	system "$dna_pattern_command $parameters | $mail_command $address &";
#     } else {
# 	if ($query->param('user_email') eq "") {
# 	    print "<B>ERROR: you did not enter your e-mail address<P>";
# 	} else {
# 	    print "<B>ERROR: the e-mail address you entered is not valid<P>";
# 	    print "$query->param('user_email')</B><P>";      
# 	}
#     } 
}

print $query->end_html;


exit(0);

sub PipingForm {

    ### prepare data for piping
    $title = $query->param("title");
    $title =~ s/\"/\'/g;
    print <<End_of_form;
<HR SIZE = 3>
<CENTER>
<TABLE>
<TR><TD>
<H3>Next step</H3>
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
