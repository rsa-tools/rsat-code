#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
$convert_seq_command = "$SCRIPTS/convert-seq";
$tmp_file_name = sprintf "convert-seq.%s", &AlphaDate;

### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("convert-seq result");
#&ListParameters;

#### update log file ####
&UpdateLogFile;

#### read parameters ####
$parameters = "";

##### input format #####
$input_format = lc($query->param('input_format'));
if ($accepted_input_seq{$input_format}) {
    $parameters .= " -from $input_format ";
} else {
    &cgiError("Invalid format for input sequence: $input_format.");
}

##### output format #####
$output_format = lc($query->param('output_format'));
if ($accepted_output_seq{$output_format}) {
    $parameters .= " -to $output_format ";
} else {
    &cgiError("Invalid format for output sequence: $output_format.");
}

### sequence file
#($sequence_file,$sequence_format) = &GetSequenceFile();
#$parameters .= " -i $sequence_file -format $sequence_format";

##### input sequence file #####
unless ($query->param('sequence') =~ /\S/) {
    &cgiError("The sequence box should not be empty.");
}
open INSEQ, ">$TMP/$tmp_file_name";
print INSEQ $query->param('sequence');
close INSEQ;
$parameters .= " -i $TMP/$tmp_file_name ";
DelayedRemoval("$TMP/$tmp_file_name");


##### add reverse-complement #####
if ($query->param('addrc') eq "yes") {
    $parameters .= " -addrc ";
}

##### line width #####
if ($query->param('line_width') =~ /\d+/) {
    $parameters .= " -lw ".$query->param('line_width');
}

### print the header
#  print <<End_Header;
#  <HEADER>
#      <TITLE>RSA-tools - convert sequence result</TITLE>
#      </HEADER><BODY BGCOLOR="#FFFFFF">
#      <H3 ALIGN=CENTER><A HREF="$WWW_RSA/RSA_home.cgi">
#      RSA-tools</A> - convert sequence result</H3>
#  End_Header


    ### execute the command ###
    if ($query->param('output') eq "display") {
	### Print the result on Web page
	open RESULT, "$convert_seq_command $parameters  & |";

#print "<PRE>$convert_seq_command $parameters</PRE> ";

	print "<PRE>";
	while (<RESULT>) {
	    print "$_";
	}
	print "</PRE>";
	close RESULT;

	print "<HR SIZE = 3>";
    } else {
	&EmailTheResult( "$convert_seq_command $parameters", , $query->param('user_email'));
	
# 	### send an e-mail with the result ###
# 	if ($query->param('user_email') =~ /(\S+\@\S+)/) {
# 	    $address = $1;
# 	    print "<B>Result will be sent to your e-mail address: <P>";
# 	    print "$address</B><P>";
# 	    system "$convert_seq_command $parameters | $mail_command $address &";
# 	} else {
# 	    if ($query->param('user_email') eq "") {
# 		print "<B>ERROR: you did not enter your e-mail address<P>";
# 	    } else {
# 		print "<B>ERROR: the e-mail address you entered is not valid<P>";
# 		print $query->param('user_email')."</B><P>";      
# 	    }
# 	} 
    }
print $query->end_html;

exit(0);


