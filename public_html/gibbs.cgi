#!/usr/bin/perl
############################################################
#
# $Id: gibbs.cgi,v 1.3 2000/11/12 10:39:43 jvanheld Exp $
#
# Time-stamp: <2000-11-12 11:39:03 jvanheld>
#
############################################################
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}

use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib.pl";
require "RSA.cgi.lib.pl";

$gibbs_command = "nice -n 30 $BIN/gibbs";
$convert_seq_command = "$SCRIPTS/convert-seq";
$tmp_file_name = sprintf "gibbs.%s", &AlphaDate;

### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("gibbs result");
#&ListParameters;

#### update log file ####
&UpdateLogFile;

#### read parameters ####
$parameters = "-v -sort -return proba ";

#### add reverse complement
if (lc($query->param("add_rc")) eq "on") {
    $add_rc = 1;
    $convert_seq_options .= "-addrc ";
} else {
    $add_rc = 0;
}

### sequence file ####
$sequence_file = &GetSequenceFile("fasta",$add_rc);

$parameters .= " $sequence_file ";

### pattern length ###
if (&IsNatural($query->param('length'))) {
    $parameters .= $query->param('length')." ";
}

### expected number of matches
if (&IsNatural($query->param($expected))) {
    $paramaters .= $query->param($expected)." ";
}

### sequence type
if (lc($query->param(seq_type)) eq "dna") {
    $parameters .= "-n ";
}

### inactivate frqgmentation
unless (lc($query->param(fragmentation)) eq "on") {
    $parameters .= "-d ";
}


### execute the command ###
if ($query->param('output') eq "display") {
    ### Print the result on Web page
    open RESULT, "$gibbs_command $parameters & |";
    
    print "<PRE>";
    print "$gibbs_command $parameters &\n";
    while (<RESULT>) {
	print;
    }
    close RESULT;
    print "</PRE>";
    
} else {
  ### send an e-mail with the result ###
    if ($query->param('user_email') =~ /(\S+\@\S+)/) {
	$address = $1;
	print "<B>Result will be sent to your account: <P>";
	print "$address</B><P>";
	system "$gibbs_command $parameters | $mail_command $address &";
    } else {
	if ($query->param('user_email') eq "") {
	    print "<B>ERROR: you did not enter your e-mail address<P>";
	} else {
	    print "<B>ERROR: the e-mail address you entered is not valid<P>";
	    print $query->param('user_email')."</B><P>\n";      
	}
    } 
}
print "<HR SIZE = 3>";
print $query->end_html;

exit(0);
