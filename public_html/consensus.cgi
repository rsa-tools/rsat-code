#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
require "cgi-lib.pl";
require "RSA.lib.pl";
require "RSA_cgi_lib.pl";
$consensus_command = "$BIN/consensus";
$convert_seq_command = "$SCRIPTS/convert-seq";
$tmp_file_name = sprintf "consensus.%s", &AlphaDate;

MAIN:
{
  ### open the Web page for the answer
  print &PrintHeader;

  ### Read the content of the form 
  &ReadParse(*input);

  #### update log file ####
  &UpdateLogFile;

  ### sequence file ####
  unless ($input{'sequence'} =~ /\S/) { ### empty matrix
      &cgiError("Error: you did not enter the sequence");
  }
  $sequence_format = lc($input{'sequence_format'});
  $sequence_file = "$TMP/$tmp_file_name.seq";
  if (($sequence_format ne "wc") && 
      ($sequence_format ne "wconsensus")) { 
      if ($accepted_input_seq{$sequence_format}) {  ### use convert-seq
	  if (open(SEQ, "| $convert_seq_command -from $sequence_format -to wconsensus -o $sequence_file")) {
	      print SEQ $input{'sequence'};
	      close SEQ;
	  }
      } else {
	  &cgiError("Error: $sequence_format invalid sequence format");
    }
  } else { ### wconsensus format
      if (open SEQ, ">$sequence_file") {
	  print SEQ $input{'sequence'};
	  close SEQ;
      }
  }
  $parameters .= "-f $sequence_file ";
  &DelayedRemoval($sequence_file);

  ### strands ###
  if ($input{'strands'} =~ /ignore/i) {
    $parameters .= "-c0 ";
  } elsif ($input{'strands'} =~ /separate/i) {
    $parameters .= "-c1 ";
  } elsif ($input{'strands'} =~ /single/i) {
    $parameters .= "-c2 ";
  }

  ### pattern length ###
  if (&IsNatural($input{'length'})) {
    $parameters .= "-L $input{'length'} ";
  }

  ### alphabet ###
  $parameters .= "-A $input{'alphabet'} ";

  ### use designated prior frequencies ###
  if ($input{'prior_freq'} eq "on") {
    $parameters .= "-d ";
  }

  ### seed with first sequence and proceed linearly ###
  if ($input{'seed'} eq "on") {
    $parameters .= "-l ";
  } elsif (&IsNatural($input{'repeats'})) {
  ### expected matches ###
      if ($input{'one_per_seq'} eq "on") {
	  $parameters .= "-n $input{'repeats'} ";
      } else {
	  $parameters .= "-N $input{'repeats'} ";
      }
  }

  ### result file
  $result_file = "$TMP/$tmp_file_name.res";
  $error_file  = "$TMP/$tmp_file_name.err";


  ### print the header
  print <<End_Header;
<HEADER>
<TITLE>CONSENSUS result</TITLE>
</HEADER><BODY BGCOLOR="#FFFFFF">
<H3 ALIGN=CENTER>Matrix extraction (consensus) result $input{'set_name'}</H3>
End_Header



  ### execute the command ###
  if ($input{'output'} eq "display") {
  ### Print the result on Web page
#    open ERR, STDERR;
    open RESULT, "$consensus_command $parameters &|";
#    print "$consensus_command $parameters &\n";

#    `($consensus_command $parameters  > $result_file ) >& $error_file`;
#    open RESULT, $result_file;
# print "($consensus_command $parameters  > $result_file ) >& $error_file\n";

    print "<PRE>";
    while (<RESULT>) {
	print;
    }
    close RESULT;

    #### note: this should recuperate the error message but it does not work. 
    #### actually the error file is not created by wwwrun, I don't understand 
    #### why. I still have to fix it.
#    while (<STDERR>) {
#	print "ERROR\t<B>$_</B>\n";
#    }
#    close ERR;

    print "</PRE>";

  } else {
  ### send an e-mail with the result ###
    if ($input{'user_email'} =~ /(\S+\@\S+)/) {
      $address = $1;
      print "<B>Result will be sent to your account: <P>";
      print "$address</B><P>";
      system "$consensus_command $parameters | $mail_command $address &";
    } else {
      if ($input{'user_email'} eq "") {
        print "<B>ERROR: you did not enter your e-mail address<P>";
      } else {
        print "<B>ERROR: the e-mail address you entered is not valid<P>";
        print "$input{'user_email'}</B><P>";      
      }
    } 
  }
  print "<HR SIZE = 3>";
  print &HtmlBot;

  exit(0);
}





