#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
require "cgi-lib.pl";
require "RSA.lib.pl";
require "RSA_cgi_lib.pl";
$gibbs_command = "nice -n 30 $BIN/gibbs";
$convert_seq_command = "$SCRIPTS/convert-seq";
$tmp_file_name = sprintf "gibbs.%s", &AlphaDate;

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
  $convert_seq_options = "-o $sequence_file -from  $sequence_format -to fasta ";

  if (lc($input{add_rc}) eq "on") {
      $add_rc = 1;
      $convert_seq_options .= "-addrc ";
  }

  if (($sequence_format ne "fasta") || ($add_rc)) {
      if ($accepted_input_seq{$sequence_format}) {  ### use convert-seq
#print "$convert_seq_command  $convert_seq_options";
	  if (open(SEQ, "| $convert_seq_command $convert_seq_options")) {
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
  $parameters .= " $sequence_file ";
  &DelayedRemoval($sequence_file);

  ### pattern length ###
  if (&IsNatural($input{'length'})) {
    $parameters .= " $input{'length'} ";
  }

  ### expected number of matches
  if (&IsNatural($input{$expected})) {
      $paramaters .= "$input{$expected} ";
  }

  ### sequence type
  if (lc($input{seq_type}) eq "dna") {
      $parameters .= "-n ";
  }

  ### inactivate frqgmentation
  unless (lc($input{fragmentation}) eq "on") {
      $parameters .= "-d ";
  }

  ### print the header
  print <<End_Header;
<HEADER>
<TITLE>GIBBS result</TITLE>
</HEADER><BODY BGCOLOR="#FFFFFF">
<H3 ALIGN=CENTER>Matrix extraction (gibbs) result $input{'set_name'}</H3>
End_Header



  ### execute the command ###
  if ($input{'output'} eq "display") {
  ### Print the result on Web page
    open RESULT, "$gibbs_command $parameters & |";

    print "<PRE>";
#print "$gibbs_command $parameters &\n";
    while (<RESULT>) {
	print;
    }
    close RESULT;
    print "</PRE>";

  } else {
  ### send an e-mail with the result ###
    if ($input{'user_email'} =~ /(\S+\@\S+)/) {
      $address = $1;
      print "<B>Result will be sent to your account: <P>";
      print "$address</B><P>";
      system "$gibbs_command $parameters | $mail_command $address &";
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





