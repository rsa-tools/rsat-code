#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
require "cgi-lib.pl";
require "RSA.lib.pl";
$random_seq_command = "$SCRIPTS/random-seq";
&UpdateLogFile;

MAIN:
{
### Read the content of the form 
  &ReadParse(*input);
  
  #### read parameters ####
  $parameters = "";
  $parameters .= "-l $input{'length'} ";
  if ($input{'repet'} > 0) {
    $parameters .= "-r $input{'repet'} ";
  }
  unless ($input{'lw'} eq "") {
    $parameters .= "-lw $input{'lw'} ";
  }
  CheckInputSeqFormat($input{'format'});
  $parameters .= "-format $input{'format'} ";
   
  if ($input{'proba'} eq "alphabet") {
    $parameters .= "-a a:t $input{'ATfreq'} c:g $input{'CGfreq'} ";
  } elsif ($input{'proba'} eq "expfreq") {
      $oligo_length = $input{'oligo_size'};
      $freq_file = "${oligo_length}nt";
    
    if ($input{'seq_type'} eq "genomic") {
      $freq_file .= ".genomic.freq";
    } elsif ($input{'seq_type'} eq "coding") {
      $freq_file .= ".coding.freq";
    } elsif ($input{'seq_type'} eq "non coding") {
      $freq_file .= ".non-coding.freq";
    }
    
    $parameters .= "-expfreq $RSA/data/yeast/oligo-frequencies/$freq_file ";
  } 
  
  ### send the command
#  open RESULT, "$random_seq_command $parameters |";

  ### Print the result
#  print &PrintHeader;
#  print &HtmlTop ("<CENTER>Random sequence generator result</CENTER>");
#  print "<HR SIZE = 3>";

#  print "<PRE>";
#  while (<RESULT>) {
#    print "$_";
#  }
#  print "</PRE>";
#  close RESULT;

#  print "<HR SIZE = 3>";
#  print &HtmlBot;  
#  exit(0);
#}
  ### open the Web page for the answer
  print &PrintHeader;

    print <<End_Header;
<HEADER>
<TITLE>RSA-tools - random sequence result</TITLE>
</HEADER><BODY BGCOLOR="#FFFFFF">
<H3 ALIGN=CENTER><A HREF="$WWW_RSA/RSA_home.cgi">
RSA-tools</A> - random sequence result</H3>
End_Header


  ### execute the command ###
  if ($input{'output'} eq "display") {
  ### Print the result on Web page
    open RESULT, "$random_seq_command $parameters  & |";
    print "<PRE>";
    while (<RESULT>) {
      print "$_";
    }
    print "</PRE>";
    close RESULT;

  } else {
  ### send an e-mail with the result ###
    if ($input{'user_email'} =~ /(\S+\@\S+)/) {
      $address = $1;
      print "<B>Result will be sent to your account: <P>";
      print "$address</B><P>";
      system "$random_seq_command $parameters | $mail_command $address &";
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

