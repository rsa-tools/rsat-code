#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
require "cgi-lib.pl";
require "RSA.lib.pl";
require "RSA_cgi_lib.pl";
$convert_seq_command = "$SCRIPTS/convert-seq";
$tmp_file_name = sprintf "convert-seq.%s", &AlphaDate;

MAIN:
{
  ### open the Web page for the answer
  print &PrintHeader;

  ### Read the content of the form 
  &ReadParse(*input);

  #### read parameters ####
  $parameters = "";

  ##### input format #####
  if ((lc($input{'input_format'}) eq "raw") || 
      (lc($input{'input_format'}) eq "wc") ||
      (lc($input{'input_format'}) eq "wconsensus") ||
      (lc($input{'input_format'}) eq "ig") ||
      (lc($input{'output_format'}) eq "multi") ||
      (lc($input{'input_format'}) eq "fasta")){
      $parameters .= "-from $input{'input_format'} ";
  }

  ##### output format #####
  if ((lc($input{'output_format'}) eq "raw") || 
      (lc($input{'output_format'}) eq "wc") ||
      (lc($input{'output_format'}) eq "wconsensus") ||
      (lc($input{'output_format'}) eq "ig") ||
      (lc($input{'output_format'}) eq "multi") ||
      (lc($input{'output_format'}) eq "fasta")){
      $parameters .= "-to $input{'output_format'} ";
  }

  ##### input sequence file #####
  unless ($input{'sequence'} =~ /\S/) {
      &cgiError("Error: the sequence box should not be empty.<P>Read on-line manual for more information.");
  }
  open INSEQ, ">$TMP/$tmp_file_name";
  print INSEQ $input{'sequence'};
  close INSEQ;
  $parameters .= "-i $TMP/$tmp_file_name ";
  DelayedRemoval("$TMP/$tmp_file_name");


  ##### add reverse-complement #####
  if ($input{'addrc'} eq "yes") {
      $parameters .= "-addrc ";
  }

  ##### line width #####
  if ($input{'line_width'} =~ /\d+/) {
      $parameters .= "-lw $input{'line_width'} ";
  }

  ### print the header
  print <<End_Header;
<HEADER>
<TITLE>RSA-tools - convert sequence result</TITLE>
</HEADER><BODY BGCOLOR="#FFFFFF">
<H3 ALIGN=CENTER><A HREF="$WWW_RSA/RSA_home.cgi">
RSA-tools</A> - convert sequence result</H3>
End_Header


  ### execute the command ###
  if ($input{'output'} eq "display") {
  ### Print the result on Web page
    open RESULT, "$convert_seq_command $parameters  & |";

print "$convert_seq_command $parameters ";

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
      system "$convert_seq_command $parameters | $mail_command $address &";
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


