#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
require "cgi-lib.pl";
require "RSA.lib.pl";
require "RSA_cgi_lib.pl";
$neighbour_orfs_command = "$SCRIPTS/neighbour-orfs";
$add_yeast_link_command = "$SCRIPTS/add-yeast-link";
$tmp_file_name = sprintf "neighbour-orfs.%s", &AlphaDate;

MAIN:
{
  ## open web page for the answer
  print &PrintHeader;

  ### Read the content of the form 
  &ReadParse(*input);

  #### update log file ####
  &UpdateLogFile;
  
  #### read parameters ####
  $parameters = "-v ";

  ### query format ###
  if ($input{'query_format'} =~ /orf/i) {
    $format = "orf";
  } elsif (($input{'query_format'} =~ /genome search/i) || ($input{'query_format'} =~ /dna-pattern/i)){
      $format = "dna-pattern";
  } else {
      $format = "pos";
  }
  $parameters .= "-format $format ";

  ### output format ###
  if ($input{'query_format'} =~ /ORF/) {
    $parameters .= "-ud ";
  } else {
    $parameters .= "-rl ";
  }
    

  #### queries ####
  unless ($input{'query'} =~ /\S/) {
      &cgiError("Error: the query box should not be empty.<P>Read on-line manual for more information.");
  }
  open QUERY, ">$TMP/$tmp_file_name";
  print QUERY $input{'query'};
  close QUERY;
  $parameters .= "-i $TMP/$tmp_file_name ";

  $command = "$neighbour_orfs_command $parameters";

  ### link to database ###
  if ($input{'database_link'} =~ /mips/i) {
      $command .= " | $add_yeast_link_command -db mips ";
  } elsif ($input{'database_link'} =~ /ypd/i) {
      $command .= " | $add_yeast_link_command -db ypd ";
  } elsif (($input{'database_link'} =~ /standford/i) || ($input{'database_link'} =~ /sgd/i)) {
      $command .= " | $add_yeast_link_command -db sgd ";
  }

  ### Print the header
  print <<End_Header;
<HEADER>
<TITLE>RSA-tools - neighbour-orfs result</TITLE>
</HEADER><BODY BGCOLOR="#FFFFFF">
<H3 ALIGN=CENTER><A HREF="$WWW_RSA/RSA_home.cgi">
RSA-tools</A> - neighbour-orfs result $input{'set_name'}</H3>
End_Header



  #### execute the command #####
  if ($input{'output'} eq "display") {

#    $result_file = "$tmp_file_name.res";
    open RESULT, "$command |";
    PrintHtmlTable(RESULT, $result_file);
    close RESULT;

  } else {
  ### send an e-mail with the result ###
    if ($input{'user_email'} =~ /(\S+\@\S+)/) {
      $address = $1;
      print "<B>Result will be sent to your account: <P>";
      print "$address</B><P>";
      system "$neighbour_orfs_command $parameters | $mail_command $address";
    } else {
      if ($input{'user_email'} eq "") {
        print "<B>ERROR: you did not enter your e-mail address<P>";
      } else {
        print "<B>ERROR: the e-mail address you entered is not valid<P>";
        print "$input{'user_email'}</B><P>";      
      }
    } 
  }

  print '<HR SIZE=3>';
  print &HtmlBot;  
}

