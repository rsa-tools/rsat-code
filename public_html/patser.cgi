#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
require "cgi-lib.pl";
require "RSA.lib.pl";
require "RSA.cgi.lib.pl";
$patser_command = "$BIN/patser";
$matrix_from_transfac_command = "$SCRIPTS/matrix-from-transfac";
$convert_seq_command = "$SCRIPTS/convert-seq";
$features_from_patser_cmd = "$SCRIPTS/features-from-patser";
$add_orf_function_command = "$SCRIPTS/add-orf-function";
$add_yeast_link_command = "$SCRIPTS/add-yeast-link";
$tmp_file_name = sprintf "patser.%s", &AlphaDate;

MAIN:
{
  ### open the Web page for the answer
  print &PrintHeader;


  #### Read the content of the form 
  &ReadParse(*input);

  #### update log file ####
  &UpdateLogFile;

  #### read parameters ####
  $parameters = "-A a:t 0.325 c:g 0.175 ";

  ### matrix ####
  unless ($input{'matrix'} =~ /\S/) { ### empty matrix
    &cgiError("Error: you did not enter the matrix");
  }

  $matrix_file = "$TMP/$tmp_file_name.matrix";

  $matrix_format = lc($input{'matrix_format'});
  if ($matrix_format =~ /transfac/i) {
    open MAT, "| $matrix_from_transfac_command > $matrix_file";
  } elsif ($matrix_format =~ /consensus/i) {
    open MAT, "> $matrix_file";
  } else {
      &cgiError("Error: invalid matrix format");
  }
  print MAT $input{'matrix'};
  close MAT;
  &DelayedRemoval($matrix_file);
  $parameters .= "-m $matrix_file ";

  ### sequence file ###
  unless ($input{'sequence'} =~ /\S/) { ### empty matrix
      &cgiError("Error: you did not enter the sequence");
  }
  $seq_format = lc($input{'seq_format'});
  $sequence_file = "$TMP/$tmp_file_name.seq";
  if (($seq_format ne "wc") && 
      ($seq_format ne "wconsensus")) { 
      if ($accepted_input_seq{$seq_format}) {  ### use convert-seq
	  if (open(SEQ, "| $convert_seq_command -from $seq_format -to wconsensus -o $sequence_file")) {
	      print SEQ $input{'sequence'};
	      close SEQ;
	  }
      } else {
	  &cgiError("Error: $seq_format invalid sequence format");
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
  if ($input{'strands'} =~ /both/i) {
    $parameters .= "-c ";
  }

  ### top value only ###
  if ($input{'return'} =~ /top/i) {
    $parameters .= "-t ";
  }

  ### thresholds ###
  if (&IsReal($input{'lthreshold'})) {
      $parameters .= "-l $input{'lthreshold'} ";
      $parameters .= "-M $input{'lthreshold'} ";
  } 

  if (&IsReal($input{'uthreshold'})) {
      $parameters .= "-u $input{'uthreshold'} ";
  }


  ### parameters for the piping to the feature map ###
  $feature_file =  "$TMP/$tmp_file_name.ft";
  $features_from_patser_cmd .= " -seq $sequence_file";
  $features_from_patser_cmd .= " -o $feature_file";

  print <<End_Header;
<HEADER>
<TITLE>RSA-tools - matrix search (patser) result</TITLE>
</HEADER><BODY BGCOLOR="#FFFFFF">
<H3 ALIGN=CENTER><A HREF="$WWW_RSA/RSA_home.cgi">
RSA-tools</A> - matrix search (patser) result $input{'set_name'}</H3>
End_Header



  ### execute the command ###
  if ($input{'output'} eq "display") {
  ### Print the result on Web page
    open RESULT, "$patser_command $parameters & |";
    open FEATURES, "| $features_from_patser_cmd";

    ### prepare data for piping
    print <<End_of_form;
<CENTER>
<TABLE>
<TR>
  <TD>
    <H4>Next step</H4>
  </TD>
  <TD>
    <FORM METHOD="POST" ACTION="feature-map.cgi">
    <INPUT type="hidden" NAME="feature_file" VALUE="$feature_file">
    <INPUT type="hidden" NAME="format" VALUE="feature-map">
    <INPUT type="hidden" NAME="fill_form" VALUE="on">
    <INPUT type="submit" value="feature map">
    </FORM>
  </TD>
</TR>
</TABLE>
</CENTER>
End_of_form

    print "<PRE>";
    while (<RESULT>) {
	print;
	print FEATURES;
    }
    close FEATURES;
    close RESULT;
    print "</PRE>";
    print "<HR SIZE=3>\n";

  } else {
  ### send an e-mail with the result ###
    if ($input{'user_email'} =~ /(\S+\@\S+)/) {
      $address = $1;
      print "<B>Result will be sent to your account: <P>";
      print "$address</B><P>";
      system "$patser_command $parameters | $mail_command $address &";
    } else {
      if ($input{'user_email'} eq "") {
        print "<B>ERROR: you did not enter your e-mail address<P>";
      } else {
        print "<B>ERROR: the e-mail address you entered is not valid<P>";
        print "$input{'user_email'}</B><P>";      
      }
    } 
    print "<HR SIZE = 3>";
  }
  print &HtmlBot;

  exit(0);
}





