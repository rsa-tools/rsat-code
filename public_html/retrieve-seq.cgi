#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
  push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib.pl";
require "RSA.cgi.lib.pl";

$retrieve_seq_command = "$SCRIPTS/retrieve-seq";
$tmp_file_name = sprintf "retrieve-seq.%s", &AlphaDate;

### Read the CGI query
$query = new CGI;

### print the header
&RSA_header("retrieve-seq result");


#### update log file ####
&UpdateLogFile;

#&ListParameters();

#### organism
if (defined($supported_organism{$query->param('organism')})) {
  $organism = $supported_organism{$query->param('organism')}->{'name'};
  $parameters .= " -org ".$query->param('organism');
} else {
    &cgiError("Organism '",
	      $query->param('organism'),
	      "' is not supported on this web site.");
}

### sequence type
if ($query->param('sequence_type')) {
  $parameters .= " -type ".$query->param('sequence_type');
}

### output format ###
$out_format = lc($query->param('format'));
if ($accepted_output_seq{$out_format}) {
  $parameters .= " -format $out_format";
}

### sequence label
$seq_label = lc($query->param('seq_label'));
if (($seq_label =~ /gene/) && 
    ($seq_label =~ /orf/)) {
  $parameters .= " -label orf_gene";
} elsif ($seq_label =~ /gene/) {
  $parameters .= " -label gene";
} elsif ($seq_label =~ /orf/) {
  $parameters .= " -label orf";
} elsif ($seq_label =~ /full/) {
  $parameters .= " -label full";
} else {
  &cgiError("Invalid option for sequence label '$seq_label'");
}

### limits ###
if (&IsInteger($query->param('from'))) {
  $parameters .= " -from ".$query->param('from');
}  
if (&IsInteger($query->param('to'))) {
  $parameters .= " -to ".$query->param('to');
}

### orf overlap ###
unless (lc($query->param('orf_overlap')) eq "on") {
  $parameters .= " -noorf ";
}

#### queries ####
if ($query->param('genes') eq "all") {
  ### take all genes as query
  $parameters .= " -all ";
} elsif ($query->param('uploaded_file')) {
    $upload_file = $query->param('uploaded_file');
    $gene_list_file = "${TMP}/${tmp_file_name}.genes";
    if ($upload_file =~ /\.gz$/) {
	$gene_list_file .= ".gz";
    }
    $type = $query->uploadInfo($upload_file)->{'Content-Type'};
    open SEQ, ">$gene_list_file" ||
	&cgiError("Cannot store gene list file in temporary directory");
    while (<$upload_file>) {
	print SEQ;
    }
    close SEQ;
    $parameters .= " -i $gene_list_file ";

} else {
  $gene_selection = $query->param('gene_selection');
  unless ($gene_selection =~ /\S/) {
    &cgiError("You should enter at least one gene identifier in the query box. <P>Read the on-line manual for more information.");
  }
  $gene_selection .= "\n";			### make sure there is a carriage return at the end
  
  
  @query_lines = split("\n",$gene_selection);
  foreach $line (@query_lines) {
    if ($line =~ /(\S+)/) {
      $parameters .= " -q '$1' "; 
    }
  }
}
  
#### execute the command #####
if ($query->param('output') =~ /display/i) {
  open RESULT, "$retrieve_seq_command $parameters |";
  #    print "<B>\$RSA</B>$RSA<P>";
#      print  "<PRE><B>Command :</B> $retrieve_seq_command $parameters</PRE><P>";
  
  ### open the mirror file ###
  $mirror_file = "$TMP/$tmp_file_name.res";
  if (open MIRROR, ">$mirror_file") {
    $mirror = 1;
    &DelayedRemoval($mirror_file);
  }
  
  ### prepare data for piping
  print <<End_of_form;
<TABLE CELLSPACING=0 CELLPADDING=10 BORDER=0 NOWRAP BGCOLOR= #FFEEDD>

<TR VALIGN="top" ALIGN="center">
    <TD COLSPAN=5 BGCOLOR=	#FFEEDD>
	<FONT SIZE=+1><B>Next step</B></FONT>
    </TD>

</TR>

<TR VALIGN="top" ALIGN="center">

    <TD BGCOLOR=		#FFEEDD>
	<B>Pattern discovery</B><BR>
	(unknown patterns)
    </TD>

    <TD>
	<FORM METHOD="POST" ACTION="oligo-analysis_form.cgi">
	<INPUT type="hidden" NAME="organism" VALUE="$organism">
	<INPUT type="hidden" NAME="from" VALUE="$from">
	<INPUT type="hidden" NAME="to" VALUE="$to">
	<INPUT type="hidden" NAME="sequence_file" VALUE="$mirror_file">
	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
	<INPUT type="submit" value="oligonucleotide analysis">
	</FORM>
    </TD>

    <TD>
	<FORM METHOD="POST" ACTION="dyad-analysis_form.cgi">
	<INPUT type="hidden" NAME="organism" VALUE="$organism">
	<INPUT type="hidden" NAME="from" VALUE="$from">
	<INPUT type="hidden" NAME="to" VALUE="$to">
	<INPUT type="hidden" NAME="sequence_file" VALUE="$mirror_file">
	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
	<INPUT type="submit" value="dyad analysis">
	</FORM>
    </TD>

    <TD>
	<FORM METHOD="POST" ACTION="fill-form.consensus.cgi">
	<INPUT type="hidden" NAME="organism" VALUE="$organism">
	<INPUT type="hidden" NAME="from" VALUE="$from">
	<INPUT type="hidden" NAME="to" VALUE="$to">
	<INPUT type="hidden" NAME="sequence_file" VALUE="$mirror_file">
	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
	<INPUT type="submit" value="consensus">
	</FORM>
    </TD>

    <TD>
	<FORM METHOD="POST" ACTION="gibbs_form.cgi">
	<INPUT type="hidden" NAME="organism" VALUE="$organism">
	<INPUT type="hidden" NAME="from" VALUE="$from">
	<INPUT type="hidden" NAME="to" VALUE="$to">
	<INPUT type="hidden" NAME="sequence_file" VALUE="$mirror_file">
	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
	<INPUT type="submit" value="gibbs sampler">
	</FORM>
    </TD>

</TR>

<TR VALIGN="top" ALIGN="center">

    <TD BGCOLOR=		#FFEEDD>
	<B>Pattern search</B><BR>
	(known patterns)
    </TD>

    <TD>
	<FORM METHOD="POST" ACTION="dna-pattern_form.cgi">
	<INPUT type="hidden" NAME="organism" VALUE="$organism">
	<INPUT type="hidden" NAME="from" VALUE="$from">
	<INPUT type="hidden" NAME="to" VALUE="$to">
	<INPUT type="hidden" NAME="sequence_file" VALUE="$mirror_file">
	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
	<INPUT type="submit" value="dna-pattern (IUPAC)">
	</FORM>
    </TD>

    <TD>
	<FORM METHOD="POST" ACTION="patser_form.cgi">
	<INPUT type="hidden" NAME="organism" VALUE="$organism">
	<INPUT type="hidden" NAME="from" VALUE="$from">
	<INPUT type="hidden" NAME="to" VALUE="$to">
	<INPUT type="hidden" NAME="sequence_file" VALUE="$mirror_file">
	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
	<INPUT type="submit" value="patser (matrices)";
	</FORM>
    </TD>

    <TD>
    &nbsp;
    </TD>

    <TD>
    &nbsp;
    </TD>

</TR>

</TABLE>
End_of_form
  
  ### print the result ### 
  print '<H4>Result</H4>';
  print "<PRE>";
  while (<RESULT>) {
    print "$_";
    print MIRROR $_ if ($mirror);
  }
  print "</PRE>";
  close RESULT;
  
} else {
  ### send an e-mail with the result ###
  if ($query->param('user_email') =~ /(\S+\@\S+)/) {
    $address = $1;
    print "<B>Result will be sent to your account: <P>";
    print "$address</B><P>";
    system "$retrieve_seq_command $parameters | $mail_command $address";
  } else {
    if ($query->param('user_email') eq "") {
      print "<B>ERROR: you did not enter your e-mail address<P>";
    } else {
      print "<B>ERROR: the e-mail address you entered is not valid<P>";
      print "$query->param('user_email')</B><P>";      
    }
  } 
}

print "<HR SIZE = 3>";
print $query->end_html;

exit(0);

