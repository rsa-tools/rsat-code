#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
  push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

$oligo_analysis_command = "$SCRIPTS/oligo-analysis";
$convert_seq_command = "$SCRIPTS/convert-seq";
$purge_sequence_command = "$SCRIPTS/purge-sequence";
$tmp_file_name = sprintf "oligo-analysis.%s", &AlphaDate;

### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("oligo-analysis result");
#&ListParameters;

#### update log file ####
&UpdateLogFile;

#### read parameters ####
$parameters = "";
$parameters .= " -sort";

#### purge sequence option
$purge = $query->param('purge');

### sequence file
($sequence_file,$sequence_format) = &GetSequenceFile();
if ($purge) {
    $command= "$convert_seq_command -i $sequence_file -from $sequence_format -to fasta | $purge_sequence_command | $oligo_analysis_command -format fasta ";
} else {
    $command= "$oligo_analysis_command -i $sequence_file  ";
}


### fields to return
$return_fields = "";

### occurrences
if ($query->param('occ')) {
  $return_fields .= "occ,";
  ### threshold on occurrences
  if ($query->param('occurrence_threshold') =~ /^\d+$/) {
    $parameters .= " -tho ".$query->param('occurrence_threshold');
  }
} 

### frequencies
if ($query->param('freq')) {
  $return_fields .= "freq,";
} 

### matching sequences
if ($query->param('mseq')) {
  $return_fields .= "mseq,";
  ### threshold on matching sequences
  if ($query->param('ms_threshold') =~ /^\d+$/) {
    $parameters .= " -thms ".$query->param('ms_threshold');
  }  
} 

### observed/expected ratio
if ($query->param('ratio')) {
  $return_fields .= "ratio,";
} 

### rank
if ($query->param('rank')) {
  $return_fields .= "rank,";
} 

### z-score
if ($query->param('zscore')) {
  $return_fields .= "zscore,";
} 

### binomial probabilities
if ($query->param('proba')) {
  $return_fields .= "proba,";
  ### threshold on probabilities
  if ($query->param('proba_occ_threshold') =~ /^[\d\.\-+e]+$/i) {
    $parameters .= " -thpo ".$query->param('proba_occ_threshold');
  }
  ### threshold on significance
  if ($query->param('occ_significance_threshold') =~ /^-{0,1}[\d\.\-+e]+$/i) {
    $parameters .= " -thosig ".$query->param('occ_significance_threshold');
  }
} 

### positions
if ($query->param('pos')) {
  $return_fields .= "pos,";
} 

$return_fields =~ s/,$//;

if ($return_fields eq "") {
  &cgiError("You should select at least one option in the \"Return\" box.");
} else {
  $parameters .= " -return $return_fields";
}

### single strand search
if ($query->param('strand') =~ /single/) {
  $parameters .= " -1str";
} else {
  $parameters .= " -2str";
}

### prevent overlapping matches of the same pattern
unless ($query->param('grouprc')) {
  $parameters .= " -nogrouprc";
} 

### prevent overlapping matches of the same pattern
if ($query->param('noov')) {
  $parameters .= " -noov";
} 

### verbose
$parameters .= " -v";

#### sequence type
$parameters .= " -seqtype ".$query->param("sequence_type");

#### oligo size ####
if ($query->param('oligo_size') =~ /\d/) {
    $oligo_length = $query->param('oligo_size') ;
} 
$parameters .= " -l $oligo_length";

#### expected frequency estimation ####
if ($query->param('freq_estimate') =~ /oligo freq.* non-coding regions/i) {
    if (($query->param("proba")) ||
	($query->param("ratio")) ||
	($query->param("zscore"))
	) {
	
	### check organism
	unless ($organism = $query->param('organism')) {
	    &cgiError("You should specify an organism to use non-coding frequency calibration");
	}
	unless (defined(%{$supported_organism{$organism}})) {
	    &cgiError("Organism $org is not supported on this site");
	}
	$freq_option = " -ncf -org $organism";
    }
} elsif ($query->param('freq_estimate') =~ /residue frequenc/i) {
  $freq_option = " -a input";
} elsif ($query->param('freq_estimate') =~ /markov/i) {
  $freq_option = " -markov";
  if (&IsNatural($query->param('markov_order'))) {
      $freq_option .= " ".$query->param('markov_order');
  }
} elsif ($query->param('freq_estimate') =~ /lexicon/i) {
  $freq_option = " -lexicon";
} else {
  $freq_option = "";
} 
$parameters .= "$freq_option";

#### pseudo weight
if (&IsReal($query->param('pseudo_weight'))) {
    my $pseudo = $query->param('pseudo_weight');
#    if ($pseudo > 1) {
#	&cgiError("Pseudo-weight must be <= 1.");
#    } elsif ($pseudo < 0) {
#	&cgiError("Pseudo-weight must be >= 0.");
#    } 
    $parameters .= " -pseudo $pseudo";
}

#### neighborhood ####
if ($query->param('neighborhood') =~ /N at one position/i) {
  $parameters .= " -oneN";
} elsif ($query->param('neighborhood') =~ /one degenerated position/i) {
  $parameters .= " -onedeg ";
}


print "<PRE>command: $command $parameters<P>\n</PRE>" if $ECHO;

if ($query->param('output') =~ /display/i) {

    &PipingWarning();
    
    ### execute the command ###
    $result_file = "$TMP/$tmp_file_name.res";
    open RESULT, "$command $parameters |";
    
    ### Print result on the web page
    print '<H2>Result</H2>';
    &PrintHtmlTable(RESULT, $result_file, true);
    close(RESULT);
    
    #### oligonucleotide assembly ####
    if (&IsReal($query->param('occ_significance_threshold'))) {
	$fragment_assembly_command = "$SCRIPTS/pattern-assembly -v";
	if ($query->param('strand') =~ /single/) {
	    $fragment_assembly_command .= " -1str";
	} else {
	    $fragment_assembly_command .= " -2str";
	}

	unless ($ENV{RSA_ERROR}) {
	    print "<H2>Fragment assembly</H2>\n";
	    open CLUSTERS, "$fragment_assembly_command -i $result_file |";
	    print "<PRE>\n";
	    while (<CLUSTERS>) {
		print;
	    }
	    print "</PRE>\n";
	    close(CLUSTERS);
	}
    }

    &PipingForm();
  
} else {
  #### send e-mail with the result
  if ($query->param('user_email') =~ /(\S+\@\S+)/) {
    $address = $1;
    print "<B>Result will be sent to your account: <P>";
    print "$address</B><P>";
    system "$command $parameters | $mail_command $address &"; 
  } else {
    if ($query->param('user_email') eq "") {
      &cgiError("You did not enter your e-mail address");
    } else {
      &cgiError("The e-mail address you entered is not valid");
      print $query->param('user_email')."</B><P>";      
    }
  }
  print '<HR SIZE=3>';
}

print $query->end_html;

exit(0);


sub PipingForm {
  ### prepare data for piping
  $title = $query->param('title');
  $title =~ s/\"/\'/g;
  print <<End_of_form;
<HR SIZE = 3>
<TABLE>
<TR>
<TD>
<H3>Next step</H3>
</TD>
<TD>
<FORM METHOD="POST" ACTION="dna-pattern_form.cgi">
<INPUT type="hidden" NAME="title" VALUE="$title">
<INPUT type="hidden" NAME="pattern_file" VALUE="$result_file">
<INPUT type="hidden" NAME="sequence_file" VALUE="$sequence_file">
<INPUT type="hidden" NAME="sequence_format" VALUE="$sequence_format">
<INPUT type="submit" value="pattern matching (dna-pattern)">
</FORM>
</TD>
</TR>
</TABLE>
End_of_form

}



