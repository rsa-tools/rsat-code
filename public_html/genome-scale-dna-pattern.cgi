#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
require "$RSA/public_html/genome-scale.lib.pl";

$dna_pattern_command = "$SCRIPTS/dna-pattern";
$add_linenb_command = "$SCRIPTS/add-linenb";
$add_orf_function_command = "$SCRIPTS/add-orf-function";
$link_command = "$SCRIPTS/add-yeast-link -db all ";
$tmp_file_name = sprintf "genome-scale-dna-pattern.%s", &AlphaDate;

### Read the CGI query
$query = new CGI;

### print the header of the result page
&RSA_header("dna-pattern result ".$query->param("title"));

#&ListParameters;


#### update log file ####
&UpdateLogFile;

&ReadRetrieveSeqParams();

################################################################
#
# dna-pattern parameters
#

#### read parameters ####
$parameters_dna_pattern = "";
$parameters_dna_pattern .= " -v";
$parameters_dna_pattern .= " -format ".$query->param('sequence_format');

### pattern file ####
unless ($query->param('patterns') =~ /\S/) {
  &cgiError("The pattern box should not be empty.<P>Read on-line manual for more information.");
}
$pattern_file = "$TMP/$tmp_file_name.pat";
if (open PAT, ">$pattern_file") {
  print PAT $query->param('patterns');
  close PAT;
  &DelayedRemoval($pattern_file);
}
$parameters_dna_pattern .= " -pl $pattern_file";

### return match count ###
if ($query->param('return') =~ /count/i) {
  $parameters_dna_pattern .= " -c";
  if (($query->param('threshold') =~ /^\d+$/) && ($query->param('threshold') > 0)) {
    $parameters_dna_pattern .= " -th ".$query->param('threshold');
  }


### return match count table
} elsif ($query->param('return') =~ /table/i) {
  $parameters_dna_pattern .= " -table";
  ### add a rwo and a column with the totals
  if (lc($query->param('total')) eq "on") {
    $parameters_dna_pattern .= " -total";
  }
  
### return matching positions
} elsif ($query->param('return') =~ /positions/) { 
  ### origin ###
  if ($query->param('origin') =~ /end/i) {
    $parameters_dna_pattern .= " -origin -0";
  }
  
  if ($query->param('flanking') =~ /^\d+$/) {
    $parameters_dna_pattern .= " -N ".$query->param('flanking');
  }
}

### prevent overlapping matches
if (lc($query->param('noov')) eq "on") {
  $parameters_dna_pattern .= " -noov";
}


### strands
if ($query->param('strands') =~ /direct/i) {
  $parameters_dna_pattern .= " -D";
} elsif  ($query->param('strands') =~ /reverse/i) {
  $parameters_dna_pattern .= " -R";
}

### substitutions
if ($query->param('subst') =~ /^\d+$/) {
  $parameters_dna_pattern .= " -subst ".$query->param('subst');
}


$command = "$retrieve_seq_command $retrieve_seq_parameters ";
$command .= "| $dna_pattern_command $parameters_dna_pattern ";
$command .= "| $add_orf_function_command -org $org ";
$command .= "| $add_linenb_command ";
if ($org eq "Saccharomyces_cerevisiae") { #### not yet supported for other organisms
    $command .= "| $link_command  ";
}

### execute the command ###
if ($query->param("output") =~ /display/i) {

    ### execute the command ###
    $result_file = "$TMP/$tmp_file_name.res";
    print "<PRE>$command</PRE>" if ($ECHO);
    open RESULT, "$command & |";
    
    ### Print the result on Web page
    &PrintHtmlTable(RESULT, $result_file,,,,1);
    close RESULT;
    
    $export_genes = `grep -v '^;' $result_file | cut -f 1 | sort -u '`;
    
    &PipingForm () ;
    
    
} else {
    ### send an e-mail with the result ###
    if ($query->param('user_email') =~ /(\S+\@\S+)/) {
	$address = $1;
	print "<B>Result will be sent to your e-mail address: <P>";
	print "$address</B><P>";
	system "$command | $mail_command $address &";
    } else {
	if ($query->param('user_email') eq "") {
	    print "<B>ERROR: you did not enter your e-mail address<P>";
	} else {
	    print "<B>ERROR: the e-mail address you entered is not valid<P>";
	    print "$query->param('user_email')</B><P>";      
	}
    } 
    print "<HR SIZE = 3>";
}

print $query->end_html;


exit(0);

sub PipingForm {
  ### prepare data for piping
  $title = $query->param("title");
  $title =~ s/\"/\'/g;
  $organism = $org;
  $organism =~ s/_/ /g;
## if ($query->param('return') =~ /positions/) {
## if ($org eq "Saccharomyces_cerevisiae") {

  print <<End_of_form;
<HR SIZE = 3>
<CENTER>
<TABLE>

<TR>

<TD>
<H3>Next step</H3>
</TD>

<TD>
<FORM METHOD="POST" ACTION="feature-map_form.cgi">
<INPUT type="hidden" NAME="title" VALUE="$title">
<INPUT type="hidden" NAME="feature_file" VALUE="$result_file">
<INPUT type="hidden" NAME="format" VALUE="dna-pattern">
<INPUT type="hidden" NAME="fill_form" VALUE="on">
<INPUT type="submit" VALUE="feature map">
</FORM>
</TD>

<TD>
</TD>

</TR>


<TR>

<TR>
<TD colpsan=2>
<h3>External servers</h3>
</TD>
</tr>

<tr>
<TD>
<a href="http://www.biologie.ens.fr/fr/genetiqu/puces/publications/ymgv_NARdb2002/index.html" target=_blank>yMGV transcription profiles</a>
</TD><td>
<FORM METHOD="POST" ACTION="http://www.transcriptome.ens.fr/ymgv/list_signatures.php3" target=_blank>
<INPUT type="hidden" NAME="generequest" VALUE="$export_genes">
<INPUT type="submit" VALUE="Send">
</FORM>
</TD>
</tr>

<!--
<tr>
<TD>
<a href="http://www.genome.ad.jp/kegg/kegg2.html#pathway" target=_blank>KEGG pathway coloring</a>
</TD><td>
<FORM METHOD="POST" ACTION="http://www.genome.ad.jp/kegg-bin/search_pathway_multi_www" target=_blank>
<INPUT type="hidden" NAME="org_name" VALUE="$organism">
<INPUT type="hidden" NAME="unclassified" VALUE="$export_genes">
<INPUT type="submit" VALUE="Send">
</FORM>
</TD>
</TR>
-->


</TABLE>
</CENTER>
End_of_form
}
