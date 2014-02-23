#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
#### redirect error log to a file
BEGIN {
    $ERR_LOG = "/dev/null";
#    $ERR_LOG = &RSAT::util::get_pub_temp()."/RSA_ERROR_LOG.txt";
    use CGI::Carp qw(carpout);
    open (LOG, ">> $ERR_LOG")
	|| die "Unable to redirect log\n";
    carpout(*LOG);
}
require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
require "$ENV{RSAT}/public_html/genome-scale.lib.pl";

### Read the CGI query
$query = new CGI;

### print the header of the result page
&RSA_header("dna-pattern result ".$query->param("title"), "results");

&ListParameters() if ($ENV{rsat_echo} >= 2);

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

@result_files = ();
$dna_pattern_command = "$SCRIPTS/dna-pattern -nolimits";
$add_linenb_command = "$SCRIPTS/add-linenb";
$add_orf_function_command = "$SCRIPTS/add-gene-info -info descr";
$link_command = "$SCRIPTS/add-yeast-link -db all ";
$prefix = "gs-dna-pattern";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); ($tmp_file_dir, $tmp_file_name) = &SplitFileName($tmp_file_path);


################################################################
#
# Organism
#

$org = $query->param("organism");


################################################################
#
# Sequence retrieval parameters
#
#  if ($query->param("sequence_type") =~ /chromosome/) {
#      $seq_format = $supported_organism{$org}->{'seq_format'};
#      $command = "dna-pattern -i ".$supported_organism{$org}->{'genome'};
#  #      #### DEBUGGING
#  #      print "<PRE>";
#  #      print `ls -lt $supported_organism{$org}->{'genome'}`;
#  #      print `cat $supported_organism{$org}->{'genome'}`;
#  #      print "</PRE>";

#  } else {
#    &ReadRetrieveSeqParams();
#    $command = "$retrieve_seq_command $retrieve_seq_parameters | $dna_pattern_command ";
#    $seq_format = "fasta";
#}


################################################################
#
# dna-pattern parameters
#

#### read parameters ####
$seq_format = "fasta";
$parameters_dna_pattern = "";
$parameters_dna_pattern .= " -v ";
$parameters_dna_pattern .= " -format $seq_format ";

### pattern file ####
unless ($query->param('patterns') =~ /\S/) {
  &cgiError("The pattern box should not be empty.<P>Read on-line manual for more information.");
}
$pattern_file = $tmp_file_path.".pat";
push @result_files, ("patterns",$pattern_file);
if (open PAT, ">$pattern_file") {
  print PAT $query->param('patterns');
  close PAT;
  &DelayedRemoval($pattern_file);
}
$parameters_dna_pattern .= " -pl $pattern_file";

### return match count ###
if ($query->param('return') =~ /count/i) {
  $parameters_dna_pattern .= " -return counts";
  if (($query->param('threshold') =~ /^\d+$/) && ($query->param('threshold') > 0)) {
    $parameters_dna_pattern .= " -th ".$query->param('threshold');
  }


### return match count table
} elsif ($query->param('return') =~ /table/i) {
  $parameters_dna_pattern .= " -return table";
  ### add a rwo and a column with the totals
  if (lc($query->param('total')) eq "on") {
    $parameters_dna_pattern .= " -return total";
  }
  
### return matching positions
} elsif ($query->param('return') =~ /positions/) { 
    $parameters_dna_pattern .= " -return sites";
    
    ### origin ###
    if ($query->param('origin') =~ /end/i) {
	$parameters_dna_pattern .= " -origin -0";
    }
    
    if ($query->param('flanking') =~ /^\d+$/) {
	$parameters_dna_pattern .= " -N ".$query->param('flanking');
    }
    
    
    #### match format
    if ($query->param('match_format') eq "fasta") {
	$parameters_dna_pattern .= " -match_format fasta";
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



#### pattern matching parameters
&ReadRetrieveSeqParams();
$command = "$retrieve_seq_command $retrieve_seq_parameters | $dna_pattern_command ";
$command .= " $parameters_dna_pattern ";

################################################################'
#
# Additional information
#
if ($query->param('return') =~ /positions/) {
    $orf_col = 4;
} else {
    $orf_col = 1;
}

unless (($query->param("sequence_type") =~ /chromosome/) ||
	($query->param("match_format") eq "fasta")) {
    $command .= "| $add_orf_function_command -org $org ";
    $command .= "| $add_linenb_command ";

    #### linking to external databases
#    if ( ($query->param("output") =~ /display/i) &&
#	 ($org eq "Saccharomyces_cerevisiae")) { #### not yet supported for other organisms
#	$command .= "| $link_command -col $orf_col  ";
#    }
}

&ReportWebCommand($command);

################################################################
### execute the command ###
if ($query->param("output") =~ /display/i) {


    ### execute the command ###
    $result_file = $tmp_file_path.".dnapat";
    push @result_files, ("result",$result_file);
    open RESULT, "$command & |";

    ### Print the result on Web page
    &PrintHtmlTable(RESULT, $result_file, "", 100000);
    close RESULT;

    my $gene_file = "$result_file.genes";
    system "grep -v '^;' $result_file | cut -f $orf_col | sort -u > $gene_file";
    $export_genes = `cat $result_file.genes`;
    &DelayedRemoval($gene_file);

    &PrintURLTable(@result_files);

    unless ($query->param("match_format" eq "fasta")) {
#    if ($export_genes =~ /\S/) {
	&PipingForm () ;
#    }
    }

    print "<HR SIZE = 3>";
} else {
    &EmailTheResult($command, $query->param('user_email'), $pattern_file);
}

print $query->end_html;


exit(0);



### prepare data for piping
sub PipingForm {
  $title = $query->param("title");
  $title =~ s/\"/\'/g;
  $organism = $org;
  $organism =~ s/_/ /g;
## if ($query->param('return') =~ /positions/) {
## if ($org eq "Saccharomyces_cerevisiae") {


  print <<part1;
<HR SIZE = 3>
<TABLE class = 'nextstep'>
part1


if ($query->param('return') =~ /positions/) {
#### pipe to feature-map
    print <<part2;
<TR>
<TD>
<H3>Next step</H3>
</TD>
</tr>
<tr>
<TD>
<FORM METHOD="POST" ACTION="feature-map_form.cgi">
<INPUT type="hidden" NAME="title" VALUE="$title">
<INPUT type="hidden" NAME="feature_file" VALUE="$result_file">
<INPUT type="hidden" NAME="format" VALUE="dna-pattern">
<INPUT type="hidden" NAME="fill_form" VALUE="on">
<INPUT type="submit" VALUE="feature map">
</FORM>
</TD>
</TR>
part2
    
}
  
  if (($org eq "Saccharomyces_cerevisiae") && 
      !($query->param("sequence_type") =~ /chromosome/) &&
      !($query->param("match_format") eq "fasta") 
      ) {
#### pipe to external servers
      print <<part3

<tr>
<td><h3>External servers</h3></td>
</tr>

<!--
<tr>
<TD>
<a href="http://www.biologie.ens.fr/fr/genetiqu/puces/publications/ymgv_NARdb2002/index.html" target=_blank>yMGV transcription profiles</a>
</TD></tr><tr><td align = 'left'>
<FORM METHOD="POST" ACTION="http://www.transcriptome.ens.fr/ymgv/list_signatures.php3" target=_blank>
<INPUT type="hidden" NAME="generequest" VALUE="$export_genes">
<INPUT type="submit" VALUE="Send">
</FORM>
</TD>
</tr>
-->

<tr>
<TD>
<a href="http://www.genome.ad.jp/kegg/kegg2.html#pathway" target=_blank>KEGG pathway coloring</a>
</TD>
</tr><tr><td>
<FORM METHOD="GET" target=_blank ACTION='http://www.genome.ad.jp/kegg-bin/search_pathway_multi_www'>
<INPUT type="hidden" NAME=org_name VALUE=sce>
<INPUT type="hidden" NAME=unclassified VALUE="$export_genes">
<INPUT type="submit" VALUE="Send">
</FORM>
</TD>
</TR>

part3
}
  
  print <<End_of_form;
</TABLE>
End_of_form


}
