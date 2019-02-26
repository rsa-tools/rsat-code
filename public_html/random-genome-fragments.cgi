#!/usr/bin/env perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
#require "cgi-lib.pl";
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
#### redirect error log to a file
#BEGIN {
#    $ERR_LOG = "/dev/null";
#    $ERR_LOG = &RSAT::util::get_pub_temp()."/RSA_ERROR_LOG.txt";
#    use CGI::Carp qw(carpout);
#    open (LOG, ">> $ERR_LOG")
#	|| die "Unable to redirect log\n";
#    carpout(*LOG);
#}
require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

$command = $SCRIPTS."/random-genome-fragments";
$prefix="random-genome-fragments";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); $tmp_file_name = &ShortFileName($tmp_file_path);

@result_files = ();

### Read the CGI query
$query = new CGI;

### print the header
&RSA_header("Random genome fragments result", "results");


## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

&ListParameters() if ($ENV{rsat_echo} >= 2);


############################################################
#### read parameters ####
$parameters = "";

############################################################
## Random fragments

## Query type
my $query_type = $query->param('fragment_sizes');
if ($query_type eq "template") {

  ## Template format (we cannot take it from MultiGetSequenceFile because
  ## we also support bed format).
  $template_format = $query->param('template_format');
  
  ## template file (optional)
  ($template_file) = &MultiGetSequenceFile(1, $tmp_file_path."_template.".$template_format, 0);
  
  
  ## a template file has been given
  my $length_file = "";
#  if ($template_file) {
  push @result_files, ("Template file ($template_format)",$template_file);
    
#  ## Compute sequence lengths from the template sequence file
#  $length_file = $tmp_file_path.".lengths";
#  push @result_files, ("Sequence lengths",$length_file);

#  my $seqlength_cmd = $SCRIPTS."/sequence-lengths -v 1 -i ".$template_file;
#  $seqlength_cmd .= " -in_format ".$template_format;
#  $seqlength_cmd .= " -o ".$length_file;
#  system($seqlength_cmd);


  ## Add the sequence length file as template for random-genome-fragments
#  $parameters .= " -template_format len -i ".$length_file;
  $parameters .= " -template_format ".$template_format." -i ".$template_file;

} else {
  #### number of fragments
  $frag_nb = $query->param('frag_nb');
  if (&IsNatural($frag_nb)) {
    $parameters .= " -n $frag_nb ";
  } else {
    &FatalError("Fragment number must be a natural number");
  }

  #### length of fragments
  $frag_length = $query->param('frag_length');
  if (&IsNatural($frag_length)) {
    $parameters .= " -l $frag_length ";
  } else {
    &FatalError("Fragment length must be a natural number");
  }
}

############################################################
## Organim
local $organism;
if ($query->param('org_select')) {
  ## RSAT organism
  if ($query->param('org_select') eq "rsat_org"){
    unless ($organism = $query->param('organism')) {
      &FatalError("You should specify an organism");
    }
    $organism = &CheckOrganismAvail($organism);
    #if (%{$supported_organism{$organism}}) {
    if(! ($organism eq "")){
        $parameters .= " -org $organism ";
    } else {
      &FatalError("Organism ".$query->param('organism'). " is not supported on this site");
    }

    ## EnsEMBL organism
  } elsif ($query->param('org_select') eq "ensembl_org") {
    unless ($organism_ens = $query->param('organism_ens')) {
      &FatalError("You should specify an Ensembl organism");
    }
    $parameters .= " -org_ens $organism_ens ";
  }
}


############################################################
## Output
if ($query->param('outputformat')) {

  ## return sequence
  if ($query->param('outputformat') eq "outputseq"){
    $output_format = "fasta";

    ## Sequence output is only valid for RSAT organisms
    if ($query->param('org_select') ne "rsat_org") {
      &FatalError("Sequence output is only compatible with RSAT organisms. Select a RSAT organism or choose as output format 'genomic coordinates' ");
    } else {
      $parameters .= " -return seq ";
    }

    ## Return coordinates
  } elsif ($query->param('outputformat') eq "outputcoord") {
    if ($query->param('coord_format')) {
      $output_format = $query->param('coord_format');
      $parameters .= " -return coord -coord_format ".$output_format;
      $parameters .= " -v 1 ";
    }
  }
}

## Mask repeated sequences
if ($query->param('rm') =~ /on/) {
  $parameters .= " -rm ";
}



## Output file
$result_file = $tmp_file_path."_fragments.".$output_format;
$parameters .= " -o ".$result_file;
#&RSAT::message::Info("result_file", $result_file) if ($echo >= 0);
push @result_files, ("Genome fragments ($output_format)",$result_file);

## Error log
$err_file = $tmp_file_path."_error_log.txt";
$parameters .= " 2> ".$err_file;
push @result_files, ("Error log (text)",$err_file);

############################################################
## Report the command
&ReportWebCommand($command." ".$parameters);

################################################################
## Run the command
# open RESULT, "$command $parameters |";

## open RESULT, "perl /export/space7/rsa-tools/perl-scripts/random-genome-fragments -org Saccharomyces_cerevisiae -n 10 -l 10 |";

if (($query->param('output') =~ /display/i) ||
    ($query->param('output') =~ /server/i)) {
  &PipingWarning();

  ## Run the command
  &doit("$command $parameters"); 


  if ($query->param('output') =~ /display/i) {
      open RESULT, "$result_file"; 

      ## Print the output on the screen
      print '<h4>Result</h4>';
      print '<pre>';
      while (<RESULT>) {
	  print $_;
      }
      print '</pre>';
      close(RESULT);
  }
  # ## Print the result
  # print # '<H4>Result</H4>';
  # 
  # ## Open the sequence file on the server
  # if (open MIRROR, ">$result_file") {
  #   $mirror = 1;
  #   &DelayedRemoval($result_file);
  # }
  # print "<PRE>";
  # while (<RESULT>) {
  #     print $_;
  #   print "$_" unless ($query->param('output') =~ /server/i);
  #   print MIRROR $_ if ($mirror);
  # }
  # print "</PRE>";
  # close RESULT;
  # close MIRROR if ($mirror);

#die(join ("\n", "HELLO\t", $template_file, $length_file, $command." ".$parameters, "\n"));

  ## Print table with links to the result files
  &PrintURLTable(@result_files);

  ## Prepare data for piping
  if ($query->param('outputformat') eq "outputseq"){
    $out_format = "fasta"; ## Fasta is the only supported format, but it is necessary to specify it for the piping form
    &PipingFormForSequence($result_file, "fasta");

  } elsif ($query->param('coord_format') eq "bed") {
    &PipingForm();
  }

  print "<HR SIZE = 3>";

} else {
  &EmailTheResult("$command $parameters", $query->param('user_email'), $result_file);
}

print $query->end_html;

exit(0);

############################################
sub PipingForm {
	my $assembly = `grep Ensembl ${tmp_file_path}.res `;
	$assembly =~ s/.*assembly:(.*)$/$1/;
    ### prepare data for piping
    print <<End_of_form;
	<hr>
	    <table class = "nextstep">
	    <tr><td colspan = 5><h3>next step</h3></td></tr>


	    <tr valign="top" align="center">


	    <td align=center>
	    <form method="POST" action="retrieve-seq-bed_form.cgi">
	    <input type="hidden" name="input_file" value="$result_file">
	    <input type="hidden" name="organism" value="$organism">
	    <input type="submit" value="retrieve-seq-bed">
	    </form>
	    Get sequences (fasta) from genomic coordinates (bed).
	    </td>
	    </td>
	    </tr>

	    <td align=center>
	    <form method="POST" action="fetch-sequences_form.php">
	    <input type="hidden" name="bedfile" value="$result_file">
	    <input type="submit" value="fetch sequences from UCSC">
	    </form>
	    Fetch sequences corresponding to the coordinates
	    </td>
	    </td>
	    </tr>

	    </table>
End_of_form
}



