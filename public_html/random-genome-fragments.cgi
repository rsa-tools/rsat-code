#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
#require "cgi-lib.pl";
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
#### redirect error log to a file
BEGIN {
    $ERR_LOG = "/dev/null";
#    $ERR_LOG = "$TMP/RSA_ERROR_LOG.txt";
    use CGI::Carp qw(carpout);
    open (LOG, ">> $ERR_LOG")
	|| die "Unable to redirect log\n";
    carpout(*LOG);
}
require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
$command = "$SCRIPTS/random-genome-fragments";
$tmp_file_name = sprintf "random-genome-fragments.%s", &AlphaDate;
$result_file = "$TMP/$tmp_file_name.res";

### Read the CGI query
$query = new CGI;

### print the header
&RSA_header("Random genome fragments result", "results");

#### update log file ####
&UpdateLogFile;

&ListParameters() if ($ENV{rsat_echo} >= 2);


############################################################
#### read parameters ####
$parameters = "";

############################################################
## Random fragments

#### template file (optional)
($template_sequence_file, $template_sequence_format) = &MultiGetSequenceFile(1, "$TMP/$tmp_file_name"."_template.fa", 0);

## a template file has been given
if ($template_sequence_file ne '') { 
	## calculates the lengths of this file
	my $length_file = "$TMP/$tmp_file_name".".lengths";
	my $seqlength_cmd = "$SCRIPTS/sequence-lengths -v 1 -i ".$template_sequence_file." -o ".$length_file;
	`$seqlength_cmd`;
	
	$parameters .= " -lf $length_file ";
	
} else {
	#### number of fragments
	$frag_nb = $query->param('frag_nb');
	if (&IsNatural($frag_nb)) {
    	$parameters .= " -r $frag_nb ";
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
## Organims

#### organism 

 if ($query->param('org_select')) {
 	
 	## rsat organism
 	if ($query->param('org_select') eq "rsat_org"){
 		
 		unless ($organism = $query->param('organism')) {
    		&FatalError("You should specify an organism");
			}
		if (defined(%{$supported_organism{$organism}})) {
    		$parameters .= " -org $organism ";
		} else {
    	&FatalError("Organism $organism is not supported on this site");
	}
 	
 	## ensembl organism
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
 		
 		## not compatible with non-RSAT organisms
 		if ($query->param('org_select') ne "rsat_org") {
    		&FatalError("Sequence output is only compatible with RSAT organisms. Select a RSAT organism or choose as output format 'genomic coordinates' ");
		} else {
			$parameters .= " -oformat fasta ";
			}	
			
	## return coordinates
 	} elsif ($query->param('outputformat') eq "outputcoord") {
 		if ($query->param('coord_format')) {
 			$parameters .= " -oformat ".$query->param('coord_format');
 		}
 	}
 }

## repeats
if ($query->param('rm') =~ /on/) {
 	$parameters .= " -rm ";
 }

############################################################
## Command

print "<PRE>command: $command $parameters <P>\n</PRE>"  if ($ENV{rsat_echo} >= 1);

### execute the command ###
if ($query->param('output') eq "display") {
    #&PipingWarning();

    ### prepare data for piping
    open RESULT, "$command $parameters |";

    print '<H2>Result</H2>';
    print '<PRE>';
    while (<RESULT>) {
	print $_;
	$genes .= $_;
    }
    print '</PRE>';
    close(RESULT);


    print "<HR SIZE = 3>";

} elsif ($query->param('output') =~ /server/i) {
    &ServerOutput("$command $parameters", $query->param('user_email'));
} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'));
}
print $query->end_html;

exit(0);



