#!/usr/bin/env perl
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

### Read the CGI query
$query = new CGI;

### print the header
&RSA_header("retrieve-variation-seq result", 'results');


## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

&ListParameters() if ($ENV{rsat_echo} >= 2);

$prefix = "retrieve-variation-seq -v 0";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1,0); $tmp_file_name = &ShortFileName($tmp_file_path);
@result_files = ();


################
## Initalize command


$parameters = "";
$command = "$C_SCRIPTS/retrieve-variation-seq";

################
## Parameters

## Define species 
$organism = $query->param('organism');
if($organism eq ""){ $organism = $query->param('organism_bg')};
$organism = &CheckOrganismAvail($organism);
if(!($organism eq "")){
#if (defined($supported_organism{$organism})) {
    $organism_name = $supported_organism{$organism}->{'name'};


    if ($ENV{plant_server} == 1){

	$parameters .= " -source plants ";

        @org_name_split=split(/\./,$organism_name) ;
        $species=$org_name_split[0] ;
        $species=~ s/ /_/ ;
        $assembly= $org_name_split[1] ;

    } else {



    @org_name_split=split(" ",$organism_name);
    $species=join("_", $org_name_split[0], $org_name_split[1]);
    $assembly =$org_name_split[2];
    }
    
    $parameters .= " -species ".$species; ## Specied ID is the first two parts of the organims ID
    $parameters .= " -assembly ".$assembly; ## Assembly is the third part of species ID
    if (scalar (@org_name_split)>=4){
	if (scalar (@org_name_split)>4){
	    $species_suffix=join("_",@org_name_split[3..$#org_name_split]);
	}else {
	    $species_suffix=$org_name_split[3];
	}
	$parameters .= " -species_suffix ".$species_suffix; ## 
    }
} else {
    &cgiError("Organism '".$query->param('organism')."' is not supported on this web site.");
}

## Get input

unless ($input_set_file=$query->param('variants_file')){
    $input_set_file= $tmp_file_path."retrieve-variation-seq_input";
    if ($query->param('uploaded_file')) {
	$upload_file = $query->param('uploaded_file');
	if ($upload_file =~ /\.gz$/) {
	    $input_set_file .= ".gz";
	}
	$type = $query->uploadInfo($upload_file)->{'Content-Type'};
	open INPUT_FILE, ">". $input_set_file ||
	    &cgiError("Cannot store input file in temporary directory");
	while (<$upload_file>) {
	    print INPUT_FILE;
	}
	close INPUT_FILE;
    } else {
	my $input_var = $query->param('input');
	$input_varc =~ s/\r/\n/g;
	my @input_var = split (/[\n\r]/, $input_var);
	if ($input_var =~ /\S/) {
	    open QUERY, ">".$input_set_file;
	    foreach my $row (@input_var) {
		next unless $row =~ /\S/; ## Skip empty rows
		chomp($row); ## Suppress newline character
		$row =~ s/ +/\t/; ## replace white spaces by a tab for the multiple genomes option
		print QUERY $row, "\n";
	    }
	    close QUERY;
	} else {
	    &cgiError("You should enter at least one variant ID, genomic regions or variant in rsat format in the query box.");
	}
    }
}
push @result_files ,("input", $input_set_file) ;
$parameters .= " -i ".$input_set_file;
&DelayedRemoval($input_set_file);


### Input format
if ($query->param('input_type')) {
    ($input_type) = split " ", $query->param('input_type'); ### take the first word
    $parameters .= " -format ".$input_type;
}


### Lenght of the longest matrix, sequence window around the variant ###
if (&IsInteger($query->param('mml'))) {
    $parameters .= " -mml ".$query->param('mml');
}

&ReportWebCommand($command." ".$parameters);
$var_file = "$tmp_file_path.varSeq";
push @result_files, ("sequences", $var_file);

#### execute the command #####
if (($query->param('output') =~ /display/i) ||
    ($query->param('output') =~ /server/i)) {

    open RESULT, "$command $parameters |";

    ### print the result
    &PipingWarning();

    ### open the sequence file on the server
    if (open MIRROR, ">$var_file") {
	$mirror = 1;
	&DelayedRemoval($var_file);
    }

    print "<PRE>";
    while (<RESULT>) {
	print "$_" unless ($query->param('output') =~ /server/i);
	print MIRROR $_ if ($mirror);
    }
    print "</PRE>";
    close RESULT;
    close MIRROR if ($mirror);

    &PrintURLTable(@result_files);

    ### prepare data for piping
    &PipingForm();

    print "<HR SIZE = 3>";

} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'), $var_file);
}

print $query->end_html;

exit(0);


################
## Send result file to variation-scan form
sub PipingForm {
    ### prepare data for piping
    print <<End_of_form;
    <CENTER>
	<TABLE class="nextstep">
	<TR>
	<TD colspan=2>
	<H3>Next step</H3>
	</TD>
	</TR>
	<TR>
	<TD valign=top>
	<FORM METHOD="POST" ACTION="variation-scan_form.cgi">
        <INPUT type="hidden" NAME='variants_seq_file' VALUE="$var_file">
	<INPUT type="submit" value="variation scan">
	</FORM>
	</TD>
 
End_of_form
print '   
</TR>
</TABLE>
</CENTER>';
}
