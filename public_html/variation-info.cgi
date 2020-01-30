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
&RSA_header("variation-info result", 'results');


## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

&ListParameters() if ($ENV{rsat_echo} >= 2);

$prefix = "variation-info";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1,0); $tmp_file_name = &ShortFileName($tmp_file_path);
@result_files = ();


################
## Initalize command


$parameters = "";
$command = "$SCRIPTS/variation-info";

################
## Parameters

## Define species 
$organism = $query->param('organism');

if (defined($supported_organism{$organism})) {
    $organism_name = $supported_organism{$organism}->{'name'};

    my ($species, $assembly) = ('','');
 
    if($organism_name =~ /^([A-Za-z]+_[A-Za-z]+_*[A-Za-z]*)\.(\S+?).\d+/){ # as in plants, installed with ensemblgenomes_FTP_client.mk 
        $species = $1;
        $assembly = $2;
    } else {
        @org_name_split=split(" ",$organism_name);
        $species=join("_", $org_name_split[0], $org_name_split[1]);
        $assembly = $org_name_split[2];
    }

    &RSAT::message::Debug("organism: ".$organism, 
			  "organism_name: ".$organism_name,
			  "species: ".$species,
			  "assembly: ".$assembly,
	);

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
    &cgiError("Organism '",
	      $organism,
	      "' is not supported on this web site.");
}

## Get input

my $input_set_file= $tmp_file_path."variation-info_input";
push @result_files ,("input", $input_set_file) ;
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
} elsif ($query->param('input_url') =~ /\S/) {
    my $url = $query->param('input_url');
    &RSAT::message::Info("Fetching variants from URL ".$url) if ($ENV{rsat_echo} >= 1);
    my $var = "";
    if (open SEQ, ">$input_set_file") {
      $var = get($url);
      if ($var =~ /\S/) {
	print SEQ $var;
	close SEQ;
      } else {
	&RSAT::error::FatalError("Input could not  be downloaded from the URL ".$url);
      }
    }

    ## Check sequence file
    my $file_type = `file $input_set_file`;
    if ($file_type =~ "gzip") {
      &RSAT::message::TimeWarn("Uncompressing input file", $input_set_file);
      my $cmd = "mv ".$input_set_file." ".$input_set_file.".gz";
      $cmd .= " ; gunzip ".$input_set_file.".gz";
      &doit($cmd);
    }

    &RSAT::message::Debug("Variants file=", &RSAT::util::hide_RSAT_path($input_set_file), 
			  "<p>File type=", &RSAT::util::hide_RSAT_path($file_type),
	) if ($main::verbose >= 5);


    ### Read sequence from the textarea "input"
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
	&cgiError("You should enter at least one variant ID or genomic regions in the query box.");
    }
}
$parameters .= " -i ".$input_set_file;
&DelayedRemoval($input_set_file);

### Input format
if ($query->param('input_type')) {
    ($input_type) = split " ", $query->param('input_type'); ### take the first word
    $parameters .= " -format ".$input_type;
}

## Output file
$var_file = "$tmp_file_path.varBed";
push @result_files, ("Variations in varBed format", $var_file);
$parameters .= " -o ".$var_file;

## Error log
$err_file = $tmp_file_path."_error_log.txt";
$parameters .= " 2> ".$err_file;
push @result_files, ("Error log (text)",$err_file);

&ReportWebCommand($command." ".$parameters);

#### execute the command #####

		    
if (($query->param('output') =~ /display/i) ||
    ($query->param('output') =~ /server/i)) {

  ################################################################
  ## Run the command
  system($command." ".$parameters);
    
  ## Print the result
  print '<H4>Result</H4>';

  print "<PRE>";
  open RESULT, $result_file;
  while (<RESULT>) {
    print "$_" unless ($query->param('output') =~ /server/i);
  }
  print "</PRE>";
  close RESULT;
    # open RESULT, "$command $parameters |";

    # ### print the result
    # &PipingWarning();

    # ### open the sequence file on the server
    # if (open MIRROR, ">$var_file") {
    # 	$mirror = 1;
    # 	&DelayedRemoval($var_file);
    # }

    # print "<PRE>";
    # while (<RESULT>) {
    # 	print "$_" unless ($query->param('output') =~ /server/i);
    # 	print MIRROR $_ if ($mirror);
    # }
    # print "</PRE>";
    # close RESULT;
    # close MIRROR if ($mirror);

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
	    <FORM METHOD="POST" ACTION="retrieve-variation-seq_form.cgi">
	    <INPUT type="hidden" NAME="variants_file" VALUE="$var_file">
	    <INPUT type="hidden" NAME="input_type" VALUE="varBed">
	<INPUT type="submit" value="retrieve variants sequence">
	</FORM>
	</TD>
 
End_of_form
print '   
</TR>
</TABLE>
</CENTER>';
    
}
