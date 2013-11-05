#!/usr/bin/perl
############################################################
#
# $Id: util.pl,v 1.13 2011/02/17 05:07:46 rsat Exp $
#
# Time-stamp: <2003-07-10 11:46:32 jvanheld>
#
############################################################
### util.pl
### utilities for the AMAZE project

$delivery_date =`date +%Y%m%d`;
chomp $delivery_date;
$null = "<NULL>";

################################################################
#### Initialize all config variables to the default value
foreach my $k (keys %main::default) {
    $$k = $default{$k};
}

unless (defined($parsed_dir)) {
    $parsed_dir = "~/parsed_data";
}

################################################################
#### read generic options
$generic_option_message = "	-h	detailed help
	-help	short list of options
	-test	fast parsing of partial data, for debugging
	-outdir output directory
	-parseddir (default $parsed_dir)
		Directory where parsed data must be stored. The output
		directory is calculated by adding the database name
		and the current dte to this directory.
	-v #	verbosity level
	-obj	export data in object format (.obj file)
		which is human-readable (with some patience and
		a good cup of coffee)
	-nocomp skip data compression
	-clean	remove all files from the output directory before
		parsing
";

#	-schema	db schema (default=$main::default{schema}, current=$main::schema)
#	-user	db user (default=$main::default{user}, current=$main::user)
#	-pass	db password (default=$main::default{password}, current=$main::password)
#	-host	db host (default=$main::default{host}, current=$main::host)


sub ReadGenericOptions {	### warn level
    my ($a) = @_;

    #### level of verbosity
    if (($ARGV[$a] eq "-v" ) && 
	($ARGV[$a+1] =~ /^\d+$/)){
	$main::verbose = $ARGV[$a+1];
	
	### output dir
    } elsif ($ARGV[$a] eq "-outdir") {
	$a++;
	$main::dir{output} = $ARGV[$a];
	
	### parsed data dir
    } elsif ($ARGV[$a] eq "-parseddir") {
	$a++;
	$main::parsed_data = $ARGV[$a];
	$dir{output} = join "/", $main::parsed_data, $main::export_subdir, $main::delivery_date;
	
	### clean
    } elsif ($ARGV[$a] eq "-clean") {
	$main::clean = 1;
	
	### quick test
    } elsif ($ARGV[$a] eq "-test") {
	$test = 1;

	### skip data compression
    } elsif ($ARGV[$a] eq "-nocomp") {
	$no_compression = 1;

	### database host
    } elsif ($ARGV[$a] eq "-host") {
	$main::host = $ARGV[$a+1];

	### database schema
    } elsif ($ARGV[$a] eq "-schema") {
	$main::schema = $ARGV[$a+1];

	### database user
    } elsif ($ARGV[$a] eq "-user") {
	$main::user = $ARGV[$a+1];

	### database password
    } elsif ($ARGV[$a] eq "-pass") {
	$main::password = $ARGV[$a+1];

	### export .obj file
    } elsif ($ARGV[$a] eq "-obj") {
	$export{obj} = 1;

	### help
    } elsif (($ARGV[$a] eq "-h") ||
	     ($ARGV[$a] eq "-help")) {
	&PrintHelp();
	exit(0);
    }
}

#### compress parsed data
sub CompressParsedData {
    warn "; compressing parsed data\n" if ($verbose >= 1);
    unless ($no_compression) {
	system "gzip -f $dir{output}/*.tab $dir{output}/*.txt";
	system "gzip -f $dir{output}/*.obj" if ($export{obj});
    }
}


#### check for the existence of a file
sub checkfile {
    my ($file) = @_;
    unless (-e $file) {
	die "Error: file '".$file."' does not exist\n";
    }
}

sub unquote {
    my ($input_string) = @_;
    my $output_string = $input_string;
    $output_string =~ s/^\"(.*)\"$/$1/g;
    $output_string =~ s/^\'(.*)\'$/$1/g;
    warn "input_string", "\t", $input_string, "\tunquoted\t", $output_string, "\n" if ($verbose >=10);
    return $output_string;
}

sub standardize {
    my ($input_string) = @_;
    my $output_string = lc(&trim($input_string));
    $output_string =~ s/ +/\-/g;
    $output_string =~ s/\-+/\-/g;
    #### specific suppression of single quotes because 
    #### for compound names,  single quotes hve a meaning
    $output_string =~ s/^\"//;
    $output_string =~ s/\"$//;
    warn "input_string", "\t", $input_string, "\tstandardized\t", $output_string, "\n" if ($verbose >=10);
    return $output_string;
}


sub intersection {
    my ($list1_ref, $list2_ref) = @_;
    my @result = ();
    foreach my $elem1 (@$list1_ref) {
	foreach my $elem2 (@$list2_ref) {
	    push @result, $elem1 if ($elem1 eq $elem2);
	}
    }
    return @result;
}

# sub IsInteger {
#   my ($query) = @_;
#   if ($query =~ /^[\+\-]{0,1}\d+$/) {
#     return 1;
#   } else {
#     return 0;
#   }
# }

################################################################
#### Print a message in th error log file
sub ErrorMessage {
  my @messages = @_;
  print ERR join ("\n", @messages), "\n";
  warn @messages if ($main::verbose >= 2);
}


################################################################
#### Increment a counter nd return the value with a fixed number of
#### digits (8)
sub next_count {
  my $count = sprintf "%8d", ++$main::object_count;
  $count =~ s/ /0/g;
  return $count;
}


# ################################################################
# #### print the command-line arguments
# sub PrintArguments {
#     my $local_out = $_[0];
# #    unless ($local_out)  {
# #	$local_out = STDOUT;
# #    }
#     my $argument_string = "";

#     foreach my $a (@main::ARGV) {
# 	if (($a =~ /\s+/)  ||
# 	    ($a !~ /\S+/) ||
# 	    ($a =~ /[\(\)\>\<\&]/)) {
# 	    $argument_string .= " '$a'";
# 	} else {
# 	    $argument_string .= " $a";
# 	}
#     }
#     print $local_out $argument_string, "\n" if ($local_out);
#     return $argument_string;
# }

################################################################
#### Default verbose message
sub DefaultVerbose {
    warn join (" ", "; $0 ", &PrintArguments(STERR)), "\n";
    warn "; directories\n";
    while (($key, $value) = each %dir) {
	warn sprintf ";\t%-23s\t%s\n", $key, $value;
    }
    warn "; input files\n";
    while (($key, $value) = each %in_file) {
	warn sprintf ";\t%-23s\t%s\n", $key, $value;
    }
    warn "; output files\n";
    while (($key, $value) = each %outfile) {
	warn sprintf ";\t%-23s\t%s\n", $key, $value;
    }
}

# ################################################################
# #### Return date and time in alphabetical order (year - month - day)
# sub AlphaDate {
#   ### usage : $alpha_date = &AlphaDate;
#   my $date = `date +%Y_%m_%d_%H%M%S`;
#   $date =~ s/\s+$//;
#   return $date;
# }


# ################################################################
# #### Trim the leading and trailing spaces of a string
# sub trim {
#   ### remove leading and trailing spaces from a string
#   ### usage $trimmed_string = &trim($string);
#   local $string = $_[0];
#   $string =~ s/\s*$//;
#   $string =~ s/^\s*//;
#   return $string;
# }

################################################################
#### Trim the leading and trailing spaces of a string
sub my_trim {
  my ($in_string) = @_;
  ### spaces
  $in_string =~ s/^\s+//;
  $in_string =~ s/\s+$//;
  ### double quotes
  $in_string =~ s/^\"//;
  $in_string =~ s/\"$//;
  return $in_string;
}



### split a line into tab-separated fields ###
sub MySplit {
  @fields = ();
  next if (/^;/);        ### skip comment lines
  next unless (/\S/);    ### skip empty lines
  s/\r//g;               ### fix problem with DOS carriage return
  chomp;                 ### remove carriage return 
  @fields = split "\t";  ### split
  #### trim
  foreach my $f (0..$#fields) {
      $fields[$f] = &my_trim($fields[$f]);
  }
  return @fields;
}

sub PrologString {
### usage
### $prolog_string = &PrologString($string);
  local($string) = $_[0];
  $string =~ s/^\'//;
  $string =~ s/\'$//;
  $string =~ s/\'/prime/g;
  $string =~ s/\"//g;
  return "'$string'";
}

return 1;
