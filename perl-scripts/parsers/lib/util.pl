#!/usr/bin/perl
############################################################
#
# $Id: util.pl,v 1.3 2003/11/30 07:46:17 jvanheld Exp $
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

################################################################
#### read generic options
$supported_dbms = join "," , keys %supported_dbms;
$generic_option_message ="	-h	detailed help
	-help	short list of options
	-test	fast parsing of partial data, for debugging
	-outdir output directory
	-parseddir (default $parsed_dir)
		Directory where parsed data must be stored. The output
		directory is calculated by adding the database name
		and the current dte to this directory.
	-v #	warn level
		Warn level 1 corresponds to a restricted verbose
		Warn level 2 reports all object instantiations
		Warn level 3 reports failing get_attribute()
	-obj	export data in object format (.obj file)
		which is human-readable (with some patience and
		a good cup of coffee)
	-nocomp skip data compression
	-clean	remove all files from the output directory before
		parsing
	-dbms	database management system (default $default{dbms})
		supported: $supported_dbms
	-db	database schema (default=$main::default{schema}, current=$main::schema)
";


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

	### dbms
    } elsif ($ARGV[$a] eq "-dbms") {
	$main::dbms = $ARGV[$a+1];
	unless ($supported_dbms{$main::dbms}) {
	    die "Error: this dbms is not supported\n";
	}

	### database schema
    } elsif ($ARGV[$a] eq "-schema") {
	$main::schema = $ARGV[$a+1];

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

sub IsInteger {
  my ($query) = @_;
  if ($query =~ /^[\+\-]{0,1}\d+$/) {
    return 1;
  } else {
    return 0;
  }
}

sub ErrorMessage {
  my @messages = @_;
  print ERR @messages;
  warn @messages if ($main::verbose >= 2);
}

sub next_count {
  my $count = sprintf "%8d", ++$main::object_count;
  $count =~ s/ /0/g;
  return $count;
}
sub DefaultVerbose {
  warn "; directories\n";
  while (($key, $value) = each %dir) {
    warn sprintf ";\t%-23s\t%s\n", $key, $value;
  }
  warn "; input files\n";
  while (($key, $value) = each %in_file) {
    warn sprintf ";\t%-23s\t%s\n", $key, $value;
  }
  warn "; output files\n";
  while (($key, $value) = each %out_file) {
    warn sprintf ";\t%-23s\t%s\n", $key, $value;
  }
}

sub AlphaDate {
  ### usage : $alpha_date = &AlphaDate;
  my $date = `date +%Y_%m_%d_%H%M%S`;
  $date =~ s/\s+$//;
  return $date;
}


sub trim {
  ### remove leading and trailing spaces from a string
  ### usage $trimmed_string = &trim($string);
  local $string = $_[0];
  $string =~ s/\s*$//;
  $string =~ s/^\s*//;
  return $string;
}


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
