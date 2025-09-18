################################################################
##
## A class for util handling
##
package RSAT::util;

use POSIX;
use File::Temp qw/ tempdir /;
use RSAT::GenericObject;
use RSAT::message;
use RSAT::error;
#use RSAT::email;
@ISA = qw( RSAT::GenericObject );

=pod

=head1 NAME

    RSAT::util

=head1 DESCRIPTION

This method class rassembles various utilities for RSAT.

=head1 METHODS

=cut


################################################################

=pod

=item IsNatural

Returns 1 if the first argument is a natural number, 0 otherwise

Usage: &RSAT::util::IsNatural($x);

=cut

sub IsNatural {
    if ($_[0] =~ /^\s*\d+\s*$/) {
	return 1;
    } else {
	return 0;
    }
}




################################################################

=pod

=item IsInteger

Returns 1 if the first argument is an integer number, 0 otherwise

Usage: &RSAT::util::IsInteger($x);

=cut
sub IsInteger {
    if ($_[0] =~ /^ *[\+\-]{0,1}\d+ *$/) {
	return 1;
    } else {
	return 0;
    }
}

################################################################

=pod

=item IsReal

Returns 1 if the first argument is a real number, 0 otherwise

Usage: &RSAT::util::IsReal($x);

=cut
sub IsReal {
    my ($value) = @_;
    return if ($#_ < 0);
    return 0 unless defined($value);

    if (($value =~ /^ *[\+\-]{0,1}\d+(\.\d*){0,1} *$/) ||
        ($value =~ /^ *[\+\-]{0,1}\d+(\.\d*){0,1}e[\+\-]{0,1}\d+ *$/i)) {
	return 1;
    } else {
	return 0;
    }
}

################################################################

=pod

=item round

Returns the closest integer value from a real number.

Usage : my $i = &RSAT::util::round($r);

=cut
sub round {
    my ($my_real) = @_;
    my $my_int = sprintf "%.0f", $my_real;
#     if (abs(($my_real - int($my_real))) < 0.5) {
# 	return int($my_real);
#     } elsif ($my_real < 0) {
# 	return int($my_real) - 1;
#     } else {
# 	return int($my_real) + 1;
#     }
    return($my_int);
}


################################################################

=pod

=item number_with_zeros

Usage:

my $number_string = &RSAT::util::number_with_zeros($value, digits=>$digits);

my $number_string = &RSAT::util::number_with_zeros($value, maxval=>$maxval);

Print a Natural number with a fixed number of digits. If neccessary, append
leading zeros. This ensures consistency between alphabetical and numerical
order. It can be useful for naming a succession of files.

=cut

sub number_with_zeros {
  my ($value, %args) = @_;
  my $digits;

  if (defined($args{digits})) {
    $digits = $args{digits};
  } elsif (defined($args{maxval})) {
    $digits = POSIX::ceil(log($args{maxval}+1)/log(10));
  } else {
    &RSAT::error::FatalError("&RSAT::util::number_with_zeros()", "arguments must include either digist=>\$digits or maxval=>\$maxval")
  }
#    &RSAT::message::Debug("&RSAT::util::number_with_zeros()", $digits) if ($main::verbose >= 10);

  &RSAT::error::FatalError("&RSAT::util::number_with_zeros()", $digits, "Invalid value for digit, should be Natural.")
      unless (&IsNatural($digits));
  my $digits_number = sprintf "%${digits}s", $value;
  $digits_number =~ s/ /0/g;
  return($digits_number);
}


=pod

=item trim

Suppress space characters at the start and the end of a string

Usage: $trimmed = &RSAT::util::trim($to_trim);

=cut
sub trim {
  my ($string) = @_;
  if ($string) {
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
  }
  return $string;
}



################################################################

=pod

=item AlphaDate

Returns a date in the inverted form, so that the choronological and
alphabetical sorting ar identical.

Usage: $alpha_date = &RSAT::util::AlphaDate();

=cut
sub AlphaDate {
  my ($sec, $min, $hour,$day,$month,$year) = localtime(time());
  my $my_date = sprintf("%02d-%02d-%02d.%02d%02d%02d", 1900+$year,$month+1,$day,$hour, $min, $sec);
  return $my_date;
}


################################################################

=pod

=item Date

Returns a date in human-readable way.

Usage: my $date = &RSAT::util::Date();

=cut
sub Date {
    my $my_date = "";
    my ($sec, $min, $hour,$day,$month,$year) = localtime(time());
    $my_date = sprintf("%02d/%02d/%04d", $day, $month+1, 1900+$year);
#    $my_date = `date +%Y_%m_%d.%H%M%S`;
    chomp $my_date;
    return $my_date;
}

################################################################

=pod

=item DateAndTime

Returns a date and time in human-readable way.

Usage: my $date_time_string = &RSAT::util::DateAndTime();

=cut
sub DateAndTime {
    my $my_date = "";
    my ($sec, $min, $hour,$day,$month,$year) = localtime(time());
    $my_date = sprintf("%02d/%02d/%04d %02d:%02d", $day, $month+1, 1900+$year, $hour, $min);
#    $my_date = `date +%Y_%m_%d.%H%M%S`;
    chomp $my_date;
    return $my_date;
}


################################################################
## Update the start log file
## Usage:
##  local $start_time = &StartScript();
sub StartScript {
  my $start_time = &AlphaDate();

  ## If specified in the server configuration, report start time
  ## in a specific log file.
  if ($ENV{start_time}) {
    my $script_name = &RSAT::util::ShortFileName($0) || 'undefined';
    my $command = join (" ", $script_name, @ARGV);
    my $login = getlogin || getpwuid($<) || "Kilroy";

    ## Get hostname
    my $hostname = `hostname`;
    chomp($hostname);

    ## Get remote address (only for Web queries)
    my $remote_addr = "";
    if (defined($ENV{REMOTE_ADDR})) {
      $remote_addr = $ENV{REMOTE_ADDR};
    }

    &RSAT::message::TimeWarn("Updating start script log file", $main::start_time_log_file)
      if ($main::verbose >= 4);

    ## Write header of the exec time log file if required
    unless (-e $main::start_time_log_file) {
	if (open LOG, ">".$main::start_time_log_file) {
	    print LOG join ("\t",
			    "#start_date.time",
			    "hostname",
			    "PID",
			    "username",
			    "script_name",
			    "command",
			    "remote_addr",
		), "\n";
	    close LOG;
	    chmod 0666, $main::start_time_log_file;
	}
    }

    if (open LOG, ">>".$main::start_time_log_file) {
	print LOG join ("\t",
			$start_time,
			$$,
			$login,
			$script_name,
			$command,
			$remote_addr,
	    ), "\n";
	close LOG;
	chmod 0666, $main::start_time_log_file;
    }

  }
  return($start_time);
}

################################################################

=pod

=item B<ReportTime()>

Report execution time.

Usage:
  print $out &RSAT::util::ReportExecutionTime($start_time) if ($main::verbose >= 1);

=cut
sub ReportExecutionTime {
  my ($start_time) = @_;
  my $time_report;
  my $done_time = &AlphaDate();
  my $elapsed = times;
  my $hostname = `hostname`;
  chomp($hostname);

  ## Report the execution time string only if verbosity >= 1.
  $time_report  = "; Host name\t".$hostname."\n";
  $time_report .= "; Job started\t".$start_time."\n";
  $time_report .= "; Job done\t".$done_time."\n";
  $time_report .= "; Seconds\t".$elapsed."\n";
  my ($user,$system,$cuser,$csystem) = times;
  $time_report .= ";\tuser\t".$user."\n";
  $time_report .= ";\tsystem\t".$system."\n";
  $time_report .= ";\tcuser\t".$cuser."\n";
  $time_report .= ";\tcsystem\t".$csystem."\n";

  ## If specified in the server configuration, report task + execution
  ## time in log file.
  &RSAT::server::UpdateExecTimeLogFile($start_time, $done_time, $elapsed)
    if $ENV{exec_time};

  return($time_report);
}

################################################################

=pod

=item SplitFileName

Splits a full file name into directory and short file name

Usage: my ($dir, $short_file_name) = &RSAT::util::SplitFileName($full_file_name);

=cut
sub SplitFileName {
    my ($full_file_name) = @_;
    my $short_file_name = $full_file_name;
    my $dir = "";
    if ($full_file_name =~ /[^(\/)]+$/) {
      $dir = $`;
      $short_file_name = $&;
      $dir =~ s|/+|/|g; #### suppress multiple /
      $dir =~ s|/$||; #### suppress the trailing /
    }
    return ($dir, $short_file_name);
}

################################################################

=pod

=item ShortFileName

Removes the directory from a full file name, and returns the short
file name.

Usage: $file_name = &RSAT::util::ShortFileName($full_file_name);

=cut
sub ShortFileName {
  my $full_file_name = $_[0];
  my ($dir, $short_file_name) = &SplitFileName($full_file_name);
  return $short_file_name;
}

################################################################

=pod

=item ExtractPath

Returns the path from a file name.

Usage: my $path = &RSAT::util::ExtractPath($long_file_name);

=cut
sub ExtractPath {
    my $l_file_name = $_[0];
    my $l_path = "";
    while ($l_file_name =~ /\/(.*)/) {
	$l_path .= "$`\/";
	$l_file_name = $1;
    }
    return $l_path;
}


=pod

=item RelativePath

Return a relative path between a referring file (e.g. a HTML report
summarizing the results and images returned by a program) and a
referred file (e.g. an image that has to be displayed in the HTML
file).

Usage:
  my ($link, $shared_path) = &RSAT::util::RelativePath($referring_file, $referred_file);

=cut

sub RelativePath {
  my ($referring_file, $referred_file) = @_;

  $referring_file =~ s|/+|/|g;
  my ($referring_dir, $referring_basename) = &RSAT::util::SplitFileName($referring_file);
  $referring_dir =~ s/\/+$//;
  my @referring_path = split /\/+/, $referring_dir;

  $referred_file =~ s|/+|/|g;
  my ($referred_dir, $referred_basename) = &RSAT::util::SplitFileName($referred_file);
  $referred_dir =~ s/\/+$//;
  my @referred_path = split /\/+/, $referred_dir;

  my @shared_path = ();
  my $shared_path = "";
  my $link = "";

  #   &RSAT::message::Debug("&RSAT::util::RelativePath()",
  # 			"\n\treferring_file", $referring_file,
  # 			"\n\treferring_dir", $referring_dir, join("; ", @referring_path),
  # 			"\n\treferred_file", $referred_file,
  # 			"\n\treferred_dir", $referred_dir, join("; ", @referred_path))
  #     if ($main::verbose >= 10);


  ## Particular cases
  if ($referring_file eq $referred_dir) {
    $link = $referred_basename;
    $shared_path = $referring_dir;

  } elsif ($referred_file =~ /${referring_file}\/(.*)/) {
    ## The referring "file" is a directory in the path of the referred file
    $link = $1 || ".";
    $shared_path = $referring_file;


  } elsif (($referred_file =~ /^\//)
	   && ($referring_file !~ /^\//)) {
    ## Referred file path is absolute
    $link = $referred_file;
    $shared_path = "";

  } else {

    while ((scalar(@referring_path) > 0) && (scalar(@referred_path) > 0)) {
      if ($referring_path[0] eq $referred_path[0]) {
	my $shared_folder = shift(@referring_path);
	shift(@referred_path);
	push @shared_path, $shared_folder;
      } else {
	last;
      }
    }

    my $up_levels = 0;
    $shared_path = join "/", @shared_path;
    if ((defined($referring_path[0])) && ($referring_path[0] eq $referred_basename)) {
      ## The referred path is the folder of the referring file
      $link = ".";
    } else {
      $up_levels = scalar(@referring_path);
      $link = "../" x $up_levels;
      $link .= join ("/", @referred_path, $referred_basename);
    }

    #  &RSAT::message::Debug("Computed relative path", $referring_file, $referred_file,
    #			$shared_path, $link,
    #			)
    #    if ($main::verbose >= 10);

  }
  if (wantarray) {
    return ($link, $shared_path);
  } else {
    return ($link);
  }
}

=pod

=item hide_RSAT_path()

Hide the full RSAT path to avoid exposing it on the Web server
interface. The full path is replaced by $RSAT.

# C'est une idee fantastique, ca n'a absolument pas bousille tout
# NeAT... un grand merci au concepteur genial de ce 'TRUC' qui n'a
# fait aucun test après l'implémentation, et qui m'a fait perdre au
# moins une matinee en debugage!!!!!!  Sylvain

=cut
sub hide_RSAT_path {
  my ($string) = @_;
  $string =~ s|$ENV{RSAT}\/+|\$RSAT\/|g; ## Avoid providing the full path of the RSAT package
  return $string;
}

=pod

=item rsat_path_to_url()

Convert the full path of a document in ${RSAT}/public_html to an URL.

=cut
sub rsat_path_to_url {
  my ($string) = @_;
  $string =~ s|$ENV{RSAT}\/public_html\/|$ENV{rsat_www}\/|g;
  return $string;
}


################################################################

=pod

=item CheckOutDir()

Check for the existence of an output directory, and, if it does not
exist, create it.

Usage: &RSAT::util::CheckOutDir($my_dir, $umask, $chmod);

=cut
sub CheckOutDir {
  my ($output_dir, $umask, $chmod) = @_;

  unless ($output_dir) {
    &RSAT::message::Warning("&RSAT::util::CheckOutDir()", "No directory has been specified.", "Command ignored.") if ($main::verbose >= 1);
    return();
  }


  ## Specify a mask for the new directory
  $umask = 0002 unless ($umask);
  $chmod = '0775' unless ($chmod);
  umask($umask);
  if ($main::verbose >= 4) {
    my $wd = $ENV{PWD}; #`pwd`;
    &RSAT::message::Info("Current directory", $wd);
  }


  if ($output_dir) {
    if (-d $output_dir) {
      &RSAT::message::Warning("&RSAT::util::CheckOutDir()", "Directory $output_dir already exists") if ($main::verbose >= 4);
      return();
    }


    &RSAT::message::Info("Creating directory", $output_dir) if ($main::verbose >= 3);
    mkdir ($output_dir, $chmod);


    unless (-d $output_dir) {
      &RSAT::message::Info("Creating directory with all parents", $output_dir) if ($main::verbose >= 3);
      system("mkdir -p $output_dir"); ## create output directory with all parents
    }
    chmod $chmod, $output_dir; ## Not sure the $chmod argument works with mkdir)

    unless (-d $output_dir) {
      &RSAT::error::FatalError("Could not create output directory $output_dir");
    }


  } else {
    $output_dir = ".";
  }


  ## Change access mode if required
  if ((defined($chmod)) && ($chmod =~ /\d{3}/)) {
    &RSAT::message::Info("Changing access mode", $chmod, $output_dir) if ($main::verbose >= 5);
    system("chmod ".$chmod." ".$output_dir);
  }
}



=pod

=item B<OutFileName>

Define one output file name by concatenating arguments.

The variable $main::dir{output} and $main::param{prefix} must have
been defined beforehand.

Usage:
 $outfile{key} = &OutFileName($subdir, $extension, @name_elements);

=cut
sub OutFileName {
  my ($subdir, $extension, @name_elements) = @_;
  my $dir = $main::dir{output};
  if ($subdir) {
    $dir = &CheckSubDir($subdir);
  }
  my $out_file_name = $dir;
  if (scalar(@name_elements) > 0) {
    $out_file_name .= "/";
    $out_file_name .= join ("_", $main::param{prefix}, @name_elements);
    if ($extension) {
      $out_file_name .= $extension;
    }
  }
  if (wantarray) {
    return ($out_file_name, $dir);
  } else {
    return ($out_file_name);
  }
}

=pod

Define the name of an output subdirectory by concatenating arguments.

A variable $main::dir{output} must have been defined beforehand.

Usage:
   $dir = &CheckSubDir($subdir);

=cut
sub CheckSubDir {
  my ($subdir) = @_;
  my $dir = $main::dir{output};
  $dir = join ("/", $main::dir{output}, $subdir);
  &RSAT::util::CheckOutDir($dir,"",777);
  return($dir);
}


=pod

=item OpenInputFile

Open an input stream. If a file name is specified, open this file. If
not, read input from the STDIN. If the input file has extension .gz,
or .Z, it is decopmpressed on the flight.

Usage
   ($in,$input_dir) = &RSAT::util::OpenInputFile($filename);
returns a file handle

=cut
sub OpenInputFile {
  my ($filename) = @_;
  &RSAT::message::Debug("Opening input stream", $filename) if ($main::verbose >= 10);
    my $input_stream;
    my $input_dir = "";
    my $short_file_name = "";
    if ($filename) {
	($input_dir, $short_file_name) = &RSAT::util::SplitFileName($filename);
	unless (-e $filename) {
	    if (-e "${filename}.gz") {
		$filename .= ".gz";
	    } elsif (-e "${filename}.Z") {
		$filename .= ".Z";
	    } else {
		&RSAT::error::FatalError("OpenInputFile: File $filename does not exist.");
	    }
	}
	if ($filename =~ /\.gz$/) { ### gzip file -> decompress it on the fly
	    warn "; Uncompressing .gz file\n" if ($main::verbose >= 4);
	    $filename = "gunzip -c $filename |";
				  } elsif ($filename =~ /\.Z$/) { ### gzip file -> decompress it on the fly
				      warn "; Uncompressing .Z file\n" if ($main::verbose >= 4);
				      $filename = "uncompress -c $filename |";
								}
	open (my $input_fh, $filename) ||
	    &RSAT::error::FatalError ("Cannot read input file ".$filename);
	return ($input_fh, $input_dir)
    } else {
	my $input_fh = *STDIN;
	return ($input_fh, $input_dir);
    }

#    &RSAT::message::Debug("Input stream has been open", $filename, $input_dir, $input_stream) if ($main::verbose >= 10);

   return ($input_stream, $input_dir);
}

################################################################

=pod

=item OpenOutputFile

Open an output stream. If a file name is specified, open this file. If
not, write output to the STDOUT. If the file name has extension .gz,
or .Z, it is copmpressed on the flight.

Usage
   $out = &RSAT::util::OpenOutputFile($filename);
returns a file handle

=cut
sub OpenOutputFile {
  ### usage $output_stream = &OpenOutputFile($filename);
  my ($filename) = @_;
  my $to_open = "";
  my $fh;

  if ($filename) {
    if ($filename =~ /\.gz$/) { ### gzip file -> decompress it on the fly
      $to_open = "| gzip -c > $filename";
    } else {
      $to_open = ">$filename";
    }
#    warn("HELLO\t", $to_open, "\n");
    open($fh, $to_open) ||
	&RSAT::error::FatalError ("Cannot write output file ".$to_open);
  } else {
    $fh = *STDOUT;
  }
  return ($fh);
}

################################################################

=pod

=item AppendOutputFile

Append to an output stream. If a file name is specified, open and
append to this file. If not, write output to the STDOUT. If the
file name has extension .gz,or .Z, it is compressed on the flight.

Usage
   $out = &RSAT::util::AppendOutputFile($filename);
returns a file handle

=cut
sub AppendOutputFile {
  ### usage $output_stream = &OpenOutputFile($filename);
  my ($filename) = @_;
  my $to_open = "";
  my $fh;

  if ($filename) {
    if ($filename =~ /\.gz$/) { ### gzip file -> decompress it on the fly
      $to_open = "| gzip -c > $filename";
    } else {
      $to_open = ">>$filename";
    }
#    warn("HELLO\t", $to_open, "\n");
    open($fh, $to_open) ||
	&RSAT::error::FatalError ("Cannot append to output file ".$to_open);
  } else {
    $fh = *STDOUT;
  }
  return ($fh);
}
################################################################

=pod

=item ConvertStrand

Convert different strand formats into "D" (direct) and "R" (reverse)

=cut
sub ConvertStrand {
    my ($strand_to_convert) = @_;
    ### Heli's notation: F (forward) or R (reverse)
    $strand_to_convert =~ s/F/D/;
    $strand_to_convert =~ s/R/R/;

    ### > OR <
    $strand_to_convert =~ s/>/D/;
    $strand_to_convert =~ s/</R/;

    ### Yeast notation: W (Watson) or C (Crick)
    $strand_to_convert =~ s/W/D/;
    $strand_to_convert =~ s/C/R/;
    return $strand_to_convert;
}

################################################################
##  get Background Color from a score value.  this subroutine was
##  created by Rekin's and adapted to the RSA-tools by Sylvain.  It
##  converts a score value into a color (in RGB) according to 5
##  different color gradients: -blue, red, green, grey and fire (from
##  yellow to red)
sub getBgColorFromOneScore {
    my($score,$min,$max,$log,$gradient)=@_;
    my $rval=0;
    my $gval=0;
    my $bval=0;
    my $new_score = "";

    &RSAT::message::Debug("&RSAT::util::getBgColorFromOneScore", $score, $min, $max, $log, $gradient) if ($main::verbose >= 5);
    $score = 0 unless ($score > 0);
    $min = 0 unless ($min > 0);
    $max = 1 unless ($max > 0);

    ## Avoid illegal divisions by 0 when $min == $ax (this can happen if all the edges have the same weight)
    unless ($max > $min) {
      $min = 0;
    }

    if ($log) {
      $new_score = (log ($score-$min+1))/(log($max-$min+1)); # + pseudo weight (1) to avoid log(0)
    } else {
      $new_score = ($score-$min) /($max-$min);
    }
    if ($gradient eq "green"){ # white-green
      $gval=255;
      $bval=(1-$new_score)*255;
      $rval=(sqrt($gval*$bval));
    }
    if ($gradient eq "red"){ # white-red
      $rval=255;
      $gval=(1-$new_score)*255;
      $bval=(1-$new_score)*255;
    }
    if ($gradient eq "blue"){# white-blue
      $bval=255;
      $gval=(1-$new_score)*255;
      $rval=(1-$new_score)*255;
    }
    if ($gradient eq "fire"){ # yellow-red
      $rval=255;
      $gval=(1-$new_score)*255;
#      $bval=(sqrt($gval*$bval));
      $bval=sqrt($gval*$bval);
    }
    if ($gradient eq "grey"){ # grey
      $rval=(1-$new_score)*255;
      $gval=(1-$new_score)*255;
      $bval=(1-$new_score)*255;
    }

    &RSAT::message::Info(join(" ","Score: ",$new_score,"Red:",$rval,"Green:",$gval,"Blue:",$bval)) if ($main::verbose >= 4);
    my $hex = &rgb2hex($rval,$gval,$bval);

    &RSAT::message::Info(join(" ","Hexcolor :",$hex))  if ($main::verbose >= 4);
    return($hex);
}

######################################################
## Convert RGB code to hexadecimal
sub rgb2hex{
  my(@rgb)=@_;
  my $hex = "#";
  my $i=0;
  for($i=0; $i<3; $i++) {
    if( ($rgb[$i] > 255) || ($rgb[$i] < 0) ) {
      return(0);
    }
    $tmp = sprintf("%x", $rgb[$i]);
    if(length($tmp) < 2){
      $hex .= "0". $tmp;
    }
    else{
      $hex .= $tmp;
    }
  }
  return $hex;
}

######################################################
## Convert hexadecimal code (#RRGGBB) to RGB colors (R, G, B) where
## each channel is a number between 0 and 255
sub hex2rgb {
  my($hex_color)=@_;
  my ($r,$g,$b) = ();
  if ($hex_color =~ /^#(\S\S)(\S\S)(\S\S)$/) {
    $r = hex($1);
    $g = hex($2);
    $b = hex($3);
  } else {
    &RSAT::error::FatalError("&RSAT::util::hex2rgb", $hex_color, "is not a valid hexadecimal color specification (format: #RRGGBB)");
  }

  return ($r,$g,$b);
}

################################################################
## Return the root of the public temporary directory
sub get_pub_temp {
  my $public_temp_dir = $ENV{RSAT}."/public_html/tmp";
  return ($public_temp_dir);
}

################################################################
## Return a user-specific directory for storing temporary files.
##
## By default, temporary files are stored in a hidden folder
## $HOME/.rsat_tmp_dir. If $HOME not writable call newdir (docker)
##
## For the Web server, temporary files are stored in
## $RSAT/public_html/tmp, in order to be accessible to web browsers.
sub get_temp_dir {
  my ($sec, $min, $hour,$day,$month,$year) = localtime(time());
  my $login = getpwuid($<) || "temp_user";
  my $tmp_base;

  if ((defined($ENV{RSA_OUTPUT_CONTEXT})) &&
      (($ENV{RSA_OUTPUT_CONTEXT}eq "cgi") || ($ENV{RSA_OUTPUT_CONTEXT} eq "RSATWS"))) {
    $tmp_base = &get_pub_temp()."/".$login;
  } else {
    $tmp_base = $ENV{HOME}."/.rsat_tmp_dir";
    if( ! -w $ENV{HOME}) {
      $tmp_base = tempdir( CLEANUP => 1 );
    }
  }

  my $tmp_dir = sprintf("%s/%04d/%02d/%02d", $tmp_base, 1900+$year,$month+1,$day);

  &RSAT::message::Info("&RSAT::util::get_temp_dir()", $tmp_dir) if ($main::verbose >= 5);

  return($tmp_dir);
}

################################################################
## Return a unique name for a temporary file in the $TMP directory
## Usage:
##  my $temp_file_name = &RSAT::util::make_temp_file($tmp_dir, $tmp_prefix, $add_date, $make_dir);
##    $tmp_dir: if not specified, the default RSAT temporary dir is used
##    $tmp_prefix: prefix for the file name
##    $add_date (value 0 or 1): if 1, the date is added to the suffix
##    $make_dir (value 0 or 1): if 1, create a temporary directory rather than temporary file
sub make_temp_file {
  my ($tmp_dir, $tmp_prefix, $add_date, $make_dir, $protect) = @_;
#    &RSAT::message::Debug("&RSAT::util::make_temp_file()",
#  			"\n\ttmp_dir=".$tmp_dir,
#  			"\n\ttmp_prefix=".$tmp_prefix,
#  			"\n\tadd_date=".$add_date,
#  			"\n\tmake_dir=".$make_dir,
#  		       ) if ($main::verbose >= 10);
  my $prefix_dir = "";

  ## Check $tmp_dir and $tmp_prefix
  if ($tmp_prefix) {
    ($prefix_dir, $tmp_prefix) = &SplitFileName($tmp_prefix);
    if (($tmp_dir) && ($prefix_dir)) {
      $tmp_dir = $tmp_dir."/".$prefix_dir;
    } elsif ($prefix_dir) {
      $tmp_dir = $prefix_dir;
    }
  } else {
    $tmp_prefix = 'tmp';
  }

  $tmp_prefix =~ s/[\(\)\|\:\s]/_/g;
#  &RSAT::message::Debug("&RSAT::util::make_temp_file()", "tmp_prefix", $tmp_prefix) if ($main::verbose >= 10);

  ## Check that temp dir is defined and create it if required
  unless ($tmp_dir) {
    $tmp_dir = &get_temp_dir();
  }

  # die(join("\n\t",
  # 	   "ENV{RSAT}\t".$ENV{RSAT},
  # 	   "&get_pub_temp\t".&get_pub_temp(),
  # 	   "&get_temp_dir\t".&get_temp_dir(),
  # 	   "TMP\t".$TMP,
  # 	   "tmp_dir\t".$tmp_dir,
  # 	   "tmp_prefix\t".$tmp_prefix,
  # 	   "prefix_dir\t".$prefix_dir,
  # 	    "OK"));
  # die "HEREAMI";

  &CheckOutDir($tmp_dir, "", 755); ## temporary dir and all of its parents must be writable by all users

  ## Create an index file in the new directory to prevent Web users
  ## from seing its whole content
  if ($protect) {
    if ((defined($ENV{RSA_OUTPUT_CONTEXT})) &&
	($ENV{RSA_OUTPUT_CONTEXT}eq "cgi")) {
#  if ($ENV{RSA_OUTPUT_CONTEXT} eq "cgi") {
      my $index_file = $tmp_dir."/index.html";
      unless (-e $index_file) {
	open INDEX, ">".$index_file;
	print INDEX "<html>";
	print INDEX "<b>Access forbidden</b>";
	print INDEX "</html>";
	close INDEX;
      }
    }
  }

  ## Add date if required
  if ($add_date) {
    $tmp_prefix .= "_";
    $tmp_prefix .= &AlphaDate();
  }

  ## request the temporary file to the system
  my $mktmp_cmd = "mktemp";
  $mktmp_cmd .= " -u"; ## Ensure that temp file is not created (we just need to return the file name)
  $mktmp_cmd .= " -d " if ($make_dir);
  $mktmp_cmd .= " ".$tmp_dir."/".$tmp_prefix."_XXXXXX";
  my $temp_file = `$mktmp_cmd`;
  chomp($temp_file);

  ## Ensure that everyone can read the temporary file
#  system("chmod a+r $temp_file");
#   &RSAT::message::Debug("&RSAT::util::make_temp_file()",
#  			"\n\ttmp_dir=".$tmp_dir,
#  			"\n\ttmp_prefix=".$tmp_prefix,
#  			"\n\tprefix_dir=".$prefix_dir,
  			"\n\ttemp_file\t".$temp_file,
#  			"\n\tmktmp_cmd=".$mktmp_cmd,
#  		       ) if ($main::verbose >= 10);
  return ($temp_file);
}

################################################################
## echo a command and send it to the system
## Usage:
##   &doit($command, $dry, $die_on_error, $verbose, $batch, $job_prefix, $log_handle, $err_handle, $cluster_queue);
##
sub doit {
  my ($command, $dry, $die_on_error, $verbose, $batch, $job_prefix, $log_handle, $err_handle, $cluster_queue) = @_;

  ## Print the command to the log file
  if ($log_handle) {
    print $log_handle "\n\n", $command, "\n";
  }

  &RSAT::message::Debug("&RSAT::util::doit()", "Command:", $command) if ($main::verbose >= 5);
  unless ($command) {
    &RSAT::message::Warning("&RSAT:::util::doit() was called with empty command. I don't do anything.");
    return();
  }

  ## Define current working directory
  my $wd = $ENV{PWD};

  ## Fix scope problem with the variable $main::verbose
  unless(defined($verbose)) {
    $verbose = $main::verbose;
  }

  if ($batch) {

    ## command used to send the actual command as a script to the job scheduler
    my $qsub_command;

    ## Define the shell
    my $shell = $ENV{CLUSTER_SHELL} || $ENV{SHELL};

    unless ($ENV{CLUSTER_SHELL}) {
      $shell = &RSAT::server::GetProgramPath($shell);
    }


    ## Store the command in a sh script (the job)
    my $job_dir = $wd."/jobs";
    $job_dir .= "/".`date +%Y%m%d`;
    chomp($job_dir);
    &CheckOutDir($job_dir, "", 777);

    my $job_prefix = $job_prefix || "doit";
    &RSAT::util::CheckOutDir($job_dir);
    $job_file = `mktemp ${job_dir}/${job_prefix}.XXXXXX`;
    chomp $job_file;
    my $job_script = "";
    $job_script .=  "#!".$shell."\n";
    $job_script .= "(cd ".$wd;
    $job_script .= "; date > ".$job_file.".started"; ## Write a file called [job].started indicating the time when the job was started
    $job_script .= "; hostname >> ".$job_file.".started";
    $job_script .= "; source ".$ENV{RSAT}."/RSAT_config.bashrc"; ## Required for PERL5LIB
    $job_script .= "; ".$command;
    $job_script .= "; date > ".$job_file.".done"; ## Write a file called [job].done indicating the time when the job was done
    $job_script .= "; hostname >> ".$job_file.".done"; ## Write a file called [job].done indicating the time when the job was done
    $job_script .= " )"."\n";
    open JOB, ">$job_file" || die "Cannot write job file ".$job_file;
    print JOB $job_script;
    close JOB;
    &RSAT::message::TimeWarn("wd", $wd, "Job written in file", $job_file) if ($verbose >= 3);

    my $job_name = $job_file;
    $job_name =~ s/\//_/g;
    my $job_log = $job_file.".log";

    &RSAT::message::Debug("\n\twd", $wd, 
			  "\n\tshell", $shell,
			  "\n\tjob_dir", $job_dir,
			  "\n\tjob_file", $job_file,
			  "\n\tjob_log", $job_log
	) if ($main::verbose >= 10);

    ################################################################
    ## Check that the cluster parameters are well defined before
    ## trying to send jobs to a PC cluster

    ## Identify the job scheduler
    my $queue_manager="";
    if (defined($ENV{QUEUE_MANAGER})) {
      $queue_manager=$ENV{QUEUE_MANAGER};

    } elsif (defined($ENV{QSUB_MANAGER})) {
      $queue_manager=$ENV{QSUB_MANAGER};
      &RSAT::message::Warning("Parameter QSUB_MANAGER is obsolete.\nPlease rename if to QUEUE_MANAGER in ".$RSAT."/RSAT_config.props") if ($main::verbose >= 1);

    } else {
      ## If job scheduler has not beed specified
      $ENV{QUEUE_MANAGER} = "batch";
      $queue_manager=$ENV{QUEUE_MANAGER};
      &RSAT::message::Warning("Cluster job scheduler not defined, using the  default value 'torque'.") if ($verbose >= 2);
    }

    ## Check if qsub options have been specified in RSAT properties or
    ## as environment variable.
    my $qsub_options="";
    if (defined($ENV{QSUB_OPTIONS})) {
      $qsub_options = $ENV{QSUB_OPTIONS};
    }


    ## Cluster queue
    unless ($cluster_queue) {
      if ($ENV{CLUSTER_QUEUE}) {
	$cluster_queue = $ENV{CLUSTER_QUEUE};
      } else {
	&RSAT::error::FatalError("In order to send jobs to a PC cluster, you need to define an environment variable CLUSTER_QUEUE.");
      }
    }

    ## Treatment of the user feed-back
    my $batch_mail = $ENV{BATCH_MAIL} || "a";

    ## optional: restrict the jobs to selected nodes
    ## Example: -l nodes=1:k2.6
    my $selected_nodes =$ENV{NODES};
#    my $selected_nodes =$ENV{NODES} || " -l nodes=1:k2.6 "; ##

    ################################################################
    ## Choose the job scheduler depending on the local configuration
    if (lc($queue_manager) eq "torque") {
#      $wd = $ENV{PWD}; #`pwd`;
#      chomp($wd);

      ## qsub command functionning using torque
      $qsub_command = "qsub";
      $qsub_command .= " ".$selected_nodes if ($selected_nodes);
      $qsub_command .= " -V";                     ## Pass environment variables to qsub !
      $qsub_command .= " -d ".$wd;                ## Working directory used for the job
      $qsub_command .= " -m ".$batch_mail;        ## Email options (error, completion)
      $qsub_command .= " -q ".$cluster_queue;     ## Queue
      $qsub_command .= " -N ".$job_file;          ## Job name
      $qsub_command .= " -e ".$job_file.".err";   ## STDER log file
      $qsub_command .= " -o ".$job_file.".log";   ## STDOUR log file
      $qsub_command .= " ".$job_file;             ## Script
      &RSAT::message::Debug("qsub command for torque", $qsub_command) if ($main::verbose >= 2);

    } elsif (lc($queue_manager) eq "sge") {
      ## qsub command functionning using Sun Grid Engine
      $qsub_command = join(" ", "qsub",
		       "-m", $batch_mail,
		       "-q ", $cluster_queue,
		       " -j y ",
		       "-N ", $job_name,
		       "-o ".$job_log,
		       $qsub_options,
		       $job_file);

    } elsif (lc($queue_manager) eq "slurm") {
	## Send command to slurm via sbatch command
	$qsub_command = "sbatch ";
	$qsub_command .= $job_file;
	
    } elsif (lc($queue_manager) eq "batch") {
      ## qsub command functionning using Sun Grid Engine
      $qsub_command = &RSAT::server::GetProgramPath("batch");
      $qsub_command = $qsub_command;
      $qsub_command .= " -f ".$job_file;

    } else {
      &RSAT::error::FatalError($queue_manager,
			       "Invalid job scheduler. Supported: torque | SGE.",
			       "Please define the job scheduler by setting the variable QUEUE_MANAGER in RSAT_config.props.");
    }

    &doit($qsub_command, $dry, $die_on_error,$verbose,0);


  } else {
    ## Verbose: report command
    &RSAT::message::TimeWarn("Working dir", $wd) if ($verbose >= 4);
    &RSAT::message::TimeWarn(&hide_RSAT_path($command))  if ($verbose >= 3);
    ## Send the command to the queue
    unless ($dry) {
#	&RSAT::message::Debug("Running command", $command) if ($main::verbose >= 3);

      eval {
	my $error = system($command);
	1;
      } or do {
	## Prepare error message
	if ($error == -1) {
	  $error_message = "Could not execute the command\n\t".$command;
	} elsif ($error) {
	  $error_message = join("\t", $0, "Error", $error , "occurred during execution of the command:", "\n", $command);
	}

	## Print error message in the error log file if specified
	if ($error) {
	  if ($err_handle) {
	    print $err_handle "\n\n", $error_message, "\n";
	  }
	}
	## Die if requested
	if ($die_on_error) {
	  &RSAT::error::FatalError($error_message);
	} else {
	  &RSAT::message::Warning($error_message);
	  return ($error);
	}
      }
    }
  }
}


################################################################
## Run a Web server query and send two notification email (submission
## and completion, respectively), plus display the URL where results
## will become available after completion.
##
## index=>$index_file  index generated by the called program (e.g. compare-matrices)
sub EmailTheResult {
    my ($command, $recipient, $tmp_file_path, %args) = @_;
    my $delay = "72 hours";
    my $subject = $args{title};
    my $no_email = $args{no_email} || 0;

    my $script_name = &RSAT::util::hide_RSAT_path($0);
#    my $script_name = $0;
#    ## Avoid providing the full path of the RSAT package
#    $script_name =~ s|$ENV{RSAT}\/+public_html\/+||g;
#    $script_name = &RSAT::util::hide_RSAT_path($script_name);
    $script_name =~ s|.*public_html\/+||g; ## Suppress the public_html folder
    $script_name =~ s|\.cgi$||; ## Suppress cgi extension


    ## Check the subject
    if ($subject) {
	$subject =~ s/'/\'/g;
    } else {
      $subject = join " ", "[RSAT]", $script_name, &AlphaDate();
    }

    ## Avoid exposing the full RSAT path in the mail subject
    $subject = &RSAT::util::hide_RSAT_path($subject);

    ## Check the email address
    &RSAT::email::CheckEmailAddress($recipient)  unless ($no_email);;

    ## Temporary file for storing the result
    unless ($tmp_file_path) {
      $tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1);
    }
    my $result_URL = $ENV{rsat_www}."/tmp/";
    if ($args{index}) {
      $result_URL .= &RSAT::util::RelativePath(&RSAT::util::get_pub_temp(), $args{index});
    } else {
      $result_URL .= &RSAT::util::RelativePath(&RSAT::util::get_pub_temp(), $tmp_file_path);
    }
    $result_URL =~ s|//tmp|/tmp|g;

    ## Run the task in nice mode
    $command = "nice -n 19 ".$command;
    my $clean_command = &RSAT::util::hide_RSAT_path($command);

    ## Add escape to '$RSAT' before calling the mail command, to avoid replacement by its actual value
#    $clean_command =~ s|\$|\\\$|;


    ## Indicate the URL of the future result file
    my $submission_message = $subject."\n";
    $submission_message .= "\n\nCommand: ".$clean_command;
    $submission_message .= "\n\nYour task has been submitted to the RSAT server.";
    $submission_message .= "\n\nAfter completion of the task, an email will be sent to ".$recipient unless ($no_email);
    $submission_message .= "\n\nThe result will become available at\n";
    $submission_message .= "\t".${result_URL}."\n\n";

    ## Print the task submission message the Web page
    $html_submission_message = $submission_message;
    $html_submission_message =~ s|(http://\S+)|<a target='_blank' href='$1'>$1</a>|gm;
    $html_submission_message =~ s|\\\$|\$|gm;
    &RSAT::message::Info($html_submission_message);

    ## Send email to notify submission of the task
    unless($no_email)
    {
        if($ENV{starttls}){
            &RSAT::email::send_mail_STARTTLS($submission_message, $recipient, $subject.' ; Job submitted');
        } else {
            &RSAT::email::send_mail($submission_message, $recipient, $subject.' ; Job submitted');
        }
    }

    ## Concatenate the command with the email notification
    my $completion_message;
    $completion_message .= $subject."\n";
    $completion_message .= "\n\nCommand: ".$clean_command;
    $completion_message .= "\n\nYour result is available at the following URL:\n\t${result_URL}";
    $completion_message .= "\nThe result file will remain there for $delay.";

    ## Run the command and send a notification mail if the user provided an email address
    my $send_mail = &RSAT::server::GetProgramPath("send-mail");
    $command .=  " | perl -pe 's|$ENV{RSAT}/(public_html/)*||g' >> ".$tmp_file_path;
    $command .= " 2> ".$args{error_file} if ($args{error_file}); ## JvH: added 2015-03-01; to validate
    $command .= "; echo '${completion_message}' | ".$send_mail." -subject '${subject} ; Job completed' -to ".$recipient unless ($no_email);
    $command .= " &"; ## Run task in background
    # print "\n\n<pre>", &RSAT::util::hide_RSAT_path($command), "</pre>" if ($ENV{rsat_echo} >= 5);
    # &RSAT::message::Debug("&RSAT::util::EmailTheResult() command", $command) if ($main::verbose >= 10);

    system ($command);

    ## Prepare removal of the temporary file
    &RSAT::server::DelayedRemoval($tmp_file_path, $delay);
}


################################################################
## Treat one command, by either executing it, or concatenating it for
## further batch processing
##
## Usage:
##   &one_command($cmd, $print_out, $time_file, %args);
##
## If the variable $print_out is set to 1, the command is printed to
## the output file $main::out.
##
## If the variable $time_file is specified, the execution time is
## measured and stored in this file.
##
## Some additional parameters can be specified with the hash table %args:
##    task=>$depending_task
##        execute the command only if ($main::task{depending_task} != 0)
##
##    out=>$output_handle
##        specify the handle for the log file (where commands are
##        printed before execution when the option print_out is
##        active)
##
sub one_command {
#  my ($cmd, $print_out, $time_file, $err_file) = @_;
  my ($cmd, $print_out, $time_file, %args) = @_;

  &RSAT::message::Debug("RSAT::util::one_command()", $cmd) if ($main::verbose >= 3);

  
  ################################################################
  ## Define a dry option for one_command, because there are two
  ## reasons why we might avoid to actually run the command ("dry" mode):
  ##
  ## 1. Because this option was called in the main work space
  my $local_dry = $main::dry;
  ## 2. Because this command is depending on the activation of a
  ##    specific task by the user (option -task).
  if ($args{task}) {
    $required_task = $args{task};
    $local_dry = 1 unless ($main::task{$required_task});
  }

  ## Store execution time in a file
  if ($time_file) {
      ## The time command has operating-system specific parameters. To
      ## fix this, I first get the operating system and adapt the
      ## command accordingly. This is really tricky, I recognize. When
      ## I have time, I should try to find some Perl function to get
      ## elapsed time.
      my $OS = `uname -a`;
      chomp($OS);
      &RSAT::message::Debug("Adapting time command to OS-specific behaviour", $OS) if ($main::verbose >= 5);
       
      if ( $OS =~ /ubuntu/i ||
	   $OS =~ /debian/i || 
           $OS =~ /linux/i ||
#	   $OS =~ /darwin/i || 
           $OS =~ /bongcam/i ) { # personal path?
	  
	  $cmd = 'time -o '.$time_file.' '.$cmd;

      } elsif ($OS =~ /biow/) {
	  ## Do nothing
      } else {
	  $cmd = 'time ('.$cmd.') >& '.$time_file;
      }
  }

#   ## Store STDERR in a file
#   if ($err_file) {
#       $cmd = $cmd." >&".$err_file;
#   }


  ################################################################
  ## If the main method is in batch, append the command to the string
  ## $main::batch_cmd.  This serves to collect several steps of an
  ## analysis in a single command that will be sent to the cluster as
  ## a single job.
  if ($main::batch) {
    if ((defined($main::batch_cmd)) &&
	 ($main::batch_cmd =~/\S/)) {
      $main::batch_cmd .= " ; $cmd";
    } else {
      $main::batch_cmd = "$cmd";
    }
  } else {

    ## Report command in the log file
    if ($print_out) {
      my $local_out;
      if ($args{log}) {
	$local_out = $args{log};
#      } elsif ($main::out) {
#	$local_out = $main::out;
#      } else {
#	$local_out = STDOUT;
#      }
	print $local_out ("\n", "; ", &AlphaDate(), "\n", &hide_RSAT_path($cmd), "\n\n");
      }
    }
    &doit($cmd, $local_dry, $main::die_on_error, $main::verbose, $main::batch, $main::job_prefix, $log_handle);
  }
}


=pod

=item PrintArguments()

print the command-line arguments

=cut
sub PrintArguments {
  my ($local_out) = @_;
  my $argument_string = "";
  foreach my $a (@main::ARGV) {
    if (($a =~ /\s+/)  ||
	($a !~ /\S+/) ||
	($a =~ /[\(\)\>\<\&]/)) {
      $argument_string .= " '$a'";
    } else {
      $argument_string .= " $a";
    }
  }
  $argument_string = &hide_RSAT_path($argument_string);

  print $local_out $argument_string, "\n" if ($local_out);
  return $argument_string;
}


################################################################
## Print threshold values for verbosity
sub PrintThresholdValues {
    my %th = (%lth, %uth);
    my @keys = keys %th;
    my $message = "";
    if (scalar(@keys > 0)) {
	$message .= sprintf "; Threshold values\n";
	$message .= sprintf  ";\t%-15s\t%s\t%s\n", "Parameter", "Lower", "Upper";
	foreach my $key (@keys) {
	    my $lth = "none";
	    my $uth = "none";
	    if (defined($lth{$key})) {
		$lth = $lth{$key};
	    }
	    if (defined($uth{$key})) {
		$uth = $uth{$key};
	    }
	    $message .= sprintf  ";\t%-15s\t%s\t%s\n", $key, $lth, $uth;
	}
    }
    return $message;
}

sub sort_unique {
    my (@input_list) = @_;
    &RSAT::message::Debug("&RSAT::util::sort_unique()", "\@input_list", scalar(@input_list)) if ($main::verbose >= 5);
    my $tmp_list_index = ();
    foreach my $element (@input_list) {
	$tmp_list_index{$element}++;
    }
    my @sorted_list = sort(keys(%tmp_list_index));
    undef(%tmp_list_index); ## There is some strange side effect: I have to explicitly undef the index
    &RSAT::message::Debug("&RSAT::util::sort_unique()", "\@sorted_list", scalar(@sorted_list)) if ($main::verbose >= 5);
    return(@sorted_list);
}



1;

__END__
