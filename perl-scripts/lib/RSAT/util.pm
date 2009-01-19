###############################################################
#
# A class for util handling
#
package RSAT::util;

use RSAT::GenericObject;
use RSAT::message;
use RSAT::error;
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
    if ($_[0] =~ /^ *\d+ *$/) {
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

Usage : my $i = &RSAT::util::rount($r);

=cut
sub round {
    my ($my_real) = @_;
    my $my_int;
    if (abs(($my_real - int($my_real))) < 0.5) {
	return int($my_real);
    } elsif ($my_real < 0) {
	return int($my_real) - 1;
    } else {
	return int($my_real) + 1;
    }
}


################################################################
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
    my $my_date = "";
    my ($sec, $min, $hour,$day,$month,$year) = localtime(time());
    $my_date = sprintf("%02d_%02d_%02d.%02d%02d%02d", 1900+$year,$month+1,$day,$hour, $min, $sec);
#    $my_date = `date +%Y_%m_%d.%H%M%S`;
    chomp $my_date;
    return $my_date;
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
      #	$dir =~ s|/$||; #### suppress the trailing /
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

################################################################
=pod

=item CheckOutDir()

Check for the existence of an output directory, and, if it does not
exist, create it.

Usage: &RSAT::util::ChechOutDir($my_dir);

=cut
sub CheckOutDir {
    my ($output_dir, $umask) = @_;

    unless ($output_dir) {
	warn "; CheckOutDir\tNo directory has been specified\n" if ($main::verbose >= 1);
	return;
    }

    $umask = 0002 unless ($umask);
    umask($umask);
    if ($main::verbose >= 3) {
	my $wd = `pwd`;
	warn "; Current directory\t", $wd, "\n";
    }

    if ($output_dir) {
	if (-d $output_dir) {
	    warn "; Directory $output_dir already exists\n" if ($main::verbose >= 3);
	    return;
	}
	warn "; Creating directory $output_dir\n" if ($main::verbose >= 2);
	mkdir ($output_dir, 0755);
	unless (-d $output_dir) {
	    warn "Creating directory with all parents $output_dir\n" if ($main::verbose >= 2);
	    system "mkdir -p $output_dir"; ## create output directory with all parents
	}
	unless (-d $output_dir) {
	    &RSAT::error::FatalError("Cannot create output directory $output_dir");
	}
    } else {
	$output_dir = ".";
    }
}


################################################################
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
	    warn "; Uncompressing .gz file\n" if ($main::verbose >= 2);
	    $filename = "gunzip -c $filename |";
				  } elsif ($filename =~ /\.Z$/) { ### gzip file -> decompress it on the fly
				      warn "; Uncompressing .Z file\n" if ($main::verbose >= 2);
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
	open (my $fh, "$to_open") || 
	    &RSAT::error::FatalError ("Cannot write output file ".$filename);
	return ($fh)
    } else {
	$fh = *STDOUT;
	return $fh;
    }
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
### get Background Color from a score value
###  this subroutine was created by Rekin's and adapted to the RSA-tools
###  by Sylvain
###  It converts a score value into a color (in RGB) according to 5 different
###  color gradients
### -blue, red, green, grey and fire (from yellow to red)


sub getBgColorFromOneScore{
    my($score,$min,$max,$log,$gradient)=@_;
    my $rval=0;
    my $gval=0;
    my $bval=0;
    my $new_score = "";
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
      $bval=(sqrt($gval*$bval));
    }
    if ($gradient eq "grey"){ # grey
      $rval=(1-$new_score)*255;
      $gval=(1-$new_score)*255;
      $bval=(1-$new_score)*255;
    }

    &RSAT::message::Info(join(" ","Score: ",$new_score,"Red:",$rval,"Green:",$gval,"Blue:",$bval)) if ($main::verbose >= 4);
    my $hex = &rgb2hex($rval,$gval,$bval );
    
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


################################################################
## echo a command and send it to the system
## Usage:
##   &doit($command, $dry, $die_on_error, $verbose, $batch, $job_prefix);
sub doit {
  my ($command, $dry, $die_on_error, $verbose, $batch, $job_prefix) = @_;
  my $wd = `pwd`;
  chomp $wd;

  if ($batch) {
    ## Store the command in a sh script (the job)
    my $job_dir = "jobs";
    $job_dir .= "/".`date +%Y%m%d`;
    chomp($job_dir);
    my $job_prefix = $job_prefix || "doit";
    &RSAT::util::CheckOutDir($job_dir);
    my $job = `mktemp ${job_dir}/${job_prefix}.XXXXXX`;
    chomp $job;
    open JOB, ">$job";
    print JOB "( cd ", $wd;
    print JOB "; date > ", $job, ".started"; ## Write a file called [job].started indicating the time when the job was started
    print JOB "; hostname >> ", $job,".started"; 
    print JOB "; ", $command;
    print JOB "; date > ", $job, ".done"; ## Write a file called [job].done indicating the time when the job was done
    print JOB "; hostname >> ", $job, ".done"; ## Write a file called [job].done indicating the time when the job was done
    print JOB " )", "\n";
    close JOB;
    &RSAT::message::TimeWarn(join("\t", "Job queued", $wd."/".$job)) if ($main::verbose >= 2);

    my $job_name = $job;
    $job_name =~ s/\//_/g;
    my $job_log = $wd."/".$job.".log";

    ## Send the command to a batch queue on a PC cluster The default
    ## values are for internal use in the BiGRe laboratory, but alternative values
    ## can be specified by specifying the environment variables
    ## CLUSTER_QUEUE and CLUSTER_MASTER
    my $cluster_queue = $ENV{CLUSTER_QUEUE} || "short";
    my $batch_mail=$ENV{BATCH_MAIL} || "a";
    my $qsub_manager=$ENV{QSUB_MANAGER} || "sge";
    $selected_nodes = " -l nodes=1:k2.6 ";
    my $qsub_cmd;

    ################################################################
    ## Choose the queue manager depending on the local configuration
    if ($qsub_manager eq "torque") {
      ## qsub command functionning using Torque
      my $cluster_master=$ENV{CLUSTER_MASTER} || "arthur.bigre.ulb.ac.be"; ## for torque only
      $qsub_cmd = "qsub ".$selected_nodes." -m ".$batch_mail." -q ".$cluster_master." -N ${job} -j oe -o ${job}.log ${job}";

    } else {
      ## qsub command functionning using Sun Grid Engine (BiGRe)
      $qsub_cmd = join(" ", "qsub", 
		       "-m",$batch_mail,
		       "-q ", $cluster_queue, 
		       " -j y ",
		       "-N ", $job_name,
		       "-o ".$job_log, $job);
    }
    &doit($qsub_cmd, $dry, $die_on_error,$verbose,0);

  } else {
    ## Verbose
    if (($dry) || ($main::verbose >= 2)) {
      warn "\n";
      &RSAT::message::TimeWarn("Working dir", $wd) if ($main::verbose >= 3);
      &RSAT::message::TimeWarn($command);
    }

    ## Send the command to the queue
    unless ($dry) {
      my $error = system $command;
      if ($die_on_error) {
	if ($error == -1) {
	  &RSAT::error::FatalError("Could not execute the command\n\t$command");
	} elsif ($error) {
	  &RSAT::error::FatalError("Error", $error , "occurred during execution of the command", "\n\t".$command);
	}
      } else {
	return ($error);
      }
    }
  }
}

################################################################
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
  print $local_out $argument_string, "\n" if ($local_out);
  return $argument_string;
}


return 1;


__END__


