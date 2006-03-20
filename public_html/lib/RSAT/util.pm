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
    my ($sec, $min, $hour,$day,$month,$year) = localtime(time);

    ###$my_date = sprintf("%02d%02d%02d.%02d%02d%02d", $year,$month+1,$day,$hour, $min, $sec);
    $my_date = `date +%Y_%m_%d.%H%M%S`;
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
    
    if ($main::verbose >= 2) {
	my $wd = `pwd`;
	warn "; Current directory\t", $wd, "\n";
    }

    if ($output_dir) {
	if (-d $output_dir) {
	    warn "; Directory $output_dir already exists\n" if ($main::verbose >= 2);
	    return;
	}
	warn "; Creating directory $output_dir\n" if ($main::verbose >= 1);
	mkdir ($output_dir, 0755);
	unless (-d $output_dir) {
	    warn "Creating directory with all parents $output_dir\n" if ($main::verbose >= 1);
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
=pod

=item PrintArguments()

print the command-line arguments

=cut
sub PrintArguments {
    my $local_out = $_[0];
#    unless ($local_out)  {
#	$local_out = STDOUT;
#    }
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


