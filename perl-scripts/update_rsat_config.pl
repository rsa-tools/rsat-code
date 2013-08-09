#!/usr/bin/perl -w

################################################################
##

our %param = ();

if (scalar(@ARGV)) {
  &PrintHelp();
}

package main;
{


  my @extensions =  ("props", "mk");

  ## Check if the RSAT environment variable has been specified
  $rsat_path = $ENV{RSAT};

  ## Try to guess RSAT path if not specified in the environment variable
  unless ($rsat_path) {
    my $pwd = `pwd`;
    chomp($pwd);
    if ($pwd =~ /rsa\-tools\/*$/) {
      $rsat_path = $pwd;
    }
  }

  ## Prompt for the RSAT path
  print "\nAbsolute path to the RSAT package ? [", $rsat_path, "] ";
  chomp(my $answer = <>);
  if ($answer) {
    $rsat_path = $answer;
  }

  ## Compute RSAT parent path
  $rsat_parent_path = `dirname $rsat_path`;
  chomp($rsat_parent_path);
  warn "RSAT parent path\t", $rsat_parent_path, "\n";

  ## Check that the RSAT path seems correct
  unless ($rsat_path =~ /rsa\-tools\/*/) {
    warn ("\nWarning: $rsat_path does not seem to be a conventional RSAT path (should terminate by rsa-tools).", "\n\n");
  }

  ## Check that the RSAT path exists and is a directory
  unless (-d $rsat_path) {
    die ("\nError: invalid RSAT path\n\t", $rsat_path, "\nDoes not correspond to an existing directory on this computer", "\n\n");
  }

  ## Treat successively the two configuration files: .props (for Perl
  ## and php scripts) and .mk (for make scripts).

  for my $extension (@extensions) {

    ## Check that the config file exists in the RSAT path
    my $config_file = $rsat_path."/RSAT_config.${extension}";
    warn("\n\nEditing ${extension} configuration file\t", $config_file,"\n\n");

    unless (-f $config_file) {
      my $default_config_file = $rsat_path."/RSAT_config_default.${extension}";
      if (-e $default_config_file) {
	warn ("\nThe config file RSAT_config.${extension} is not found in the RSAT path\n\t", $rsat_path,
	      "\nCopying from default config file\t", $default_config_file,
	      "\n\n");
	system("cp ".$default_config_file." ".$config_file);
      } else {
	die ("\nError: the config file RSAT_config.${extension} is not found in the RSAT path\n\t", $rsat_path,
	     "\nPlease check that the RSAT package has been properly installed in this directory.",
	     "\n\n");
      }
    }

    ## Prompt for the new value
    warn "\nPLEASE CHECK THE FOLLOWING LINE BEFORE GOING FURTHER\n";
    print "\nReady to update config file\t", $config_file, " [y/n] (n): ";
    chomp($answer = <>);
    unless ($answer eq "y") {
      warn("\nWARNING: Since you did not answer 'y', the edition of config file ${config_file} is aborted.\n");
      die ("Good bye\n\n");
    }

    open CONFIG, $config_file || die "\n\nCannot read config file\t", $config_file, "\n\n";

  ## Create a copy of the config file
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    my $config_file_bk = $config_file.".bk.".($year+1900)."-".($mon+1)."-".$mday."_".$hour."-".$min."-".$sec;
    warn ("\n\nBackup of previous config file\t", $config_file_bk, "\n\n");
    system("cp ".$config_file." ".$config_file_bk);

    ## Open a new file for writing the new config
    my $new_config_file = $config_file.".updated";
    open NEW_CONF, ">".$new_config_file || die "\n\nCannot write new config file\t", $new_config_file, "\n\n";

    ## Load the RSAT config file
    while (<CONFIG>) {
      if ((/(\S+)=(.*)/) && !(/^#/)) {
	my $key = $1;
	my $value = $2;
	$param{$key} = $value;

	## Prompt for the new value
	print "\n", $key, " [", $value, "] : ";
	chomp(my $new_value = <>);
	if ($new_value) {
	  $value = $new_value;
	}

	print NEW_CONF $key, "=", $value, "\n";
      } else {
	print;			## Display comments
	print NEW_CONF;
      }
    }

    close CONFIG;
    close NEW_CONF;

    system ("mv -f ".$new_config_file." ".$config_file);
    warn ("\n\nBackup of previous config file\n\t", $config_file_bk, "\n");
    warn ("Updated config file\n\t", $config_file."\n\n");
  }

  exit(0);
}


sub PrintHelp {
  print <<End_of_help;

This script allows to update the RSAT config file in an interactive
way. It should be used by RSAT amdinistrators when the configuration
has to be changed (example: change of the IP address of the server).

Author: Jacques.van-Helden\@univ-amu.fr

usage: perl update_rsat_config.pl

End_of_help
  exit(0);
}
