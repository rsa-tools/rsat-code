#!/usr/bin/perl -w

################################################################
##

our %param = ();

if (scalar(@ARGV)) {
  &PrintHelp();
}

package main;
{

  ## Check if the RSAT environment variable has been specified
  $param{RSAT} = $ENV{RSAT};

  ## Prompt for the RSAT
  print "Absolute path to the RSAT package ? [", $param{RSAT}, "] ";
  chomp(my $answer = <>);
  if ($answer) {
    $param{RSAT} = $answer;
  }

  ## Check that the RSAT path seems correct
  unless ($param{RSAT} =~ /rsa\-tools\/*/) {
    warn ("\nWarning: $param{RSAT} does not seem to be a conventional RSAT path (should terminate by rsa-tools).", "\n\n");
  }

  ## Check that the RSAT path exists and is a directory
  unless (-d $param{RSAT}) {
    die ("\nError: invalid RSAT path\n\t", $param{RSAT}, "\nDoes not correspond to an existing directory on this computer", "\n\n");
  }

  ## Check that the config file exists in the RSAT path
  my $config_file = $param{RSAT}."/RSAT_config.props";
  unless (-f $config_file) {
    die ("\nError: the config file RSAT_config.props is not found in the RSAT path\n\t", $param{RSAT},
	 "\nPlease check that the RSAT package has been properly installed in this directory.",
	 "\n\n");
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
      print $key, " [", $value, "] ";
      chomp(my $new_value = <>);
      if ($new_value) {
	$value = $new_value;
      }

      print NEW_CONF $key, "=", $value, "\n";
    } else {
      print NEW_CONF;
    }
  }

  close CONFIG;
  close NEW_CONF;

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
