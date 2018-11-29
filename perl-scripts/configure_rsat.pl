#!/usr/bin/env perl

################################################################
## This script permits to interactively define the environment
## variables and parameters that will be used by the RSAT progams.
## 
## These parameters are stored in different files for different
## purposes:
##
## RSAT_config.props
##    config file read by RSAT programs in various languages: Perl,
##    python, java
##
## RSAT_config.mk
##    environment variables for the makefiles
##
## RSAT_config.bashrc
##    environment variables that should be loaded in the (bash) shell
##    of RSAT users. There is currently no support for csh or tcsh,
##    but the file can easily be convered to obtain a non-bash cshrc
##    file.
##
## RSAT_config.conf
##    RSAT configuration for the Apache web server.

our %prev_param = ();
our %new_param = ();


################################################################
## List of file extensions for config files. 
our @props_extensions = ("props", "mk", "bashrc", "conf");

## Indicate, for each extension of config file, whether the user
## should be prompted for variable values.
our %auto_extension = ();
$auto_extension{props} =0;
$auto_extension{mk} =0;
$auto_extension{bashrc} =1;
$auto_extension{conf} =1;

## Variables that should be the same for all extensions
our %cross_ext_variable = ();
$cross_ext_variable{rsat} = 1;
$cross_ext_variable{rsat_server_admin} = 1;
$cross_ext_variable{rsat_admin_email} = 1;
$cross_ext_variable{rsat_site} = 1;
$cross_ext_variable{rsat_bin} = 1;
$cross_ext_variable{group} = 1; ## Taxonomic group specificity
$cross_ext_variable{rsat_tmp} = 1;
$cross_ext_variable{rsat_www} = 1;
$cross_ext_variable{rsat_ws} = 1;
#$cross_ext_variable{rsat_ws_tmp} = 1;
$cross_ext_variable{rsat_bin} = 1;
$cross_ext_variable{qsub_manager} = 1;
$cross_ext_variable{queue_manager} = 1;
$cross_ext_variable{qsub_options} = 1;
$cross_ext_variable{cluster_queue} = 1;
$cross_ext_variable{cluster_sell} = 1;
$cross_ext_variable{ensembl_release} = 1;
$cross_ext_variable{ensemblgenomes_release} = 1;

## Parse command-line arguments
our @arguments = @ARGV;
our %force_param = ();
our $auto_mode = 0;
while ($arg = shift(@arguments)) {
  if (lc($arg) =~ /^\-*auto$/) {
    $auto_mode = 1;
    foreach my $extension (@props_extensions) {
      $auto_extension{$extension} = 1;
    }
  } elsif ($arg =~ /^(\S+)\=(\S+)$/) {
    my $key = lc($1);
    my $value = $2;
    $force_param{$key} = $value;
    print("Setting property\t$key\t$value\n");
  }
}


# ## First argument
# if (exists($ARGV[0])) {
#   if ($ARGV[0] eq "auto") {
#     foreach my $extension (@props_extensions) {
#       $auto_extension{$extension} = 1;
#     }
#     ## Print the help message
#   } elsif (scalar(@ARGV) > 0) {
#     warn join ("\t", "\n", "!!!!  Invalid argument  !!!!", $ARGV[0]), "\n";
#     &PrintHelp();
#   }
# } else {
#   warn "Entering manual mode...\n";
# }

package main;
{

  ## BEWARE: this script MUST be executed from the rsat directory,
  ## because the RSAT path is guessed from the current directory.
  unless ($rsat_path) {
    my $pwd = `pwd`;
    chomp($pwd);
    if (($pwd =~ /rsat\/*$/) || $pwd =~(/rsa-tools/)) {
      $rsat_path = $pwd;
    } elsif ($ENV{RSAT}) {
      $rsat_path = $ENV{RSAT};
    }
  }

  ## Prompt for the RSAT path
  unless ($auto_mode) {
    print "\nAbsolute path to the RSAT package ? [", $rsat_path, "] ";
    my $answer = <STDIN>;
    chomp($answer);
    if ($answer) {
      $rsat_path = $answer;
    }
  }

  ## Compute RSAT parent path
  $rsat_parent_path = `dirname $rsat_path`;
  chomp($rsat_parent_path);
  warn "RSAT parent path\t", $rsat_parent_path, "\n";

  ## Check that the RSAT path seems correct
  unless ($rsat_path =~ /rsat\/*/) {
    warn ("\nWarning: $rsat_path does not seem to be a conventional RSAT path (should terminate by rsat).", "\n\n");
  }

  ## Check that the RSAT path exists and is a directory
  unless (-d $rsat_path) {
    die ("\nError: invalid RSAT path\n\t", $rsat_path, "\nDoes not correspond to an existing directory on this computer", "\n\n");
  }

  ## Treat successively the two configuration files: .props (for Perl
  ## and php scripts) and .mk (for make scripts).
  warn("\n", "We will now edit configuration files in interactive mode, for the ", 
       scalar(@props_extensions), " following extensions: ", join(", ", @props_extensions), "\n");
  
  for my $extension (@props_extensions) {

    ## Check that the config file exists in the RSAT path
    my $template_file = $rsat_path."/RSAT_config_default.".$extension;
    my $config_file = $rsat_path."/RSAT_config.".$extension;

    warn("\n\n\n", "################################################################\n", 
	 "## Editing \".${extension}\" configuration file\t", $config_file,"\n",
	 "## Based on template file\t", $config_file,"\n\n");

    if (-f $config_file) {
      ## Load properties from existing configuration file
      if ($extension eq "props")  {
	open CONFIG, $config_file || die "\n\nCannot read config file\t", $config_file, "\n\n";
	while (<CONFIG>) {
	  if ((/^(\S+)=(.*)/) && !(/^#/)) {
	    my $key = $1;
	    my $value = $2;
	    
	    ## Replace the RSAT parent path if required (at first installation)
	    $value =~ s/\[RSAT_PARENT_PATH\]/${rsat_parent_path}/;
	    
	    ## Store the parameter value in a hash table
	    $prev_param{lc($key)} = $value;
##		    warn join("\t", $key, $value, $prev_param{lc($key)}), "\n"; ## DEBUG
	    
	  }
	}
	close CONFIG;
      }
    } else {
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
    unless ($auto_mode) {
      warn "\nPLEASE CHECK THE FOLLOWING LINE BEFORE GOING FURTHER\n";
      print "\nReady to update config file\t", $config_file, " [y/n] (y): ";
      chomp($answer = <STDIN>);
      $answer = "y" unless ($answer);
      unless ($answer eq "y") {
	warn("\nWARNING: Since you did not answer 'y', the edition of config file ${config_file} is aborted.\n");
	die ("Good bye\n\n");
      }
    }



    ## Create a copy of the config file
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    my $backup_dir = $rsat_path."/backups";
    mkdir($backup_dir) unless (-d $backup_dir);
    my $config_file_bk = $backup_dir."/RSAT_config_bk_".($year+1900)."-".($mon+1)."-".$mday."_".$hour."-".$min."-".$sec.".".$extension;
    warn ("\n\nBackup of previous config file\t", $config_file_bk, "\n\n");
    system("cp ".$config_file." ".$config_file_bk);


    ## Open a new file for writing the new config
    my $new_config_file = $config_file.".updated";
    open NEW_CONF, ">".$new_config_file || die "\n\nCannot write new config file\t", $new_config_file, "\n\n";

    ## Load the RSAT config file
    $from_template = 1;
    if ($from_template) {
      open CONFIG, $template_file || die "\n\nCannot read template config file\t", $template_file, "\n\n";
    } else {
      open CONFIG, $config_file || die "\n\nCannot read config file\t", $config_file, "\n\n";
    }
    while (<CONFIG>) {

      ## Treat the Apache config file
      if ($extension eq "conf") {
	## For apache, the only change is to replace[RSAT_PARENT_PATH]
	## by the actual path
	s/\[RSAT_PARENT_PATH\]/${rsat_parent_path}/;
	print NEW_CONF;

      } else {

	if ((/(\S+)=(.*)/) && !(/^#/)) {
	  my $key = $1;
	  my $value = $2;

	  ## Replace the RSAT parent path if required (at first installation)
	  $value =~ s/\[RSAT_PARENT_PATH\]/${rsat_parent_path}/;

	  ## Replace the RSAT web server path if required (at first installation)
	  if ($key eq "rsat_www") {
	    $value .= "/" if ($value =~/^http/);
	    $value =~ s|\/\/$|/|;
	  } elsif (($prev_param{rsat_www}) && ($new_param{rsat_www})
		   && ($value =~ /$prev_param{rsat_www}/)
		   && ($new_param{rsat_www} ne $prev_param{rsat_www})) {
	    $value =~ s|$prev_param{rsat_www}|$new_param{rsat_www}|;
	  }

	  if ($extension eq "props") {
	    ## Select previous parameter values if defined in the props
	    ## file
##	  warn join("\t", $key, $value, "previously", $prev_param{lc($key)}), "\n"; DEBUG
	    if (defined($force_param{lc($key)})) {
	      $value = $force_param{lc($key)};
	    } elsif (defined($prev_param{lc($key)})) {
	      $value = $prev_param{lc($key)};
	    } else {
	      $prev_param{lc($key)} = $value;
	    }
	    
	  } elsif (($extension eq "mk") || ($cross_ext_variable{lc($key)})) {
	    ## Transmit appropriate variables across extensions, if
	    ## already specified.
	    
	    if (defined($force_param{lc($key)})) {
	      $value = $force_param{lc($key)};
	    } elsif (defined($new_param{lc($key)})) {
	      $value = $new_param{lc($key)};
	    }
	  }
	  
	  ## Prompt for the new value
	  unless ($auto_extension{$extension}) {
	    print "\n", $key, " [", $value, "] : ";
	    chomp(my $new_value = <STDIN>);
	    if ($new_value) {
	      $value = $new_value;
	    }
	  }

	  ## Export the line in the new config file
	  if ($extension eq "bashrc") {
	    print NEW_CONF "export ", $key, "=", $value, "\n";
	  } else {
	    print NEW_CONF $key, "=", $value, "\n";
	  }
	  $new_param{lc($key)} = $value;

#	warn join ("\t", "key=".$key, "value=".$value, "param=".$new_param{lc($key)}, "previous=".$prev_param{lc($key)}), "\n";

	} else {
	  print;			## Display comments
	  print NEW_CONF;
	}
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

Usage: 
  cd \$RSAT; perl perl-scripts/configure_rsat.pl [auto]

The option "auto" suppresses the interactive control of parameter
values.

End_of_help
  exit(0);
}
