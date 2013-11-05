#!/usr/bin/perl -w
############################################################
#
# $Id: supported-organisms-ucsc.pl,v 1.2 2012/03/14 08:39:08 jeremy Exp $
#
############################################################

## use strict;

=pod

=head1 NAME

suuported-organism-uscs

=head1 VERSION

$program_version

=head1 DESCRIPTION

Retrieve organisn disponible on UCSC Genome browser
(http://genome.ucsc.edu/).

=head1 USAGE

supported-organims-uscs.pl [-o file] [-v #] [...]

Examples

Retrieve mammal genome disponible on UCSC.

 supported-organims-uscs.pl -taxon mammal

=cut

BEGIN {
  if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
  }
}
require "RSA.lib";

use Bio::Das; ## Required to access UCSC Genome Browser

## Main package
package main;
{

  ################################################################
  ## Initialise parameters
  our $start_time = &RSAT::util::StartScript();
  our $program_version = do { my @r = (q$Revision: 1.2 $ =~ /\d+/g); sprintf "%d."."%02d" x $#r, @r };
  #    $program_version = "0.00";

	our $taxon = "";
	our $id_ok = "";
  our %outfile = ();
  our $verbose = 0;

  ## Parameters for connecting the DAS server
  our $das_server="http://genome.cse.ucsc.edu/cgi-bin/das";


  ## Hard-coded information about taxonomic grouping.
  ##
  ## This is a temporary fix because we dd not find the way to obtain
  ## the taxonomy from UCSC.

  ##First genome for each taxon
  our %first_id  = ();
  our %last_id  = ();

  $first_id{mammal} = "hg19";
  $first_id{vertebrate} = "galGal3";
  $first_id{deuterostome} = "braFlo1";
  $first_id{insect} = "dm3";
  $first_id{nematode} = "ce6";
  $first_id{other} = "aplCal1";

  ##Last genome for each taxon
  $last_id{mammal} = "ornAna1";
  $last_id{vertebrate} = "petMar1";
  $last_id{deuterostome} = "strPur1";
  $last_id{insect} = "apiMel1";
  $last_id{nematode} = "priPac1";
  $last_id{other} = "sacCer1";

  ################################################################
  ## Read argument values
  &ReadArguments();

  ################################################################
  ## Define the URL fo the DAS server
  our $das_server_url = $das_server;
  &RSAT::message::Info("DAS server", $das_server_url) if ($main::verbose >= 2);

  ## Open DAS client
  my $das = Bio::Das->new(5);  # timeout of 5 sec

  ################################################################
  ## Open output stream
  $out = &OpenOutputFile($outfile{output});

  ################################################################
  ## Define First and Last id
#  if ($taxon eq "mammal") {
#      our $first_id = $first_id_mammal;
#      our $last_id = $last_id_mammal;
#  } elsif ($taxon eq "vertebrate") {
#      our $first_id = $first_id_mammal; #mammal are vertebrate
#      our $last_id = $last_id_vertebrate;
#  } elsif ($taxon eq "deuterostome") {
#      our $first_id = $first_id_mammal; #vertebrate are deuterostome
#      our $last_id = $last_id_deuterostome;
#  } elsif ($taxon eq "insect") {
#      our $first_id = $first_id_insect;
#      our $last_id = $last_id_insect;
#  } elsif ($taxon eq "nematode") {
#      our $first_id = $first_id_nematode;
#      our $last_id = $last_id_nematode;
#  } elsif ($taxon eq "other") {
#      our $first_id = $first_id_other;
#      our $last_id = $last_id_other;
#  } else {
#      our $first_id = $first_id_mammal;
#      our $last_id = $last_id_other;
#  }

  our $first_id = $first_id{mammal};
  our $last_id = $last_id{other};
  if (defined($first_id{$taxon})) {
      $first_id = $first_id{$taxon};
  }
  if (defined($last_id{$taxon})) {
      $last_id = $last_id{$taxon};
  }

  ################################################################
  ## Print verbose
  #&Verbose() if ($main::verbose >= 1);

  ## Send request to DAS server to obtain the list of supported organisms
  &RSAT::message::TimeWarn("Sending request to DAS server", $das_server_url) if ($main::verbose >= 3);
  my @response = $das->dsn('http://genome.cse.ucsc.edu/cgi-bin/das');
  
  for my $url (@response) {
      if ($url->is_success) {
	  my @dsns = $url->results;
	  foreach (@dsns) {
	      
	      unless ($id_ok) {
		  if ($_->id eq $first_id) {
		      $id_ok=1;
		      print $out $_->id,"\t",$_->description,"\n";
		  } 
		  
	      } else {
		  print $out $_->id,"\t",$_->description,"\n";
		  
		  ##Stop loop if id is the last one of the taxon
		  if ($_->id eq $last_id) {
		      last;
		  }
	      } 
	  }
      } else {
	  &RSAT::message::Warning("DAS request returned error", $request->dsn,": ",$request->error);
      }
  }

  ################################################################
  ## Report execution time and close output stream
  my $exec_time = &RSAT::util::ReportExecutionTime($start_time); ## This has to be exectuted by all scripts
  print $out $exec_time if ($main::verbose >= 1); ## only report exec time if verbosity is specified
  close $out if ($outfile{output});

  exit(0);
}


################################################################
################### SUBROUTINE DEFINITION ######################
################################################################


################################################################
## Display full help message 
sub PrintHelp {
  system "pod2text -c $0";
  exit()
}

################################################################
## Display short help message
sub PrintOptions {
  &PrintHelp();
}

################################################################
## Read arguments 
sub ReadArguments {
  my $arg;
  my @arguments = @ARGV; ## create a copy to shift, because we need ARGV to report command line in &Verbose()
  while (scalar(@arguments) >= 1) {
    $arg = shift (@arguments);
    ## Verbosity

=pod

=head1 OPTIONS

=over 4

=item B<-v #>

Level of verbosity (detail in the warning messages during execution)

=cut
    if ($arg eq "-v") {
      if (&IsNatural($arguments[0])) {
	$main::verbose = shift(@arguments);
      } else {
	$main::verbose = 1;
      }

=pod

=item B<-h>

Display full help message

=cut
    } elsif ($arg eq "-h") {
      &PrintHelp();

=pod

=item B<-help>

Same as -h

=cut
    } elsif ($arg eq "-help") {
      &PrintOptions();

=pod

=item B<-taxon #>

Return only genome from this taxon. 

Supported Taxon.

=over

=item I<mammal>

=item I<vertebrate>

=item I<deuterostome>

=item I<insect>

=item I<nematode>

=item I<other>

=back
	

=cut
    } elsif ($arg eq "-taxon") {
      $main::taxon = shift(@arguments);

=pod

=item	B<-o outputfile>

The output file is in fasta format.

If no output file is specified, the standard output is used.  This
allows to use the command within a pipe.

=cut
    } elsif ($arg eq "-o") {
      $outfile{output} = shift(@arguments);

    } else {
      &FatalError(join("\t", "Invalid option", $arg));

    }
  }
}

