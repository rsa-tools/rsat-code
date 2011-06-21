#!/usr/bin/perl -w
############################################################
#
# $Id: inferoperon_pathway_extractor.pl,v 1.2 2011/06/21 09:11:21 rsat Exp $
#
############################################################

## use strict;

=pod

=head1 Infer Operon Pathway Extractor

=head1 VERSION 1.0

=head1 DESCRIPTION

Connection genes in a metabolic graph from infer Operon output file. 

=head1 AUTHORS

didier.croes@ulb.ac.be

=head1 USAGE

inferoperon_pathway_extractor.pl -h [-i inputfile] [-o outputdirectory] [-v verbosity] [-n name of the graph] -g graphfile -a gene2ec2reactions [-b gene2reaction]

=head1 INPUT FORMAT

=head2   1) -a gene2ec2reactions [-b gene2reaction]
  gene_id ec_number       reaction_id     species_name    taxonomy_id     gene_name
  O22340  4.2.3.- RXN-10482       Abies grandis   46611   (4S)-limonene synthase
  O22340  4.2.3.- RXN-10483       Abies grandis   46611   (4S)-limonene synthase
  O22340  4.2.3.- RXN-10566       Abies grandis   46611   (4S)-limonene synthase
  O22340  4.2.3.- RXN-10567       Abies grandis   46611   (4S)-limonene synthase
  O22340  4.2.3.- RXN-10568       Abies grandis   46611   (4S)-limonene synthase
  O22340  4.2.3.- RXN-10600       Abies grandis   46611   (4S)-limonene synthase

=cut


BEGIN {
    if ($0 =~ /([^(\/)]+)$/) {
	push (@INC, "$`lib/");
    }
}
require "RSA.lib";



################################################################
## Main package
package main;
{

    ################################################################
    ## Initialise parameters
    local $start_time = &RSAT::util::StartScript();
    $program_version = do { my @r = (q$Revision: 1.2 $ =~ /\d+/g); sprintf "%d."."%02d" x $#r, @r };
#    $program_version = "0.00";

    %main::infile = ();
    %main::outfile = ();

    $main::verbose = 0;
    $main::in = STDIN;
    $main::out = STDOUT;
    
    $main::gprfile ="METACYC_GPR_EC.tab";	# GPR Gene -> EC -> REACTION annotation file path. Default (METACYC_GPR_EC.tab)
    $main::grfile="";				# GR Gene -> REACTION annotation
    $main::graph = "";				# Graph Name 
    $main::graphfile="";				# File containing the graph
    $main::show;				# Open png image in gwenview

    ################################################################
    ## Read argument values
    &ReadArguments();

    ################################################################
    ## Check argument values
    if (!($outdir=~m/\/$/)) {
      $outdir = $outdir."/";
    }
    ################################################################
    ## Open output stream
    $main::out = &OpenOutputFile($main::outfile{output});

    ################################################################
    ## Read input
    ($main::in) = &OpenInputFile($main::infile{input});
   

    ################################################################
    ## Print verbose
    &Verbose() if ($main::verbose);

    ################################################################
    ## Execute the command
     print $outdir;
     my $oldoperonid;
     my @genelist;
     my $line;
     while ($line=<$main::in>) {
	 
	 chomp($line);
	 if (!($line=~m/^;|^#/)){
# 	   print "$line\n";
	   my @tempdata = split(/\t/,$line);
# 	   print "$tempdata[6] eq $oldoperonid\n";
	   if ($oldoperonid && !($tempdata[6] eq $oldoperonid)){
	    $oldoperonid=~s/\s+/_/g;
	    my $tempfilename = `mktemp $outdir$oldoperonid.XXXX`;
	    chomp ($tempfilename);
	    open (MYFILE, '>>'.$tempfilename);
	    print MYFILE join ("\n",@genelist). "\n";
	    close (MYFILE);
	    my $pathway_infer_cmd = "pathway_extractor.pl -i $tempfilename -g $main::graphfile  -b $main::grfile -a $main::gprfile -o $main::outdir";
	    system "$pathway_infer_cmd\n";
	    system "rm $tempfilename";
	    @genelist =();
# 	    exit(0); 
	   }
	   my $gene= $tempdata[0];
	   $gene=~ s{\.[^.]+$}{}; # removes gene version
	   push(@genelist,$gene);
	   $oldoperonid = $tempdata[6];
	 }
	  
    }
    close $main::in if ($main::infile{input});



    ################################################################
    ## Insert here output printing

    ################################################################
    ## Report execution time and close output stream
    my $exec_time = &RSAT::util::ReportExecutionTime($start_time); ## This has to be exectuted by all scripts
    print $main::out $exec_time if ($main::verbose >= 1); ## only report exec time if verbosity is specified
    close $main::out if ($main::outfile{output});

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

=item B<-hp>

Display full PathwayInference help message

=cut
    } elsif ($arg eq "-hp") {
      system "java graphtools.algorithms.Pathwayinference -h";
=pod

=item B<-show>

execute gwenview to display the pathway results in png format
=cut
    } elsif ($arg eq "-show") {
     $main::show ="true";

=pod

=item B<-i inputfile>

If no input file is specified, the standard input is used.  This
allows to use the command within a pipe.

=cut
    } elsif ($arg eq "-i") {
      $main::infile{input} = shift(@arguments);
    
=pod

=item	B<-a GPR Genes file Default (METACYC_GPR_EC.tab)>

GPR Gene -> EC -> REACTION annotation file path. Default (METACYC_GPR_EC.tab)

=cut
    } elsif ($arg eq "-a") {
      $main::gprfile = shift(@arguments);
    
=pod

=item	B<-b GR Gene -> REACTION annotation>

An gene annotation file with diredt link gene to reaction. Does not rely on the EC number annotation

=cut
    } elsif ($arg eq "-b") {
      $main::grfile = shift(@arguments);
=pod

=item	B<-n Graph name>

Name of the Graph (default: Name of the graph file) 

=cut
    } elsif ($arg eq "-n") {
      $main::graph = shift(@arguments);

=pod

=item	B<-d Graph file>

Name of the Graph (default: Name of the graph file) 

=cut
    } elsif ($arg eq "-d") {
      $main::groupdescriptor = shift(@arguments);


=pod

=item	B<-d Unique descriptor>

Unique name to differenciate output files. If not set With -i, the name of the input file will be used.  

=cut
    } elsif ($arg eq "-g") {
      $main::graphfile = shift(@arguments);
=pod

=item	B<-o output Directory>

If no output file is specified, the current directory is used. 

=cut
    } elsif ($arg eq "-o") {
      $main::outdir = shift(@arguments);

    } else {
      &FatalError(join("\t", "Invalid option", $arg));

    }
  }

=pod

=back

=cut

}

################################################################
## Verbose message
sub Verbose {
    print $main::out "; template ";
    &PrintArguments($main::out);
    printf $main::out "; %-22s\t%s\n", "Program version", $program_version;
    if (%main::infile) {
	print $main::out "; Input files\n";
	while (my ($key,$value) = each %main::infile) {
	  printf $main::out ";\t%-13s\t%s\n", $key, $value;
	}
    }
    if (%main::outfile) {
	print $main::out "; Output files\n";
	while (my ($key,$value) = each %main::outfile) {
	  printf $main::out ";\t%-13s\t%s\n", $key, $value;
	}
    }
}


__END__
