#!/usr/bin/perl -w
############################################################
#
# $Id: genome-blast.pl,v 1.2 2005/04/21 23:10:13 rsat Exp $
#
# Time-stamp: <2003-07-04 12:48:55 jvanheld>
#
############################################################
#use strict;;

=pod

=head1 NAME

genome-blast.pl

=head1 DESCRIPTION

This programs reads the BLAST output from a genome-to-genome
comparison, and creates data tables as well as the SQL scripts for
creating and loading the tables in a SQL database.

=head1 CATEGORY

util

=head1 USAGE
    
genome-blast.pl -q query_organism -db db_organism [-i inputfile] [-o outputfile] [-v #]

=head1 INPUT FORMAT

=head1 OUTPUT FORMAT


=cut


BEGIN {
    if ($0 =~ /([^(\/)]+)$/) {
	push (@INC, "$`lib/");
    }
}
require "RSA.lib";

################################################################
#### initialise parameters
my $start_time = &AlphaDate();

local %infile = ();
local %outfile = ();

local $verbose = 0;
local $in = STDIN;
local $out = STDOUT;

&ReadArguments();

################################################################
#### check argument values

unless ($query_organism) {
    &RSAT::error::FatalError("You should define a query organism");
}

unless ($db_organism) {
    &RSAT::error::FatalError("You should define a db organism");
}

unless ($infile{input}) {
    &RSAT::error::FatalError("You should define an input file");
}

################################################################
#### print verbose
&Verbose() if ($verbose);

################################################################
### open output stream
$out = &OpenOutputFile($outfile{output});

################################################################
##### read input
&RSAT::message::Info(join("\t", "Reading BLAST file", $infile{input})) if ($main::verbose >= 1);

($in) = &OpenInputFile($infile{input});
my $header = <$in>;
chomp($header);
$header =~ s/\r//;
my @header = split "\t", $header;
print $out join ("\t", "query_organism", "db_organism", @header, "rank"), "\n";
while (<$in>) {
    chomp();
    s/\r//;
    my ($query, 
	$level, 
	$subject, 
	$ident, 
	$ali_len, 
	$mismat, 
	$gap_opn, 
	$q_start, 
	$q_end,
	$s_start,
	$s_end,
	$e_value,
	$bit_sc) = split "\t";

    ## Index row per pair of sequence IDs
    $hit{$query}->{$subject} = $_;
    $score{$query}->{$subject} = $e_value;
    
#    print $out join ("\t",
#		     $query_organism,
#		     $db_organism,
#		     $query,
#		     $subject,
#		     $e_value
#		    ), "\n" if ($main::verbose >= 4);
		     
}
close $in if ($infile{input});

################################################################
###### execute the command
foreach my $query (sort keys %hit) {
    my $hit_ref = $hit{$query};
    my @hits = sort {$score{$query}->{$a} <=> $score{$query}->{$b}} keys (%$hit_ref);
    my $rank=0;
    foreach my $hit (@hits) {
	$rank++;
	print $out join("\t",
		       $query_organism,
		       $db_organism,
		       $hit{$query}->{$hit},
		       $rank), "\n";
    }
    &RSAT::message::Info(join ("\t", "Sorted hits for query", $query, scalar(@hits))) if ($main::verbose >= 3);
}

################################################################
###### print output


################################################################
###### finish verbose
if ($verbose >= 1) {
    my $done_time = &AlphaDate();
    print $out "; Job started $start_time\n";
    print $out "; Job done    $done_time\n";
}


################################################################
###### close output stream
close $out if ($outfile{output});


exit(0);


################################################################
################### subroutine definition ######################
################################################################


################################################################
#### display full help message 
sub PrintHelp {
    system "pod2text -c $0";
    exit()
}

################################################################
#### display short help message
sub PrintOptions {
    &PrintHelp();
}

################################################################
#### Read arguments 
sub ReadArguments {
    foreach my $a (0..$#ARGV) {

	## Verbosity
=pod


=head1 OPTIONS

=over 4

=item B<-v #>

Level of verbosity (detail in the warning messages during execution)

=cut
	if ($ARGV[$a] eq "-v") {
	    if (&IsNatural($ARGV[$a+1])) {
		$verbose = $ARGV[$a+1];
	    } else {
		$verbose = 1;
	    }
	    
	    ## Help message
=pod

=item B<-h>

Display full help message

=cut
	} elsif ($ARGV[$a] eq "-h") {
	    &PrintHelp();
	    
	    ## List of options
=pod

=item B<-help>

display options

=cut
	} elsif ($ARGV[$a] eq "-help") {
	    &PrintOptions();
	    

	    ## Input file
=pod

=item B<-i inputfile>

If no input file is specified, the standard input is used.  This
allows to use the command within a pipe.

The input file should be the result of the genome-to-genome BLAST,
obtained with the option blastall -m 8 (table output).

=cut
	} elsif ($ARGV[$a] eq "-i") {
	    $infile{input} = $ARGV[$a+1];
	    
	    ## Output file
=item	B<-o outputfile>

If no output file is specified, the standard output is used.  This
allows to use the command within a pipe.

=cut
	} elsif ($ARGV[$a] eq "-o") {
	    $outfile{output} = $ARGV[$a+1];
	    
	    ## Query organism
=pod

=item	B<-q query_organism>

Name of the query organism.

=cut
	} elsif ($ARGV[$a] eq "-q") {
	    $query_organism = $ARGV[$a+1];
	    
	    ## Db organism
=pod

=item	B<-db db_organism>

Name of the db organism.

=cut
	} elsif ($ARGV[$a] eq "-db") {
	    $db_organism = $ARGV[$a+1];
	    
	}

    }

=pod

=back

=cut

}

################################################################
#### verbose message
sub Verbose {
    print $out "; genome-blast.pl ";
    &PrintArguments($out);
    print $out ";\tQuery organism\t", $query_organism, "\n";
    print $out ";\tDB organism\t", $db_organism, "\n";

    if (defined(%infile)) {
	print $out "; Input files\n";
	while (($key,$value) = each %infile) {
	    print $out ";\t$key\t$value\n";
	}
    }
    if (defined(%outfile)) {
	print $out "; Output files\n";
	while (($key,$value) = each %outfile) {
	    print $out ";\t$key\t$value\n";
	}
    }
}


__END__

=pod

=head1 SEE ALSO

=cut
