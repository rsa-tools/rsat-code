#!/usr/bin/perl -w
############################################################
#
# $Id: genome-blast.pl,v 1.4 2005/04/25 09:35:08 jvanheld Exp $
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
push @INC, $RSA."/perl-scripts/parsers/" if ($RSA);
require "lib/load_classes.pl";
require RSAT::blast_hit;

################################################################
#### initialise parameters
my $start_time = &AlphaDate();

local %infile = ();
local %outfile = ();

local $verbose = 0;
local $in = STDIN;
local $out = STDOUT;

local @columns = qw(query level subject ident ali_len mismat gap_open q_start q_end s_start s_end e_value bit_sc);

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
#    $infile{input} = "q_".$query_organism."_db_".$db_organism.".tab";
    &RSAT::error::FatalError("You should define an input file");
}

#unless ($infile{recp}) {
#    $infile{recp} = "q_".$db_organism."_db_".$query_organism.".tab";
#}

################################################################
## Class factory for blast hits
my $blast_hits = classes::ClassFactory->new_class(object_type=>"RSAT::blast_hit",prefix=>"hit_");
    
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

my $h = 0;
while (<$in>) {
    chomp();
    s/\r//;
    next unless /\S/;
    $h++;
    my @fields = split "\t";

#     my ($query, 
# 	$level, 
# 	$subject, 
# 	$ident, 
# 	$ali_len, 
# 	$mismat, 
# 	$gap_opn, 
# 	$q_start, 
# 	$q_end,
# 	$s_start,
# 	$s_end,
# 	$e_value,
# 	$bit_sc) = split "\t";
#    ## Index row per pair of sequence IDs
#    $hit{$query}->{$subject} = $_;
#    $score{$query}->{$subject} = $e_value;

    ## Create a new object for the match
    my $hit = $blast_hits->new_object(id=>$query_organism."_".$db_organism."_".$h);
    foreach my $col (@columns) {
	$hit->set_attribute($col, shift @fields);
    }


    ## Index row per pair of sequence IDs
    my $query = $hit->get_attribute("query");
    my $subject = $hit->get_attribute("subject");
    
    push @{$hits_per_query{$query}}, $hit;
    push @{$hits_per_subject{$subject}}, $hit;
}
close $in if ($infile{input});


################################################################
## Calculate hit rank per query
&RSAT::message::Info("Ranking BLAST hits per query") if ($main::verbose >= 1);
foreach my $query (sort keys %hits_per_query) {
    my @sorted_hits = sort {$a->get_attribute("e_value") <=> $b->get_attribute("e_value") }  @{$hits_per_query{$query}};

    my $rank=0;
    foreach my $hit (@sorted_hits) {
	## Assign rank attribute
	$rank++;
#	&RSAT::message::Debug("Hit rank", $query, $rank, $hit) if ($main::verbose >= 10);
	$hit->set_attribute("q_rank", $rank);

	## Index best hits
	if ($rank == 1) {
	    $best_hit{$hit->get_attribute("query")} = $hit->get_attribute("subject");
	}

    }
    &RSAT::message::Info(join ("\t", "Sorted hits for query", $query, scalar(@hits))) if ($main::verbose >= 3);
}

################################################################
## Calculate hit rank per subject
&RSAT::message::Info("Ranking BLAST hits per subject") if ($main::verbose >= 1);
foreach my $query (sort keys %hits_per_subject) {
    my @sorted_hits = sort {$a->get_attribute("e_value") <=> $b->get_attribute("e_value") }  @{$hits_per_subject{$query}};

    my $rank=0;
    foreach my $hit (@sorted_hits) {
	## Assign rank attribute
	$rank++;
#	&RSAT::message::Debug("Hit rank", $query, $rank, $hit) if ($main::verbose >= 10);
	$hit->set_attribute("s_rank", $rank);

	## Index best hits
	if ($rank == 1) {
	    $best_hit{$hit->get_attribute("query")} = $hit->get_attribute("subject");
	}

    }
    &RSAT::message::Info(join ("\t", "Sorted hits for query", $query, scalar(@hits))) if ($main::verbose >= 3);
}


################################################################
###### print output
#my @header = split "\t", $header;
&RSAT::message::Info("Printing the result") if ($main::verbose >= 1);
my @header = join "\t", @columns;
print $out join ("\t", "query_organism", "db_organism", @header, "q_rank", "s_rank"), "\n";
foreach my $hit ($blast_hits->get_objects()) {
    my @fields = ();
    foreach my $col (@columns, "q_rank", "s_rank") {
	push @fields, $hit->get_attribute($col);
    }
    print $out join("\t",
		    $query_organism,
		    $db_organism,
		    @fields), "\n";
}


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
	    
	    
	    ## BLAST file
=pod

=item B<-i blast_file>

The input file should be the result of the genome-to-genome BLAST,
obtained with the option blastall -m 8 (table output). 

The input file is the result of a BLAST for all protein sequences of
the query organism against all protein sequences of the DB organism.

=cut
	} elsif ($ARGV[$a] eq "-i") {
	    $infile{input} = $ARGV[$a+1];
	    
=item	B<-o outputfile>

If no output file is specified, the standard output is used.  This
allows to use the command within a pipe.

=cut
	} elsif ($ARGV[$a] eq "-o") {
	    $outfile{output} = $ARGV[$a+1];
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
