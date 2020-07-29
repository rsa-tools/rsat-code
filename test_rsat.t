use strict;
use warnings;
use Test::More tests => 1;

ok( eval{ `perl perl-scripts/random-seq -l 10 -seed 1234 ` } =~ /GACCTCAGGC/ , 'random-seq' );

