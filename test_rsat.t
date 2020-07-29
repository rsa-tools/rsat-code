use strict;
use warnings;
use Test::More tests => 1;

# protptype with only one test, should be fully fleshed out and 
# cover all RSAT scripts

ok( eval{ `perl perl-scripts/random-seq -l 10 -seed 1234 ` } =~ /GACCTCAGGC/ , 'random-seq' );
