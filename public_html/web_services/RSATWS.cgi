#!/usr/bin/perl
# RSATWS.cgi - SOAP server for rsa-tools.

use strict;

use SOAP::Transport::HTTP;
#use lib '.';

my $RSAT = $0; $RSAT =~ s|/public_html/+web_services/.*||;
my $SCRIPTS = $RSAT.'/perl-scripts';
my $TMP = $RSAT.'/public_html/tmp';

my $server = SOAP::Transport::HTTP::CGI
  -> dispatch_to('RSATWS')
  -> handle;


package RSATWS;

sub retrieve_seq {
    my @input = @_;
    my $organism = $input[1]{"organism"};
    my $noorf = $input[1]{"noorf"};
    my $from = $input[1]{"from"};
    my $to = $input[1]{"to"};
    my $query = '';
    if ($input[1]{"query"}){
	my $i = 0;
	while ($input[1]{"query"}[$i]){
	    $query .= "-q ".$input[1]{"query"}[$i]." ";
	    $i++;
	}
    }
    my $feattype = $input[1]{"feattype"};
    my $type = $input[1]{"type"};
    my $format = $input[1]{"format"};
    my $all = $input[1]{"all"};
    my $n = $input[1]{"number"};
    my $label = $input[1]{"label"};
    my $label_sep = $input[1]{"label_sep"};
    my $nocom = $input[1]{"nocom"};
    my $result = `$SCRIPTS/retrieve-seq $organism $query -lw 0 $noorf $from $to $feattype $type $format $all $n $label $label_sep $nocom`;
    my $tmp_outfile = `mktemp $TMP/retrieve-seq.XXXXXXXXXX`;
    open TMP, ">".$tmp_outfile or die "cannot open temp file ".$tmp_outfile."\n";
    print TMP $result;
    close TMP;
#    return $tmp_outfile;
#    return $result;
    my @resultat = ($tmp_outfile,$result);
    return @resultat;
}


sub purge_seq {
    my @input = @_;
    my $sequence = $input[1]{"sequence"};
#    my $tmp_infile = $input[1]{"tmp_infile"};
    my $tmp_infile = `mktemp $TMP/purge-seq.XXXXXXXXXX`;
    chomp $tmp_infile;
    open TMP_IN, ">".$tmp_infile or die "cannot open temp file ".$tmp_infile."\n";
    print TMP_IN $sequence;
    close TMP_IN;
    my $format = $input[1]{"format"};
    my $match_length = $input[1]{"match_length"};
    my $mismatch = $input[1]{"mismatch"};
    my $str = $input[1]{"str"};
    my $delete = $input[1]{"delete"};
    my $mask_short = $input[1]{"mask_short"};
    my $result = `$SCRIPTS/purge-sequence -i $tmp_infile $format $match_length $mismatch $str $delete $mask_short`;
    my $tmp_outfile = `mktemp $TMP/purge-seq.XXXXXXXXXX`;
    open TMP_OUT, ">".$tmp_outfile or die "cannot open temp file ".$tmp_outfile."\n";
    print TMP_OUT $result;
    close TMP_OUT;
#    return $tmp_outfile;
#    return $result;
    my @resultat = ($tmp_outfile,$result);
    return @resultat;
}


sub oligo_analysis {
    my @input = @_;
    my $sequence = $input[1]{"sequence"};
    my $format = $input[1]{"format"};
    my $length = $input[1]{"length"};
    my $organism = $input[1]{"organism"};
    my $background = $input[1]{"background"};
    my $stats = $input[1]{"stats"};
    my $noov = $input[1]{"noov"};
    my $str = $input[1]{"str"};
    my $sort = $input[1]{"sort"};
    my $lth = $input[1]{"lth"};
    my $tmp_infile = $input[1]{"tmp_infile"};
#    my $tmp_infile = `mktemp $TMP/oligo_analysis.XXXXXXXXXX`;
    chomp $tmp_infile;
#    open TMP, ">".$tmp_infile or die "cannot open temp file ".$tmp_infile."\n";
#    print TMP $sequence;
#    close TMP;
    my $result = `$SCRIPTS/oligo-analysis -i $tmp_infile $format $length $organism $background $stats $noov $str $sort $lth`;
    return $result;
}

1;
