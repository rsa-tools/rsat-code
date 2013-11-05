#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
require "cgi-lib.pl";








require "RSA.lib.pl";
$form_file = "$HTML/neighbour-orfs.html";

MAIN:
{
    ### Read the content of the form 
    &ReadParse(*input);

    $query = `cat $input{'query_file'}`;

    ### Print the header
    print &PrintHeader;

    $form = `cat $form_file`;

    $form =~ s/RSA.icon.gif/$ENV{rsat_www}\/RSA.icon.gif"/;   
    $form =~ s/lablogo.gif/$ENV{rsat_www}\/lablogo.gif"/;   
    $form =~ s/href="/href="$ENV{rsat_www}\//i unless ((/http/) || (/mailto/));   
    $form =~ s/(NAME="query"[^>]*>)/${1}$query/;   
    $form =~ s/ selected>ORF/>ORF/i;   
    if ($input{'query_format'} =~ /genome/i) { 
        $form =~ s/OPTION>(result of a genome search)/OPTION SELECTED>$1/i;   
    } elsif ($input{'query_format'} =~ /positions/i) {
        $form =~ s/OPTION>(chromosomal positions)/OPTION SELECTED>$1/i;   
    } elsif ($input{'query_format'} =~ /ORF/i) {
        $form =~ s/OPTION>(ORF identifiers)/OPTION SELECTED>$1/i;   
    }
    print $form;

}

