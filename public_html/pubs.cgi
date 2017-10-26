#!/usr/bin/perl
#### this cgi script fills the HTML form for the publications page
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}

use CGI;
use CGI::Carp qw/fatalsToBrowser/;

require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

$query = new CGI;

print "Content-type: text/plain\n\n";


$default{year} = '2017';
$default{type} = 'all';

foreach $key (keys %default) {
    if ($query->param($key)) {
        $default{$key} = $query->param($key);
    }
}

my $pub_file = "publications.csv";
open(my $data, '<', $pub_file) or die "Could not open '$file' $!\n";
my %years;
while( my $line = <$data> ){
    chomp $line;
    
    my @fields = split ";", $line;
    if($years{$fields[0]}){
        push @{$years{$fields[0]}}, [$fields[1] . '. ' . $fields[2], $fields[3]];
    }else{
        $years{$fields[0]} = ();
        push @{$years{$fields[0]}}, [$fields[1] . '. ' . $fields[2], $fields[3]];
    }
}

foreach $key (sort {$b <=> $a} keys %years){
    foreach $_ (@{$years{$key}}){
        if(($default{type} ne 'all' && ${$_}[1] eq $default{type}) || ($default{type} eq 'all')){
            $hasPub = 'true';
        }
    }
    if($hasPub eq 'true' && ( ($default{year} ne 'All' && $key eq $default{year}) || ($default{year} eq 'All') )){
        print "<h3>$key</h3>";
        print "<ol>";
        foreach $_ (@{$years{$key}}){
            if(($default{type} ne 'all' && ${$_}[1] eq $default{type}) || ($default{type} eq 'all')){
                print "<p><li>";
                print ${$_}[0];
                print "</li></p>";
            }
        }
        print "</ol>";
        
    }
    $hasPub = '';
}
