#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
require "RSA.lib.pl";
use strict;

my $out_format = $ARGV[0];

print &ListSupportedOrganisms($out_format);

exit(0);

sub ListSupportedOrganisms {
  ### usage : &ListSupportedOrganisms($format);
  ### $supp =  &ListSupportedOrganisms("text");
  ### $supp =  &ListSupportedOrganisms("html_list");
  ### $supp =  &ListSupportedOrganisms("html_table");
  ### @names =  &ListSupportedOrganisms("array");
  my ($out_format) = @_;

  if ($out_format eq "html_list") {
    my $result = "<UL>\n";
    foreach my $organism (sort keys %main::supported_organism) {
      $result .= "<LI>";
      $result .= $main::supported_organism{$organism}->{name};
      $result .= "\n";
    }
    $result .= "</UL>\n";
    return $result;
  } elsif ($out_format eq "html_table") {
    my $result = "<TABLE>\n";
    foreach my $organism (sort keys %main::supported_organism) {
      $result .= "<TR>\n";
      $result .= "<TD>$organism</TD>\n";
      $result .= "<TD>";
      $result .= $main::supported_organism{$organism}->{name};
      $result .= "</TD>\n";
      $result .= "</TR>\n";
    }
    $result .= "</TABLE>\n";
    return $result;
  } elsif ($out_format eq "array") {
    my @result = ();
    foreach my $organism (sort keys %main::supported_organism) {
      push @result, $main::supported_organism{$organism}->{name};
    }
    return @result;
  }else {
    my $result = "";
    foreach my $organism (sort keys %main::supported_organism) {
      $result .= $organism;
      $result .= "\t";
      $result .= $main::supported_organism{$organism}->{name};
      $result .= "\n";
    }
    return $result;
  }
}
