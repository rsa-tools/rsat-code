#!/usr/bin/perl -w
BEGIN {
  unless ($ENV{RSAT}) {
      $ENV{RSAT} = $0; $ENV{RSAT} =~ s|/public_html/+web_services/.*||; ## Guess RSAT path from module full name
  }
  if ($0 =~ /([^(\/)]+)$/) {
	push (@INC, "$`lib/");
 	push (@INC,"$ENV{RSAT}/perl-scripts/lib/");
 	push (@INC,"$ENV{RSAT}/perl-scripts/");
 	push (@INC,"$ENV{RSAT}/public_html/web_services/modules/");
  }
}
  use SOAP::Transport::HTTP;
#   use lib '../../perl-scripts/lib/RSAT';
#   use lib '../../public_html/web_services/modules';


  my $server = SOAP::Transport::HTTP::CGI
#     -> dispatch_to('$ENV{RSAT}/public_html/web_services/modules')#
#        -> dispatch_to('PathwayExtractor_WS')
#      -> dispatch_to('modules/PathwayExtractor_WS')
       -> dispatch_to('RSATWM')

    -> handle;
