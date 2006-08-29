#!/usr/bin/perl
# RSATWS.cgi - SOAP server for rsa-tools.

use strict;

use SOAP::Transport::HTTP;
use lib '../../perl-scripts/lib/RSAT';

my $server = SOAP::Transport::HTTP::CGI
  -> dispatch_to('RSATWS')
  -> handle;

