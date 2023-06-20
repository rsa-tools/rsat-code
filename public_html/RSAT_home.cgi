#!/usr/bin/env perl

# redirects old, out-of-date home CGI to current index.php

use strict;
use warnings;
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
require "RSA.lib";

print "Content-type: text/html\n\n";
print <<"HTML";
<!DOCTYPE html>
<html>
  <head>
    <meta http-equiv="refresh" content="0; url='$ENV{rsat_www}'" />
  </head>
  <body>
    <p>redirecting to $ENV{rsat_www}</p>
  </body>
</html>
HTML
