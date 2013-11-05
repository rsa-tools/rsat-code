#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
require "RSA.lib";

&UpdateCounterFile;
&UpdateLogFile;

print "Content-type: text/plain\n\n";
print `cat $counter_file`;
print "\n";

exit(0);

