#!/usr/bin/perl

unless ($ENV{RSAT}) {
	die "ERROR: you should define the environment variable RSAT (main dir of the rsa-tools)";
}

require $ENV{RSAT}."/RSA.config";

return 1;
