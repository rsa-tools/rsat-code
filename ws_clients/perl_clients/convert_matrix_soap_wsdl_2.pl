#!/usr/bin/perl
# the first line of the script must tell us which language interpreter to use,
# in this case its perl

use strict;
use SOAP::WSDL; ## Requires version 2.0 or later of SOAP::WSDL
use lib 'RSATWS';
use MyInterfaces::RSATWebServices::RSATWSPortType;

## WSDL location
my $server = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services';

## Service call
my $soap=MyInterfaces::RSATWebServices::RSATWSPortType->new();

#Defining a few parameters
my $features = '#seq_id	ft_type	ft_name	strand	start	end	sequence	weight	proba_M	proba_B	Pval	ln_Pval	sig
GSPATG00028047001	site	matscan-matrix.1	D	85	94	TTACACGTTT	8.68.4e-03	1.6e-06	2.8e-05	-10.481	4.552
GSPATG00028047001	site	matscan-matrix.3	R	86	96	AAAAACGTGTA	10.3.6e-02	6.5e-07	4.2e-06	-12.379	5.376
GSPATG00028047001	site	matscan-matrix.2	R	127	136	AAAAACGGAT	10.6.7e-02	1.6e-06	5.2e-06	-12.162	5.282
GSPATG00032430001	site	matscan-matrix.2	R	137	146	AAAAACGAAT	10.1.8e-01	6.6e-06	1.2e-05	-11.349	4.929
GSPATG00036447001	site	matscan-matrix.1	D	37	46	TTACACGTAT	9.93.1e-02	1.6e-06	1.1e-05	-11.396	4.949
GSPATG00036447001	site	matscan-matrix.3	D	37	47	TTACACGTATA	10.2.2e-02	6.5e-07	5.8e-06	-12.059	5.237
GSPATG00036447001	site	matscan-matrix.3	R	38	48	ATATACGTGTA	8.94.3e-03	6.5e-07	2.5e-05	-10.603	4.605
GSPATG00009610001	site	matscan-matrix.1	D	185	194	TTACACGTAT	9.93.1e-02	1.6e-06	1.1e-05	-11.396	4.949
GSPATG00009610001	site	matscan-matrix.3	R	186	196	TAATACGTGTA	9.68.7e-03	6.5e-07	1.0e-05	-11.503	4.996
GSPATG00015555001	site	matscan-matrix.2	D	167	176	AAAAACGAAT	10.1.8e-01	6.6e-06	1.2e-05	-11.349	4.929
GSPATG00003426001	site	matscan-matrix.1	D	173	182	TAACACGTTT	9.11.5e-02	1.6e-06	1.9e-05	-10.858	4.715
GSPATG00003426001	site	matscan-matrix.3	R	174	184	TAAAACGTGTT	7.39.5e-04	6.5e-07	7.0e-05	-9.571	4.157
GSPATG00010335001	site	matscan-matrix.1	D	9	18	TTACACGTTT	8.68.4e-03	1.6e-06	2.8e-05	-10.481	4.552';

my $return = 'counts,parameters';
my $from = 'feature';
my $to = 'tab';

my %args = (
	'matrix' => $features,
	'return' => $return,
	'from' => $from,
	'to' => $to
	);

## Send the request to the server
print "Sending request to the server $server\n";
my $som = $soap->convert_matrix({'request' => \%args});

## Get the result
unless ($som) {
	printf "A fault (%s) occured: %s\n", $som->get_faultcode(), $som->get_faultstring();
} else {
	my $results = $som->get_response();

	## Report the remote command
	my $command = $results -> get_command();
	print "Command used on the server: ".$command, "\n";

	## Report the result
	my $server_file = $results -> get_server();
	my $result = $results -> get_client();
	print "Result file on the server: ".$server_file."\n";
	print "Retrieved sequence(s): \n".$result;
}
