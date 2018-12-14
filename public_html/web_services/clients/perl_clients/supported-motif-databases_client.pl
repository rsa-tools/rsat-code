#!/usr/bin/env perl

################################################################
## Connect a RSAT server and get a  list of supported motif dbs
##
## Usage:
##   perl supported-motif-db.pl [server_URL]

use strict;

## Import the modules we need for this test.  XML::Compile enables to
## parse the WSDL (description of the web services) on the flight. 
use XML::Compile::SOAP11;
use XML::Compile::WSDL11;
use XML::Compile::Transport::SOAPHTTP;

sub date {
  my ($sec, $min, $hour,$day,$month,$year) = localtime(time());
  my $my_date = sprintf("%02d-%02d-%02d.%02d%02d%02d", 1900+$year,$month+1,$day,$hour, $min, $sec);
  return $my_date;
}

## Specification of the server(s)
my @servers =  $ARGV[0] || qw(
        http://rsat.ulb.ac.be/rsat
    	http://prokaroyte.rsat.eu
	http://fungi.rsat.eu
	http://plants.rsat.eu
	http://metazoa.rsat.eu
	http://protists.rsat.eu
	http://teaching.rsat.eu
	http://training.rsat.eu);

my @more_servers = qw(
        http://pedagogix-tagc.univ-mrs.fr/rsat
        http://rsat-tagc.univ-mrs.fr/rsat
        http://rsat.sb-roscoff.fr
        http://www.rsat.eu
        http://www.rsat.fr
        http://embnet.ccg.unam.mx/rsa-tools
        http://rsat01.biologie.ens.fr/rsa-tools
        http://tagc.univ-mrs.fr/rsa-tools
        http://localhost/rsat
        );

## Query parameters
my $return = 'name,format,file,descr,version,url';
my %args = (
  'return'=>$return
    );


warn "Getting lists of supported motif databases from server(s)\n\t", join("\n\t", @servers), "\n\n";

foreach my $server (@servers) {
  warn "\n\n", "Querying server\t", $server, "\n";
  eval
  {
    # Retrieving and processing the WSDL
    my $wsdl_url = $server.'/web_services/RSATWS.wsdl';
    warn (&date(), "\t", "Parsing Web service description from WSDL", "\t", $wsdl_url, "\n");
    my $wsdl  = XML::LibXML->new->parse_file($wsdl_url);
    my $proxy = XML::Compile::WSDL11->new($wsdl);
    
    ## Compiling the client for supported-motif-databases
    warn (&date(), "\t", "Compiling client\n");
    my $client = $proxy->compileClient('supported_motif_databases');
    
    # Calling the service and getting the response
    warn (&date(), "\t", "Sending query to server", "\t", $server, "\n");
    my $answer = $client->( request => {%args});
    
    my $file = "databases_".$server.".txt";
    $file =~ s|http://||;
    $file =~ s|/|_|g;
    open OUT, ">$file";
    warn (&date(), "\t", "Result stored in file\t", $file, "\n");

    ## Write answer to output file
    if ( defined $answer ) {
      warn (&date(), "\t", "Server command : ".$answer->{output}->{response}->{command}."\n");
      print OUT "; Server : ", $server, "\n";
      print OUT "; WSDL URL : ", $wsdl_url, "\n";
      print OUT "; Server command : ".$answer->{output}->{response}->{command}."\n";
      print OUT "; Server file : ".$answer->{output}->{response}->{server}."\n";
      print OUT $answer->{output}->{response}->{client}."\n";
    } else {
      print OUT "No answer\n";
    }

    close OUT;
    my $nb_organisms = `grep -v '^;' ${file} | grep -v '^#' | wc -l`;
    chomp($nb_organisms);
    print join("\t", $nb_organisms, "databases at", $server, $file), "\n";
  };

  ## Report exceptions
  if ($@) {
    warn (&date(), "\t", "Caught an exception\n");
    warn ($@."\n");
    print OUT (&date(), "\t", "Caught an exception\n");
    print OUT ($@."\n");
  }  
}
