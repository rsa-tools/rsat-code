#!/usr/bin/perl

################################################################
## Connect a RSAT server and get a  list of supported organisms
##
## Usage:
##   perl supported-organisms_client_nostubb.wsdl [server_URL]

use strict;

# import the modules we need for this test; XML::Compile is included on the server
# by default.
use XML::Compile::SOAP11;
use XML::Compile::WSDL11;
use XML::Compile::Transport::SOAPHTTP;

## Specification of the server
#my $server = $ARGV[0] || "http://rsat.bigre.ulb.ac.be/rsat";
my @servers =  $ARGV[0] || qw(
			    http://rsat.ulb.ac.be/rsat
			    http://139.124.66.4/rsat
                            http://rsat.sb-roscoff.fr
			    http://www.rsat.eu
			    http://www.rsat.fr
			    http://embnet.ccg.unam.mx/rsa-tools
			    http://rsat01.biologie.ens.fr/rsa-tools
			    http://tagc.univ-mrs.fr/rsa-tools
			    http://anjie.bi.up.ac.za/rsa-tools
			    http://bongcam1.hgen.slu.se/rsat
			    http://localhost/rsat
			    );

#			    http://mamaze.ulb.ac.be/rsat
# 			    http://wwwsup.scmbb.ulb.ac.be/rsat

## Query parameters
my $taxon = 'Bacteria';
my $return = 'ID,taxonomy';
my $format = "tab";
my $depth = 5;
my %args = (
  'return'=>$return,
  'format'=>$format,
  'taxon' => $taxon,
  'depth'=>$depth,
    );


warn "Getting lists of supported organisms from server(s)\n\t", join("\n\t", @servers), "\n\n";

foreach my $server (@servers) {
  warn "\n\n", "Querying server\t", $server, "\n";
  eval
  {
    # Retrieving and processing the WSDL
    my $wsdl_url = $server.'/web_services/RSATWS.wsdl';
    warn ("Parsing Web service description from WSDL", "\t", $wsdl_url, "\n");
    my $wsdl  = XML::LibXML->new->parse_file($wsdl_url);
    my $proxy = XML::Compile::WSDL11->new($wsdl);
    
    ## Compiling the client for supported-organisms
    warn ("Compiling client\n");
    my $client = $proxy->compileClient('supported_organisms');
    
    # Calling the service and getting the response
    warn ("Sending query to server", "\t", $server, "\n");
#      warn "Getting list of supported organisms from server\t", $server, "\n";
    my $answer = $client->( request => {%args});
    #    print OUT "Answer: ".$answer."\n";
    
    my $file = "organisms_".$server.".txt";
    $file =~ s|http://||;
    $file =~ s|/|_|g;
    open OUT, ">$file";
    warn "Result stored in file\t", $file, "\n";

    ## Open output file
    # If the response arrived, look for a specific pattern
    # If the pattern is present, return 0 because the test passed.
    # If the result is something else, return 2 to indicate a warning.
    # If no answer has arrived, return 1 to indicate the test failed.
    if ( defined $answer ) {
      warn ("Server command : ".$answer->{output}->{response}->{command}."\n");
      print OUT "; Server : ", $server, "\n";
      print OUT "; WSDL URL : ", $wsdl_url, "\n";
      print OUT "; Server command : ".$answer->{output}->{response}->{command}."\n";
      print OUT "; Server file : ".$answer->{output}->{response}->{server}."\n";
      print OUT $answer->{output}->{response}->{client}."\n";
    } else {
      print OUT "No answer\n";
    }

    close OUT;
  };

  ## Report exceptions
  if ($@) {
    warn "Caught an exception\n";
    warn $@."\n";
    print OUT "Caught an exception\n";
    print OUT $@."\n";
  }  
}
