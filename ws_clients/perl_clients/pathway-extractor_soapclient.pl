#!/usr/bin/perl -w

# -- SOAP::Lite -- guide.soaplite.com -- Copyright (C) 2001 Paul Kulchenko --

#     use SOAP::Lite +autodispatch =>
#       uri => 'http://localhost/PathwayExtraction',
#       proxy => 'http://localhost/rsat/web_services/RSATWM.cgi'
#       ,
#       on_fault => sub { my($soap, $res) = @_; 
#       die ref $res ? $res->faultdetail : $soap->transport->status, "\n"};
# #  my $results = SOAP->inferpathway();
# my $results = SOAP->hi();
#   print "return :".$results."\n";

  use SOAP::Lite 
    on_action => sub {
       sprintf '%s#%s',@_
#       print "Message:".$_[0]."\n";
    },
    on_fault => sub { 
	my($soap, $res) = @_; 
#      die ref $res ? $res->faultstring : $soap->transport->status, "\n";
	print "ERROR: ". $soap->transport->status . "\n"
    };
    
  my $soap = SOAP::Lite
    -> uri('http://localhost/PathwayExtractor_WS')
    -> proxy('http://localhost/rsat/web_services/RSATWM.cgi');

  eval {  
    my $results=  $soap->inferpathway("NP_416523.1\tNP_416524.1\tNP_416525.1\tNP_416526.4\tNP_416527.1\tNP_416528.2\tNP_416529.1\tNP_416530.1"); 
    print "Results:" . $results->result ."\n"; 
  1 } or die;