use REST::Client;
use JSON qw(decode_json encode_json);
use LWP::Simple;

my $client = REST::Client->new();

my %args = ();
$args{"group"} = "Fungi";

my $arg = encode_json(\%args);
$client->POST("http://rsat-tagc.univ-mrs.fr/rest/supported-organisms", $arg, {'Content-type'=>"application/json"});
#my $ret = decode_json($client->responseContent());
print $client->responseContent();
#print "Result file URL: " . $ret->{'server'} . "\n";
