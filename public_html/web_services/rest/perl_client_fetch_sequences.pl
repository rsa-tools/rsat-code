use REST::Client;
use JSON qw/encode_json decode_json/;
use HTTP::Request::Common;

my $rsat_server = "http://rsat-tagc.univ-mrs.fr/rest/fetch-sequences";
my $client = REST::Client->new();

########## EXAMPLE FOR SENDING JSON DATA ##########
#my %args = ();
#$args{"genome"} = "mm9";
#$args{"url"} = "http://rsat-tagc.univ-mrs.fr/rsat//demo_files/fetch-sequences_Schmidt_2011_mm9_CEBPA_SWEMBL_R0.12_702peaks.bed";
#$args{"header_format"} = "galaxy";

#my $arg = encode_json(\%args);
#$client->POST($rsat_server."/fetch_sequences", $arg, {"Content-type"=>"application/json"});
#my $ret = decode_json($client->responseContent);
#print "Result file URL: " .$ret->{'server'} . "\n";
#################### END EXAMPLE

########## EXAMPLE FOR SENDING FOMR-DATA INCLUDING FILES AND OTHER TYPES OF DATA #########
my $request = HTTP::Request::Common::POST('', 'Content-Type'=>'form-data', 'Content'=> ['tmp_input_file'=>["/home/rsat/rsat/public_html/demo_files/fetch-sequences_Schmidt_2011_mm9_CEBPA_SWEMBL_R0.12_702peaks.bed"], 'genome'=>'mm9', 'header_format'=>'galaxy']);

my $headers = {'Content-type'=> $request->header('Content-Type')};
my $body_content = $request->content;

$client->POST($rsat_server, $body_content, $headers);
my $ret = decode_json($client->responseContent);
print "Result file URL: " .$ret->{'server'} . "\n";
############## END EXAMPLE
