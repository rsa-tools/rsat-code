use REST::Client;
use JSON qw/encode_json decode_json/;
use HTTP::Request::Common;

my $rsat_server = "http://rsat-tagc.univ-mrs.fr/rest";
my $client = REST::Client->new();

########## EXAMPLE FOR SENDING JSON DATA ##########
#my %args = ();
#$args{"genome"} = "mm9";
#$args{"u"} = "http://rsat-tagc.univ-mrs.fr/rsat//demo_files/fetch-sequences_Schmidt_2011_mm9_CEBPA_SWEMBL_R0.12_702peaks.bed";
#$args{"header_format"} = "galaxy";
#$args{"content_type"} = 'text';
#
#my $arg = encode_json(\%args);
#$client->POST($rsat_server."/fetch-sequences", $arg, {"Content-type"=>"application/json"});
#my $ret = $client->responseContent;
#print "Result file URL: " .$ret->{'server'} . "\n";
#print $ret;
#################### END EXAMPLE

########## EXAMPLE FOR SENDING FOMR-DATA INCLUDING FILES AND OTHER TYPES OF DATA #########
#my $request = HTTP::Request::Common::POST('', 'Content-Type'=>'form-data', 'Content'=> ['i'=>["/home/rsat/rsat/public_html/demo_files/fetch-sequences_Schmidt_2011_mm9_CEBPA_SWEMBL_R0.12_702peaks.bed"], 'genome'=>'mm9', 'header_format'=>'galaxy']);
#
#my $headers = {'Content-type'=> $request->header('Content-Type')};
#my $body_content = $request->content;
#
#$client->POST($rsat_server.'/fetch-sequences', $body_content, $headers);
#print $client->responseContent
############## END EXAMPLE

$client->setHost($rsat_server);
my $headers = {"Accept" => 'text/plan'};
my $ret = $client->GET('/fetch-sequences?genome=mm9&header_format=galaxy&u=http://rsat-tagc.univ-mrs.fr/rsat//demo_files/fetch-sequences_Schmidt_2011_mm9_CEBPA_SWEMBL_R0.12_702peaks.bed', $headers)->responseContent();
print $ret;
