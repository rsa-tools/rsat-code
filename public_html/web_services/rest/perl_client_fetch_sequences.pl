use JSON qw/encode_json decode_json/;
use HTTP::Request::Common;

use REST::Client;

my $rsat_server = "http://rsat-tagc.univ-mrs.fr/rest/fetch-sequences/";
my $client = REST::Client->new();

############## Example for sending JSON data
my %args = ();
$args{"genome"} = "mm9";
$args{"i_string"} = "http://rsat-tagc.univ-mrs.fr/rsat/demo_files/fetch-sequences_Schmidt_2011_mm9_CEBPA_SWEMBL_R0.12_702peaks.bed";
$args{"i_string_type"} = 'url';
$args{"header_format"} = 'galaxy';

my $arg =  encode_json(\%args);
$client->POST($rsat_server, $arg, {"Content-type"=>"application/json"});
#$client->GET($rsat_server.'?content-type=application/json&genome=mm9&header_format=galaxy&i_string=http://rsat-tagc.univ-mrs.fr/rsat/demo_files/fetch-sequences_Schmidt_2011_mm9_CEBPA_SWEMBL_R0.12_702peaks.bed&i_string_type=url');
my $ret = $client->responseContent;
print $ret;

############# Example for sending form-data including files
print "\nExample for sending form-data, return text/plain...\n";

my $req = HTTP::Request::Common::POST('', 'Content-type'=>'form-data', 'Content' => ['i'=>['/workspace/rsat/public_html/demo_files/fetch-sequences_Schmidt_2011_mm9_CEBPA_SWEMBL_R0.12_702peaks.bed'], 'genome'=>'mm9', 'header_format'=>'galaxy']);

my $body_content = $req->content;
my $headers = {'Content-type'=>$req->header('Content-type'), 'Accept'=>'text/plain'};
$client->POST($rsat_server, $body_content, $headers);
print $client->responseContent;
