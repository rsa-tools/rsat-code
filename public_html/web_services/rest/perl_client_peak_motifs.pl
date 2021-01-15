use REST::Client;
use JSON qw/encode_json decode_json/;
use HTTP::Request::Common;

my $rsat_server = "http://rsat-tagc.univ-mrs.fr/rest/peak-motifs";
my $client = REST::Client->new();

my %args = ();
$args{"title"} = "Human_from_Jaspar";
$args{"markov"} = "auto";
$args{"disco"} = "oligos,positions";
$args{"motif_db"} = "jaspar_core_nonredundant_vertebrates";
$args{"tmp_test_infile_URL"} = "http://rsat-tagc.univ-mrs.fr/rsat//demo_files/ChIP-seq_peaks/Oct4_peaks_top1000.fa";
$args{"task"} = "purge,seqlen,composition,disco,merge_motifs,split_motifs,motifs_vs_motifs,timelog,archive,synthesis,small_summary,motifs_vs_db,scan";

$arg = encode_json(\%args);
$client->POST($rsat_server, $arg, {'Content-type' => 'application/json'});

print $client->responseContent;
#my $ret = decode_json($client->responseContent);
#print "Result file URL: " .$ret->{'server'} . "\n";
