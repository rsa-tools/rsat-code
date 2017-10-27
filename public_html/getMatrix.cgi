#!/usr/bin/perl
############################################################
#
# $Id: getMatrix.cgi: get the matrix from the motif name
#
############################################################
#### this cgi script fills the HTML form for the program dna-pattern
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use JSON;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";

$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

$query = new CGI;

my $mode = $query->param("mode");
my $output = $query->param("output");

if($mode eq 'retrieve'){
    print "Content-type:application/json; charset=iso-8859-1\n\n";
}else{
    print "Content-type:txt/plain\n\n";
}

my $db_choice = $query->param("db_choice");
my @db_id = split(",", $query->param("db_id"));

my %matrix_db = &RSAT::server::supported_motif_databases();
my %db = %{$matrix_db{$db_choice}};
my $format = $db{format};
my $file = $db{file};

if($mode ne "retrieve" && $output ne "email"){
    print "$format</format>";
}

##### execute la commande retrieve-matrix
$file_name = $ENV{RSAT} . "/public_html/motif_databases/" . $file;
$id = $query->param("db_id");
$command = "$SCRIPTS/retrieve-matrix -i $file_name -id $id 2> /dev/null";
open RESULT, "$command |";

######## return json array of all information

    my @results;
    my %result;
    my $result_cc = "";
    my $result_all = "";
    my $mirror = 0;

    $tmp_file_path = &RSAT::util::make_temp_file("", "retrieve-matrix", 1,0);     
    $result_file = $tmp_file_path . "." . lc($format);
    if(open MIRROR, ">$result_file"){
        $mirror = 1;
        &DelayedRemoval($result_file);
    }

if($mode eq "retrieve"){
    while(my $row = <RESULT>){
        print MIRROR $row if ($mirror);

        if($row =~ /^AC\s+/){
            if(%result){
                if(! ($result_all =~ /\/\//) ){
                    $result_all .= "//";
                }
                if($result{ID} eq ""){
                    $result{ID} = "\n";
                }
                if($result{DE} eq ""){
                    $result{DE} = "\n";
                }
                if($result{CC} eq ""){
                    $result{CC} = "\n";
                }
                push @results, {"info" => "<b>Accesssion:</b> $result{AC}<br/><b>Identifier:</b> $result{ID}<br/><b>Description:</b> $result{DE}<br/><b>Additional information:</b> $result_cc", "all" => $result_all};
            }
            undef %result;
            $result_cc = "";
            $result_all = "";
            @f = split(/^AC\s+/, $row); 
            $result{AC} = $f[1];            
        }elsif($row =~ /^ID\s+/){
            @f = split(/^ID\s+/, $row);
            $result{ID} = $f[1];
        }elsif($row =~ /^DE\s+/){
            @f = split(/^DE\s+/, $row);
            $result{DE} = $f[1];
        }elsif($row =~ /^CC\s+/){
            @f = split(/^CC\s+/, $row);
            $result_cc .= "&nbsp;&nbsp;" . $f[1] . "<br/>";
        }        
        $result_all .= $row ;  
    }

    if(! ($result_all =~ /\/\//) ){
        $result_all .= "//";
    }
    if($result{ID} eq ""){
        $result{ID} = "<br/>";
    }
    if($result{DE} eq ""){
        $result{DE} = "<br/>";
    }
    if($result{CC} eq ""){
        $result{CC} = "<br/>";
    }
    push @results, {"info" => "<b>Accesssion:</b> $result{AC}<br/><b>Identifier:</b> $result{ID}<br/><b>Description:</b> $result{DE}<br/><b>Additional information:</b> $result_cc", "all" => $result_all};
    my $result_URL = $ENV{rsat_www}."/tmp/";
    $result_URL .= &RSAT::util::RelativePath(&RSAT::util::get_pub_temp(), $result_file);
    my $input_URL = $ENV{rsat_www}."/motif_databases/$file";
    print JSON::encode_json({entries => \@results, resultfile => $result_URL, inputfile => $input_URL});
}else{
    if($output eq "email"){
        while(my $row = <RESULT>){
            print MIRROR $row if ($mirror);
        }
        &EmailTheResult("$command", $query->param('user_email'), $result_file);
    }else{
        while(<RESULT>){
            print $_;
        }
    }
}
close RESULT;
exit(0);
