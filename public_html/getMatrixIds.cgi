#!/usr/bin/env perl
############################################################
#
# $Id: getMatrixIds.cgi: get the identifiants from the motif file
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

if($mode eq 'retrieveweb'){
    print "Content-type:application/json; charset=iso-8859-1\n\n";
}else{
    print "Content-type:txt/plain\n\n";
}

my @dbs_choice = split(",", $query->param("dbs_choice"));
my $db_choice = $query->param("db_choice") || "";
my $mode = $query->param("mode") || "";

my $motif_folder = "/public_html/motif_databases/";

my %matrix_db = &RSAT::server::supported_motif_databases_2();
my @dbs;
my @ids;

my %colors = (
    'Vertebrates' => 'blue',
    'Non-vertebrate Metazoa' => 'orange',
    'Plants' => 'green',
    'Fungi' => 'purple',
    'Prokaryotes' => 'teal',
    'Multi-organisms' => 'yellow',
    'RNA binding' => 'coral'
);

if($mode eq "db"){
    my %categories = ();
    foreach my $db_name (@dbs_choice){
        foreach my $key (keys %{$matrix_db{$db_name}}){
            my %db = %{$matrix_db{$db_name}{$key}};
            if(! defined $categories{$db{category}}){ $categories{$db{category}} = (); }
            push @{$categories{$db{category}}}, {"descr" => $db{descr}. " (" . $db{version} . ")", "name" => $key, };
        }
    }
    foreach my $key (sort sort {"\L$a" cmp "\L$b"} (keys %categories)){
        if($key ne ""){
            push @dbs, {"category" =>$key, "data" => \@{$categories{$key}}, "color" => $colors{$key} };
        }
    }
    print JSON::encode_json( {entries => \@dbs} );
}elsif($mode eq "id"){
    my $dbs_name = $dbs_choice[0];
    my %db = %{$matrix_db{$dbs_name}{$db_choice}};
    
    my $file = $db{file};
    
    $file_name = $ENV{RSAT} . $motif_folder . $file;
    open($fh, "<", $file_name) or die "Cannot open file $!";
    my $id = "";
    while(my $row = <$fh>){
        if($row =~ /^AC\s+/){
            my @f = split(/\s+/, $row);
            $id = $f[1];
        }
        if($row =~ /^ID\s+/){
            my @f = split(/\s+/, $row);
            my $idac = "";
            if($f[1] ne $id){
                $idac = $f[1] . " - " . $id;
            }else{
                $idac = $id;
            }
            push @ids, { "id" => $id, "idac" => $idac };
        }
        if($row =~ /^OS\s+/){
            push @ids, { "id" => $id, "idac" => $id };
        }
    }
    print JSON::encode_json( {entries => \@ids} );
    
}elsif($mode =~ /^retrieve/){
    my $dbs_name = $dbs_choice[0];
    my %db = %{$matrix_db{$dbs_name}{$db_choice}};
    my $format = $db{format};
    my $file = $db{file};
    
    if($mode eq "retrievepipe"){
        print "$format</format>";
    }
    ##### execute la commande retrieve-matrix
    $file_name = $ENV{RSAT} . $motif_folder . $file;
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
    
    if($mode eq "retrieveweb"){
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
            $result_all .= $row . "<br/>";
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
        if($mode eq "retrieveemail"){
            while(my $row = <RESULT>){
                print MIRROR $row if ($mirror);
            }
            my $title = "[RSAT] retrieve-matrix " . &RSAT::util::AlphaDate();
            &EmailTheResult("$command", $query->param('user_email'), $result_file, "title"=>$title);
        }else{
            while(<RESULT>){
                print $_;
            }
        }
    }
    close RESULT;
}

exit(0);
