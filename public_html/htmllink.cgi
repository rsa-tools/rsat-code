#!/usr/bin/perl

## CVS
## added the possibility to specify the expected frequency for each nucleotide separately

#### this cgi script fills the HTML form for the program random-seq
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;
$title = $query->param("title");
$file = $query->param("file");

print $query->header();
print sorttable_script();
    
print $query->start_html(-title=>"$title",
    -author=>'Jacques.van-Helden\@univ-amu.fr',
    -script=>[
    { -type => 'text/javascript',
        -src      => 'js/jquery.js'
    },
    { -type => 'text/javascript',
        -src      => 'RSAT_menu.js'
    },
    { -type => 'text/javascript',
        -src      => 'RSAT_tabs.js'
    }
   ],
    -style => { 	-src => ["main.css","tabs.css","chosen.css","font-awesome.min.css"],
        -type => 'text/css',
        -media => 'screen,projection,print' });
    
    if($title =~ /RSAT/ || $title eq 'data'){
    	open($fh, "menu.php");
    	while($row = <$fh>){
            if($row =~ "<!--perlscript"){ my $x = `cat $ENV{RSAT}/public_html/data/supported_organisms.tab | wc -l`; $x -= 1; print $x; }
        	print $row."\n";
    	}
    }else{
    	$neat_java = $ENV{"neat_java_host"};
        $tomcat_port = $ENV{'tomcat_port'};
        if ($tomcat_port){
            $neat_java = $neat_java . ":" . $tomcat_port;
        }

        open(my $fh, "menu_graph.php");
        while(my $row = <$fh>){
            if(index($row, "<?php") != -1 || index($row, ";?>") != -1 ){next;}
            if(index($row, '\"') != -1) { $row =~ s/\\\"/\"/g;}
            if($row =~ /\$neat_java_host/){
                $row =~ s/\$neat_java_host/$neat_java/;
            }
            print $row."\n";
        }

    }
    
    
    print "<div class='container' id='page-content-wrapper'>";


	print "<iframe src=$file frameborder='0' width='100%' height='100%' ></iframe>";


print $query->end_html();


