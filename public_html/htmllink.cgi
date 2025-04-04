#!/usr/bin/env perl

## This CGI script was initially developed to fill in the HTML form
## for the random-seq web form.

# set path to lib, fixed 24012025
use FindBin '$Bin';
use lib "$Bin/lib";

## Added note, JvH 2024-09-17: this script has been used by several
## other CGI and PHP scripts, without apparent consistency between the
## formas that use it or not. This is far from clean, it should be
## revised and re-designed.
#if ($0 =~ / ([^ (\/)]+)$/) {
#    push (@INC, "$`lib/");
#}

use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";
$ENV {RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;
$title = $query->param ("title");
$file = $query->param ("file"); #$file=~tr/[ ()'":;<>-]/[         ]/; # 20181205 ajhernan / quick fix to prevent xss

print $query->header ();
print sorttable_script ();
    
print $query->start_html (-title=>"$title",
    -author=>'Jacques.van-Helden\@univ-amu.fr',
    -script=>[
    { -type => 'text/javascript',
        -src      => 'js/jquery.js'
    },
    { -type => 'text/javascript',
        -src      => 'js/RSAT_menu.js'
    },
    { -type => 'text/javascript',
        -src      => 'js/matamo.js'
    },
    { -type => 'text/javascript',
        -src      => 'js/RSAT_tabs.js'
    }
   ],
    -style => { 	-src => ["css/main.css","css/tabs.css","css/chosen.css","css/font-awesome.min.css"],
        -type => 'text/css',
        -media => 'screen,projection,print' });
    
    if ($title =~ /RSAT/ || $title eq 'data') {
        system ("php menu.php");
    } else {
    	$neat_java = $ENV {"neat_java_host"};
        $tomcat_port = $ENV {'tomcat_port'};
        if ($tomcat_port) {
            $neat_java = $neat_java . ":" . $tomcat_port;
        }

        open (my $fh, "menu_graph.php");
        while (my $row = <$fh>) {
            if (index ($row, "<?php") != -1 || index ($row, ";?>") != -1 ) {
		next;
	    }
            if (index ($row, '\"') != -1) {
		$row =~ s/\\\"/\"/g;
	    }
            if ($row =~ /\$neat_java_host/) {
                $row =~ s/\$neat_java_host/$neat_java/;
            }
            print $row."\n";
        }

    }
    
    
    print "<div class='container' id='page-content-wrapper' style='max-width:1920px;height:1080px;'>";


	print "<iframe src=$file frameborder='0' width='100%' height='100%' ></iframe>";


print $query->end_html ();


