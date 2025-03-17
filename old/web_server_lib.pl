#!/usr/bin/env perl
use Email::Sender::Simple qw(sendmail);
use Email::Simple;
use Email::Simple::Creator;
use Email::Sender::Transport::SMTP;

################################################################
## Print a tag for Google Analytics, a tool for monitoring the access to the
## Web site.
sub google_analytics_tag {
    print <<EndGATag
<script type="text/javascript">
    var gaJsHost = (("https:" == document.location.protocol) ? "https://ssl." : "http://www.");
document.write(unescape("%3Cscript src='" + gaJsHost + "google-analytics.com/ga.js' type='text/javascript'%3E%3C/script%3E"));
</script>
<script type="text/javascript">
    try {
	var pageTracker = _gat._getTracker("xxxxxxx");
	pageTracker._trackPageview();
} catch(err) {}</script>
EndGATag
}

## ##############################################################
## Start a new HTML page and write the header of a RSAT query form or result
## page usage &RSA_header($program_name)
sub RSA_header_old {
  my $css_body_class = "form";
  my ($title) = shift;
  $title =~ s/\"//g;
  $title =~ s/\'//g;
  if (scalar @_ > 0) {
    $css_body_class = shift;
  }

  #    <link rel="stylesheet" type="text/css" href="main.css" media="screen,projection" />
#    <link rel="stylesheet" type="text/css" href="print.css" media="print" />


  #  print &html_header();
  print $query->header();

  print sorttable_script();


  ### print the header of the result page
  print $query->start_html(-title=>"RSAT : $title",
     -class => "$css_body_class",
			   -author=>'Jacques.van-Helden\@univ-amu.fr',
			   -script=>[
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
			   -style => { 	-src => ["css/main.css","css/tabs.css","css/chosen.css"],
                             	       	-type => 'text/css',
                             		-media => 'screen,projection,print' });
  print "<h3 align='center'><a href='index.php'>RSAT</a> - $title</h3>";

  ################################################################
  ## Check if client IP is blacklisted on this server.
  ##
  ## We insert this control here so it is applied to each Web form.
  &RSAT::server::DetectDeniedIP();
}

