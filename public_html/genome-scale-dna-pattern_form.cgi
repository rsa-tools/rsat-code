#!/usr/bin/perl
#### this cgi script fills the HTML form for the program genome-scale-dna-pattern
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
require "$ENV{RSAT}/public_html/genome-scale.lib.pl";

### Read the CGI query
$query = new CGI;

#### default values for dna-pattern
$default{set_name} = "";
$default{patterns} = "";
$default{strands} = "both strands";
$default{return} = "positions";
$default{noov} = "on";
$default{flanking} = "4";
$default{threshold} = "0";
$default{subst} = "0";
$default{origin} = "end";
$default{match_format} = "table";


### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
} 

### if a pattern file is specified in the query,
### read patterns from this file
if (($pattern_file = $query->param("pattern_file")) &&
     (-e $pattern_file)) {
  open PAT, $pattern_file;
  while (<PAT>) {
    $default{patterns} .= $_;
  }
  close PAT;
}

### print the form ###
&RSA_header("genome-scale dna-pattern", "form");

### head
print "<CENTER>";
print "Search a pattern (string description) within all upstream or downstream regions<P>\n";
print "</CENTER>";

#&ListParameters;


print $query->start_multipart_form(-action=>"genome-scale-dna-pattern.cgi");

################################################################
#
# retrieve-seq options
#
&DisplayRetrieveSeqOptions();


################################################################
#
# dna-pattern options
#

print $query->h4("Pattern matching options");

### text area to enter the patterns
print "<A HREF='help.dna-pattern.html#patterns'><B>\n";
print "Query pattern(s)</B></A><BR>\n";
print $query->textarea(-name=>'patterns',
		       -default=>$default{patterns},
		       -rows=>2,
		       -columns=>60);
print "<BR>\n";

### strands ###
print "<A HREF='help.dna-pattern.html#strands'><B>Search strands</B></A>&nbsp;\n";
print $query->popup_menu(-name=>'strands',
			 -Values=>['direct only',
				   'reverse complement only',
				   'both strands'],
			 -default=>$default{strands});

### prevent overlapping matches of the same pattern
print $query->checkbox(-name=>'noov',
		       -checked=>'checked',
		       -label=>'');
print "&nbsp;<A HREF='help.dna-pattern.html#noov'><B>
prevent overlapping matches
</B></A>";


### return matching positions or matching count
print "<BR>\n";
print CGI::table({-border=>0,-cellpadding=>3,-cellspacing=>0},
	       CGI::Tr({-align=>left,-valign=>MIDDLE},
		       [
		      CGI::td({-align=>left,-valign=>MIDDLE},
			      [
			       "<A HREF='help.dna-pattern.html#return'><B>Return</B></A>\n",
			       "<INPUT TYPE=RADIO NAME='return' VALUE='positions' " .
			       "CHECKED"x$default{return}=~/position/ .
			       "> match positions",
			       "<A HREF='help.all-upstream-search.html#flanking'><B> flanking residues</B></A>",
			       $query->textfield(-name=>'flanking',
						 -default=>$default{flanking},
						 -size=>2),
			       
				   "<A HREF='help.dna-pattern.html#origin'><B>Origin</B></A>",
			       $query->popup_menu(-name=>'origin',
						  -Values=>['start',
							    'end'],
						  -default=>$default{origin}),
			       "<A HREF='help.dna-pattern.html#match_format'><B>Format</B></A>",
			       $query->popup_menu(-name=>'match_format',
						  -Values=>['table',
							    'fasta'],
						  -default=>$default{match_format})

			       ]),
		      CGI::td({-align=>left,-valign=>MIDDLE},
			      [
			       '',
			       "<INPUT TYPE=RADIO NAME='return' VALUE='counts' " .
			       "CHECKED"x$default{return}=~/count/ .
			       ">match counts",
			       "<A HREF='help.all-upstream-search.html#threshold'><B> threshold on match counts</B></A>",
			       $query->textfield(-name=>'threshold',
						 -default=>$default{threshold},
						 -size=>2)
			       
			       ]),
		      CGI::td({-align=>left,-valign=>MIDDLE},
			      [
			       '',
			       "<INPUT TYPE=RADIO NAME='return' VALUE='table'>match count table",
			       $query->checkbox(-name=>'total',
						-checked=>'checked',
						-label=>'totals'),
			       ""
			       ])
			])
		 );
print "<BR>\n";


### substitutions
print "<B><A HREF='help.dna-pattern.html#subst'>Substitutions </A></B>&nbsp;\n";
print $query->popup_menu(-name=>'subst',
			 -Values=>[0..1],
			 -default=>$default{subst});
#print "<B><A HREF='help.dna-pattern.html#subst'>Substitutions </A></B>&nbsp;\n";
#print $query->textfield(-name=>'subst',
#			-default=>$default{subst},
#			-size=>2);


print "<BR>\n";



### send results by email or display on the browser
&SelectOutput;


print "<font color=red><B>Warning ! genome-scale searches can be time-consuming. If you don't obtain any result after 5 minutes,  we recommend email output.</b></font><BR>\n";


### action buttons
print "<UL><UL><TABLE class = 'formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

### data for the demo 
print $query->start_multipart_form(-action=>"genome-scale-dna-pattern_form.cgi");

$demo_patterns = "GATAAG\n";

print "<TD><B>";
print $query->hidden(-name=>'patterns',-default=>$demo_patterns);
print $query->hidden(-name=>'from',-default=>-500);
print $query->hidden(-name=>'to',-default=>-1);
print $query->hidden(-name=>'return',-default=>'match counts');
print $query->hidden(-name=>'threshold',-default=>'3');
print $query->hidden(-name=>'flanking',-default=>'0');
print $query->hidden(-name=>'organism',-default=>'Saccharomyces cerevisiae');
print $query->hidden(-name=>'set_name',-default=>'genome-scale dna-pattern');
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


print "<TD><B><A HREF='help.dna-pattern.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='tutorials/tut_genome-scale-dna-pattern.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);





