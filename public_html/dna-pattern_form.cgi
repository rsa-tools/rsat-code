#!/usr/bin/perl

#### this cgi script fills the HTML form for the program dna-pattern
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

### default values for filling the form
$default{set_name} = "";
$default{patterns} = "";
$default{sequence} = "";
$default{sequence_format} = "fasta";
$default{sequence_file} = "";
$default{strands} = "both strands";

## Return formats
$default{match_positions} = 'checked';
$default{limits} = 'checked';
$default{counts} = '';
$default{table} = '';
$default{stats} = '';
$default{notacgt} = '';

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

### print the form ###
&RSA_header("dna-pattern", "form");

### head
print "<CENTER>";
print "Search a pattern (string description) within a DNA sequence<P>\n";
print "</CENTER>";

print $query->start_multipart_form(-action=>"dna-pattern.cgi");

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

### text area to enter the patterns
print "<A HREF='help.dna-pattern.html#patterns'><B>\n";
print "Query pattern(s)</B></A><BR>\n";
print $query->textarea(-name=>'patterns',
		       -default=>$default{patterns},
		       -rows=>5,
		       -columns=>60);
print "<BR>\n";

&DisplaySequenceChoice;

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


### return 
print "<BR>\n";
print CGI::table({-border=>0,-cellpadding=>3,-cellspacing=>0},
	       CGI::Tr({-align=>left,-valign=>MIDDLE},
		       [
		      CGI::td({-align=>left,-valign=>MIDDLE},
			      [
			       ## Return matching positions
			       "<A HREF='help.dna-pattern.html#return'><B>Return</B></A>\n",
			       $query->checkbox(-name=>'match_positions',
						-checked=>$default{match_positions},
						-label=>' match positions'),

			       "<A HREF='help.dna-pattern.html#origin'><B>Origin</B></A>",
			       $query->popup_menu(-name=>'origin',
						  -Values=>['start',
							    'end'],
						  -default=>$default{origin}),

			       "<A HREF='help.dna-pattern.html#flanking'><B> flanking</B></A>",
			       $query->textfield(-name=>'flanking',
						 -default=>$default{flanking},
						 -size=>2)
			      ]),
#			       "<A HREF='help.dna-pattern.html#match_format'><B>Format</B></A>",
#			       $query->popup_menu(-name=>'match_format',
#						  -Values=>['table',
#							    'fasta'],
#						  -default=>$default{match_format})
#			       ]),
			## Sequence limits
		      CGI::td({-align=>left,-valign=>MIDDLE},
			      [
			       '',
			       $query->checkbox(-name=>'limits',
						-checked=>$default{limits},
						-label=>' sequence limits'),
			       '',
			       ''
			       ]),
			    ## notacgt
				CGI::td({-align=>left,-valign=>MIDDLE},
			      [
			       '',
			       $query->checkbox(-name=>'notacgt',
						-checked=>$default{notacgt},
						-label=>' non ACGT characters'),
			       '',
			       ''
			       ]),
		      CGI::td({-align=>left,-valign=>MIDDLE},
			      [
			       '',
			       $query->checkbox(-name=>'counts',
						-checked=>$default{counts},
						-label=>' match counts'),
			       "<A HREF='help.all-upstream-search.html#threshold'><B>min count</B></A>",
			       $query->textfield(-name=>'threshold',
						 -default=>$default{threshold},
						 -size=>2)
			       ]),

		      CGI::td({-align=>left,-valign=>MIDDLE},
			      [
			       '',
			       $query->checkbox(-name=>'table',
						-checked=>$default{table},
						-label=>' match count table'),
			       $query->checkbox(-name=>'total',
						-checked=>'',
						-label=>'totals'),
			       $query->checkbox(-name=>'scores',
						-checked=>$default{scores},
						-label=>' match scores'),
			       $query->checkbox(-name=>'rank',
						-checked=>$default{rank},
						-label=>' match rank'),
			       $query->checkbox(-name=>'sort',
						-checked=>$default{sort},
						-label=>' sort'),
			       ]),

			## Statistics
		      CGI::td({-align=>left,-valign=>MIDDLE},
			      [
			       '',
			       $query->checkbox(-name=>'stats',
						-checked=>$default{stats},
						-label=>' matching statistics'),
			       '',
			       ''
			       ])
			])
		 );
print "<BR>\n";


### substitutions
print "<B><A HREF='help.dna-pattern.html#subst'>Substitutions </A></B>&nbsp;\n";
print $query->popup_menu(-name=>'subst',
			 -Values=>[0..2],
			 -default=>$default{subst});
#print $query->textfield(-name=>'subst',
#			-default=>$default{subst},
#			-size=>2);


print "<BR>\n";



### send results by email or display on the browser
&SelectOutput;

### action buttons
print "<UL><UL><TABLE class = 'formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

### data for the demo 
print $query->start_multipart_form(-action=>"dna-pattern_form.cgi");

$demo_sequence = ">PHO5   pho5 upstream sequence, from -800 to -1
TTTTACACATCGGACTGATAAGTTACTACTGCACATTGGCATTAGCTAGGAGGGCATCCA
AGTAATAATTGCGAGAAACGTGACCCAACTTTGTTGTAGGTCCGCTCCTTCTAATAATCG
CTTGTATCTCTACATATGTTCTATTTACTGACCGAAAGTAGCTCGCTACAATAATAATGT
TGACCTGATGTCAGTCCCCACGCTAATAGCGGCGTGTCGCACGCTCTCTTTACAGGACGC
CGGAGACCGGCATTACAAGGATCCGAAAGTTGTATTCAACAAGAATGCGCAAATATGTCA
ACGTATTTGGAAGTCATCTTATGTGCGCTGCTTTAATGTTTTCTCATGTAAGCGGACGTC
GTCTATAAACTTCAAACGAAGGTAAAAGGTTCATAGCGCTTTTTCTTTGTCTGCACAAAG
AAATATATATTAAATTAGCACGTTTTCGCATAGAACGCAACTGCACAATGCCAAAAAAAG
TAAAAGTGATTAAAAGAGTTAATTGAATAGGCAATCTCTAAATGAATCGATACAACCTTG
GCACTCACACGTGGGACTAGCACAGACTAAATTTATGATTCTGGTCCCTGTTTTCGAAGA
GATCGCACATGCCAAATTATCAAATTGGTCACCTTACTTGGCAAGGCATATACCCATTTG
GGATAAGGGTAAACATCTTTGAATTGTCGAAATGAAACGTATATAAGCGCTGATGTTTTG
CTAAGTCGAGGTTAGTATGGCTTCATCTCTCATGAGAATAAGAACAACAACAAATAGAGC
AAGCAAATTCGAGATTACCA
>PHO8   pho8 upstream sequence, from -800 to -1
TCTTCACCAAATTTCTTTTTTTTTTCCTACTAGAAGAAGGCGTAGCAGATAAGAAGGAAA
AATTATATTAAGCGTGCGGGTAAAGGCAAGGAAGAATCAAGTAAGACCTCAAGAATGGCA
CTATAAGTGTGGTATTATAATCTGTGTAATCCTAATTTGAGCTCTACACAATACCATTCG
ACGGTTAACAGCTACTGCATCACCGTCCAGTCATGTCGTACAACGGAATAGGGCTCAAGT
CGGCAAAAGGGTCATCTACGTCGGGCCACGTGCAGCGATCACTTGCTAGCAACAATAGGC
GCAGACCACAGGGTAGTCAACAGCAGCGGCAACAACGACAAAATGCGATCAAAAAGGCCA
GCCATGACAAGGCAAGCAGGCCTCTTGCTGTGCAGAAACAGATAGAGACTCATATGGAGA
AACGTGAGATTGAAGTACAAGTTAGCGAGCTACGGGACCGACTGGAGGAGGAAGAAACGC
TCTCGGAAGAGCAGATTGACAAGAAATGTGAAGCGTTGAGGGCAAAACTGACGAACGAGT
GGCAAGAACAGCAGCGGATGTCCTCTTTGTACACCCCTCGTAAGGCGCGTCTAACGGAAG
AGCAGCATCGACATGAATAGCAGCATTGACGATAGCGATAAGCTTCGCGCGTAGAGGAAA
AGTAAAGGGATTTTAGTATATAAAGAAAGAAGTGTATCTAAACGTTTATATTTTTTCGTG
CTCCACATTTTGCCAGCAAGTGGCTACATAAACATTTACATATCAGCATACGGGACATTA
TTTGAACGCGCATTAGCAGC
>PHO11  pho11 upstream sequence, from -800 to -1
GCAGCCTCTACCATGTTGCAAGTGCGAACCATACTGTGGCCACATAGATTACAAAAAAAG
TCCAGGATATCTTGCAAACCTAGCTTGTTTTGTAAACGACATTGAAAAAAGCGTATTAAG
GTGAAACAATCAAGATTATCTATGCCGATGAAAAATGAAAGGTATGATTTCTGCCACAAA
TATATAGTAGTTATTTTATACATCAAGATGAGAAAATAAAGGGATTTTTTCGTTCTTTTA
TCATTTTCTCTTTCTCACTTCCGACTACTTCTTATATCTACTTTCATCGTTTCATTCATC
GTGGGTGTCTAATAAAGTTTTAATGACAGAGATAACCTTGATAAGCTTTTTCTTATACGC
TGTGTCACGTATTTATTAAATTACCACGTTTTCGCATAACATTCTGTAGTTCATGTGTAC
TAAAAAAAAAAAAAAAAAAGAAATAGGAAGGAAAGAGTAAAAAGTTAATAGAAAACAGAA
CACATCCCTAAACGAAGCCGCACAATCTTGGCGTTCACACGTGGGTTTAAAAAGGCAAAT
TACACAGAATTTCAGACCCTGTTTACCGGAGAGATTCCATATTCCGCACGTCACATTGCC
AAATTGGTCATCTCACCAGATATGTTATACCCGTTTTGGAATGAGCATAAACAGCGTCGA
ATTGCCAAGTAAAACGTATATAAGCTCTTACATTTCGATAGATTCAAGCTCAGTTTCGCC
TTGGTTGTAAAGTAGGAAGAAGAAGAAGAAGAAGAGGAACAACAACAGCAAAGAGAGCAA
GAACATCATCAGAAATACCA
>PHO81  pho81 upstream sequence, from -800 to -1
AAACGAGCATGAGGGTTACAAAGAACTTCCGTTTCAAAAATGAATATAATCGTACGTTTA
CCTTGTGGCAGCACTAGCTAACGCTACGTGGAATGAACGTACCGTGCCCTATTATTCTTG
CTTGTGCTATCTCAAGAATTGCATTTTGTAATAACAACTGCATGGGAAAAATTATATAGA
TTTTCTACTATTATGTCCGCCTAAGTCAGTTAACCATCTTTATCACAAAATATACAATTA
ACCAACTACTTAATCAATTCGGTTATATTGCTTAGTATATACGTCTTTGGCACGCGATTG
AAACGCGCTAATTGCATCAGCCTATCTTTCTATGCAAGAATGCAAGAAAAATTGATGTGA
TGTGCCTTATCACAATTCATTACCTCCTATTTCCTCTGCAGCAACAAGTTTCCTTGATTA
TAAAGGTCTTTAGCGTGAGAGGTACAGGTGTTATGGCACGTGCGAATAAGGGCAGAAATT
AATCAAATTTATCAACTATTTGGCGATGGCTCGAGACAGGTATAGAACCACTACTAGGTG
ATATTGAGGCTTTTGTACAATTTATAGCAAGTTTTTGAGAGTCCCTTCAAGTTTGTTACA
TAATCTTCTTTGTGCAACGTACAAGAGCAAAGTAGAAAAATTTGGTTTTTATTTTTTTAA
GCAACATCAGCTGCACTAGTTGAGCTTTTGACAAGACATACTGCTCAAAAAATCTTCATA
ACATTATTTTTCGGTTCCACAGTGATTGAGCTTTTTGAGAGAATAACCCTTTGGAGGCAA
CATAGATAGATAAACGTGCA
>PHO84  pho84 upstream sequence, from -800 to -1
AAAAAAAAAGATTCAATAAAAAAAGAAATGAGATCAAAAAAAAAAAAAATTAAAAAAAAA
AAGAAACTAATTTATCAGCCGCTCGTTTATCAACCGTTATTACCAAATTATGAATAAAAA
AACCATATTATTATGAAAAGACACAACCGGAAGGGGAGATCACAGACCTTGACCAAGAAA
ACATGCCAAGAAATGACAGCAATCAGTATTACGCACGTTGGTGCTGTTATAGGCGCCCTA
TACGTGCAGCATTTGCTCGTAAGGGCCCTTTCAACTCATCTAGCGGCTATGAAGAAAATG
TTGCCCGGCTGAAAAACACCCGTTCCTCTCACTGCCGCACCGCCCGATGCCAATTTAATA
GTTCCACGTGGACGTGTTATTTCCAGCACGTGGGGCGGAAATTAGCGACGGCAATTGATT
ATGGTTCGCCGCAGTCCATCGAAATCAGTGAGATCGGTGCAGTTATGCACCAAATGTCGT
GTGAAAGGCTTTCCTTATCCCTCTTCTCCCGTTTTGCCTGCTTATTAGCTAGATTAAAAA
CGTGCGTATTACTCATTAATTAACCGACCTCATCTATGAGCTAATTATTATTCCTTTTTG
GCAGCATGATGCAACCACATTGCACACCGGTAATGCCAACTTAGATCCACTTACTATTGT
GGCTCGTATACGTATATATATAAGCTCATCCTCATCTCTTGTATAAAGTAAAGTTCTAAG
TTCACTTCTAAATTTTATCTTTCCTCATCTCGTAGATCACCAGGGCACACAACAAACAAA
ACTCCACGAATACAATCCAA
";
$demo_patterns = "CACGTG\nCACGTT\n";

print "<TD><B>";
print $query->hidden(-name=>'limits',-default=>"");
print $query->hidden(-name=>'patterns',-default=>$demo_patterns);
print $query->hidden(-name=>'sequence',-default=>$demo_sequence);
print $query->hidden(-name=>'organism',-default=>'Saccharomyces cerevisiae');
print $query->hidden(-name=>'set_name',-default=>'upstream sequences from the yeast PHO genes');
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


#print "<TD><B><A HREF='demo.dna-pattern.html'>DEMO</A></B></TD>\n";
print "<TD><B><A HREF='help.dna-pattern.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='tutorials/tut_dna-pattern.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);





