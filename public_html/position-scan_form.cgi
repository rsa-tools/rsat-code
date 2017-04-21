#!/usr/bin/perl
#### this cgi script fills the HTML form for the program position-scan
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";
require "matrix_web_forms.lib.pl";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

$default{sequence_file} = ""; ### [-f <Name of sequence file---default: standard input>]
$default{sequence} = ""; ### [-f <Name of sequence file---default: standard input>]
$default{sequence_format} = "fasta"; ### automatic conversion from any format to wc
$default{origin}="center";
$default{bg_method}="bginput";
$checked{$default{bg_method}} = "CHECKED";
$default{markov_order} = "1";
$default{organism} = "";
$default{matrix_format} = "transfac";
$default{pseudo_counts} = 1;
$default{pseudo_distribution} = "pseudo_prior";
$default{pseudo_prior} = "pseudo_prior";
$checked{$default{pseudo_prior}} = "CHECKED";
$default{bg_pseudo} = "0.01";
$default{bg_format}="oligo-analysis";
$default{decimals} = "1";
$default{class_interval} = "25";
$default{html_title} = "position-scan";

## Threshold values for site detection
$default{thresh_field} = "p_val";
$default{thresh_value} = 1e-3;

### print the form ###
&RSA_header("position-scan");
&ListParameters() if ($ENV{rsat_echo} >= 2);

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
  if ($query->param($key) =~ /checked/i) {
    $checked{$key} = "CHECKED";
  }
  if ($key eq "bg_method"){
  	$checked{$query->param($key)} = "CHECKED";
  }
}

&ReadMatrixFromFile();

### head
print "<center>";
print "Scan a set of DNA sequences (e.g. ChIP-seq peaks) with a profile matrix and display the profile of putative TFBSs along the sequences.<br>\n";
print "<br>Conception<sup>c</sup>, implementation<sup>i</sup> and testing<sup>t</sup>&nbsp: ";
print "<a target='_blank'>Jaime Castro</a><sup>cit</sup>\n";
print "<a target='_blank' href='http://www.bigre.ulb.ac.be/Users/jvanheld/'>Jacques van Helden</a><sup>cit</sup>\n";
print "</CENTER>";

# print "<div align=center>";
# print "<b>Citation</b>: <a href='mailto:jturatsi\@bigre.ulb.ac.be (Jean Valery Turatsinze)'>Jean Val&eacute;ry Turatsinze</A>, <A HREF='mailto:morgane\@bigre.ulb.ac.be (Morgane Thomas-Chollier)'>Morgane Thomas-Chollier</A>, <a href='mailto:defrance@bigre.ulb.ac.be'>Matthieu Defrance</a> and <A HREF='mailto:Jacques.van-Helden\@univ-amu.fr (Jacques van Helden)'>Jacques van Helden</a> (2008).<br> Using RSAT to scan genome sequences for transcription factor binding sites and cis-regulatory modules. Nat Protoc, 3, 1578-1588. <a href='http://www.ncbi.nlm.nih.gov/pubmed/18802439'>Pubmed 18802439</a>";
# print "</p>";
# print "</div>";

## demo description
print $default{demo_descr1};

print $query->start_multipart_form(-action=>"position-scan.cgi");


################################################################
#### Title specification
print "<hr>";
print "<h2 style='margin-left: 50px;'> Title ";
print $query->textfield(-name=>'html_title',
			 -default=>$default{html_title},
			 -size=>30) ."</h2>";


################################################################
#### Matrix specification
print "<fieldset>
<legend><b><a href='help.convert-matrix.html#io_format'>1 - Matrix </a></b></legend>";
&GetMatrix("consensus"=>0,"no_pseudo"=>1);
print "<p></p>";
print "</fieldset><p/>";


################################################################
#### Sequence specification
print "<fieldset>
<legend><b><a href='help.formats.html'>2 - Sequences </a></b></legend>";
#&MultiSequenceChoice("");
&DisplaySequenceChoice();
print "</fieldset><p/>";


################################################################
#### Background model specifiaction

print "<fieldset>
<legend><b><a href='help.matrix-scan.html#markov_order'>3 - Background </a></b></legend>";
my %bg_params =("markov" => 1,
		"bg_input" => 1,
		"markov_message" => 1,
	       );
&GetBackgroundModel(%bg_params);

#print "<br/>Note: Only Bernoulli models are supported. Higher-order Markov models are converted into Markov 0 (Bernoulli).";
print "</fieldset><p/>";


################################################################
#### scanning options

print "<fieldset>
<legend><b>4 - Scanning options</b></legend>";

## Class interval
print "<B>","&nbsp"x5,"<A HREF='help.position-analysis.html#class_interval'>window size</A>&nbsp;</B>\n";
print $query->textfield(-name=>'class_interval',
			-default=>$default{class_interval},
			-size=>3);

# ## Offset for calculating positions
# print "&nbsp;"x4,  "<A HREF='help.position-analysis.html#offset'><B>Offset</B></A>\n";
# print $query->textfield(-name=>'offset',
# 			-default=>$default{offset},
# 			-size=>8);

################################################################
#### origin for calculating positions
print "&nbsp;"x4,  "<A HREF='help.matrix-scan.html#origin'><B>Sequence Origin</B></A>\n";
print $query->popup_menu(-name=>'origin',
			 -Values=>['center',
                                  ],
			 -default=>$default{origin});
print "<br/>";

################################################################
## Fields to return + thresholds
print "&nbsp;"x4,  "<A HREF='help.matrix-scan.html#return_fields'><B>Return</B></A>\n";

my %returns = ("pval" => "sites + pval");

my $Popup = "";
    $Popup .=  "<SELECT NAME='return_field' onChange=\"toggle(this.options[this.selectedIndex].value)\">";
     foreach my $f (keys %returns) {
     	if ($f eq $default{return_field}){
			$Popup .=  "<OPTION  SELECTED VALUE=$f>$returns{$f}</option>\n";
     	} else {
     		$Popup .=  "<OPTION VALUE=$f  >$returns{$f}</option>\n";
     	}
     }
    $Popup .=  "</SELECT>";
    print $Popup;
print "<br/>";

## thresholds	 
print "&nbsp;"x4,  "<A HREF='help.matrix-scan.html#thresholds'><B>Threshold</B></A>\n";
print "<br/>";

print "<table style='padding-left:50px;'>";
print "<tr>";
print "<td rowspan=2>";
print $query->textfield(-name=>'thresh_value',
			-default=>$default{thresh_value},
			-size=>5);

print "<tr>";
print "<td style='text-align:right;'>pval <=</td>";
print "<td><i>if return is '<b>sites + pval </b>', the threshold is set on the pval</i></td>";
print "</tr>";
print "</table>";


print "</fieldset>";


################################################################
### send results by email or display on the browser
print "<BR>\n";
&SelectOutput();

################################################################
### action buttons
print "<UL><UL><TABLE>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

################################################################
### data for the demo
print $query->start_multipart_form(-action=>"position-scan_form.cgi");
my $demo_title = "Pososition profiles of vertebrates TFs on ChIP-seq peaks of JunD";
my $demo_seq_file = "$ENV{RSAT}/public_html/demo_files/Jun_Chip_seq_sequences.fasta.zip";
my $demo_seq = `gunzip -c $demo_seq_file`;

$demo_matrix = "
AC  cluster_114
XX
ID  cluster_114
XX
DE  GmGgGGGcGkGgg
P0       A     C     G     T
1        1     7    32     4
2       11    22     4     7
3        1     3    40     0
4        6     9    21     8
5        4     0    40     0
6        1     0    39     4
7        0     1    42     1
8        7    26     1    10
9        0     0    44     0
10       3     5    23    13
11       8     3    31     2
12       3     6    27     8
13       2     9    30     3
XX
CC  program: transfac
CC  matrix.nb: 18
CC  accession: cluster_114
CC  AC: cluster_114
CC  id: cluster_114
CC  name: cluster_114
CC  version: 
CC  name: cluster_114
CC  description: GmGgGGGcGkGgg
CC  transfac_consensus: 
CC  matrix.nb: 18
XX
//
AC  cluster_136
XX
ID  cluster_136
XX
DE  gAAAsyrAAw
P0       A     C     G     T
1        9     3    26     2
2       36     4     0     0
3       40     0     0     0
4       36     1     0     3
5        4    14    17     5
6        1    16     6    17
7       10     2    27     1
8       37     0     1     2
9       32     8     0     0
10      27     0     3    10
XX
CC  program: transfac
CC  matrix.nb: 42
CC  accession: cluster_136
CC  AC: cluster_136
CC  id: cluster_136
CC  name: cluster_136
CC  version: 
CC  name: cluster_136
CC  description: gAAAsyrAAw
CC  transfac_consensus: 
CC  matrix.nb: 42
XX
//
AC  cluster_137
XX
ID  cluster_137
XX
DE  GAAAgcGAAAyT
P0       A     C     G     T
1        1     0     4     0
2        5     0     0     0
3        5     0     0     0
4        5     0     0     0
5        0     1     3     1
6        1     3     0     1
7        0     0     4     1
8        5     0     0     0
9        5     0     0     0
10       5     0     0     0
11       0     3     0     2
12       0     0     1     4
XX
CC  program: transfac
CC  matrix.nb: 43
CC  accession: cluster_137
CC  AC: cluster_137
CC  id: cluster_137
CC  name: cluster_137
CC  version: 
CC  name: cluster_137
CC  description: GAAAgcGAAAyT
CC  transfac_consensus: 
CC  matrix.nb: 43
XX
//
AC  cluster_143
XX
ID  cluster_143
XX
DE  gGGggGGGggGGGtGGGrg
P0       A     C     G     T
1      701   874  4198  1179
2      490   497  4969   996
3      489   526  5136   801
4      618   681  4460  1193
5      823   813  4023  1293
6      726   602  4730   894
7      503   439  5200   810
8      587   481  4977   907
9     1019   699  4103  1131
10    1473   804  3322  1353
11    1021   289  5221   421
12      40    43  6844    25
13     133   105  6663    51
14    1167  1418   340  4027
15      81    59  6744    68
16      70    72  6107   703
17     168   295  5999   490
18    1758  1021  2686  1487
19    1646   893  3135  1278
XX
CC  program: transfac
CC  matrix.nb: 50
CC  accession: cluster_143
CC  AC: cluster_143
CC  id: cluster_143
CC  name: cluster_143
CC  version: 
CC  name: cluster_143
CC  description: gGGggGGGggGGGtGGGrg
CC  transfac_consensus: 
CC  matrix.nb: 50
XX
//
AC  cluster_144
XX
ID  cluster_144
XX
DE  gggrkGyGkGGGhGGgrg
P0       A     C     G     T
1       95    96   219    90
2      104   103   223    70
3      106    72   226    96
4      135    70   234    61
5      111    37   190   162
6       58    41   378    23
7       84   242    25   149
8       21    21   404    54
9       17     4   226   253
10      80     3   402    15
11       4     0   488     8
12      16    19   457     8
13     131   226     3   140
14      90     6   375    29
15      36     6   387    71
16     118    46   267    69
17     139    49   229    83
18     109    92   245    54
XX
CC  program: transfac
CC  matrix.nb: 51
CC  accession: cluster_144
CC  AC: cluster_144
CC  id: cluster_144
CC  name: cluster_144
CC  version: 
CC  name: cluster_144
CC  description: gggrkGyGkGGGhGGgrg
CC  transfac_consensus: 
CC  matrix.nb: 51
XX
//
AC  cluster_149
XX
ID  cluster_149
XX
DE  wGcTGAcks
P0       A     C     G     T
1       15     2     2    28
2        1     0    45     1
3       11    31     4     1
4        1     0     0    46
5        0     4    39     4
6       41     3     0     3
7        3    31     7     6
8       11     5    16    15
9        7    13    19     8
XX
CC  program: transfac
CC  matrix.nb: 56
CC  accession: cluster_149
CC  AC: cluster_149
CC  id: cluster_149
CC  name: cluster_149
CC  version: 
CC  name: cluster_149
CC  description: wGcTGAcks
CC  transfac_consensus: 
CC  matrix.nb: 56
XX
//
AC  cluster_15
XX
ID  cluster_15
XX
DE  mATGACt
P0       A     C     G     T
1       27    29  18.5    15
2     68.5     5    15     1
3      1.5   0.5     1  86.5
4      4.5     0  83.5   1.5
5       89   0.5     0     0
6      1.5  71.5   2.5    14
7       19     7  21.5    42
XX
CC  program: transfac
CC  matrix.nb: 57
CC  accession: cluster_15
CC  AC: cluster_15
CC  id: cluster_15
CC  name: cluster_15
CC  version: 
CC  name: cluster_15
CC  description: mATGACt
CC  transfac_consensus: 
CC  matrix.nb: 57
XX
//
AC  cluster_150
XX
ID  cluster_150
XX
DE  tGCTGAckyar
P0       A     C     G     T
1       81   130    52   542
2       18    13   766     8
3       64   706     6    29
4       24    30     3   748
5       16     4   689    96
6      791     2     3     9
7       14   535   110   146
8      126   105   344   230
9      131   299   100   275
10     347   136   140   182
11     214    96   333   162
XX
CC  program: transfac
CC  matrix.nb: 58
CC  accession: cluster_150
CC  AC: cluster_150
CC  id: cluster_150
CC  name: cluster_150
CC  version: 
CC  name: cluster_150
CC  description: tGCTGAckyar
CC  transfac_consensus: 
CC  matrix.nb: 58
XX
//
AC  cluster_160
XX
ID  cluster_160
XX
DE  ggAAAsyGAAAsbrrra
P0       A     C     G     T
1        6     6    29     3
2       10     2    30     2
3       39     1     4     0
4       43     0     0     1
5       40     0     2     2
6        7    11    25     1
7        6    16     9    13
8        4     1    38     1
9       43     0     1     0
10      43     0     0     1
11      42     2     0     0
12       9    21    11     3
13       1    15    12    16
14      13     3    20     8
15      19     7    15     3
16      22     4    13     5
17      22     7     7     8
XX
CC  program: transfac
CC  matrix.nb: 69
CC  accession: cluster_160
CC  AC: cluster_160
CC  id: cluster_160
CC  name: cluster_160
CC  version: 
CC  name: cluster_160
CC  description: ggAAAsyGAAAsbrrra
CC  transfac_consensus: 
CC  matrix.nb: 69
XX
//
AC  cluster_163
XX
ID  cluster_163
XX
DE  GGGmGGGGssggggggGgGgGg
P0       A     C     G     T
1       31    30   430     9
2        6    14   475     5
3        4    11   484     1
4      209   261     4    26
5        3     9   487     1
6        6    10   481     3
7       66    35   380    19
8       25    23   448     4
9       55   194   230    21
10      66   216   187    31
11      73    58   324    45
12      89   101   298    12
13      91    62   323    24
14      70    81   335    14
15      70    72   333    25
16      94   108   270    28
17      62    68   351    19
18      84    98   306    12
19      39    88   358    15
20      76   112   292    20
21      37    77   366    20
22      30   117   315    38
XX
CC  program: transfac
CC  matrix.nb: 72
CC  accession: cluster_163
CC  AC: cluster_163
CC  id: cluster_163
CC  name: cluster_163
CC  version: 
CC  name: cluster_163
CC  description: GGGmGGGGssggggggGgGgGg
CC  transfac_consensus: 
CC  matrix.nb: 72
XX
//
AC  cluster_164
XX
ID  cluster_164
XX
DE  gGGGGCGGGGCcgGgsssgggs
P0       A     C     G     T
1       89    58   336    12
2       52    25   365    53
3       61    41   389     4
4       12    25   448    10
5       12     8   472     3
6       32   436     8    19
7        4    18   461    12
8       20    31   403    41
9       64    23   379    29
10      34    33   414    14
11      43   363    66    23
12      39   326    83    47
13     100    71   205   119
14      26    87   355    27
15      60   112   302    21
16      75   124   259    37
17      41   169   243    42
18      40   155   263    37
19      82   104   275    34
20      85   106   256    48
21      45   100   321    29
22      58   162   235    40
XX
CC  program: transfac
CC  matrix.nb: 73
CC  accession: cluster_164
CC  AC: cluster_164
CC  id: cluster_164
CC  name: cluster_164
CC  version: 
CC  name: cluster_164
CC  description: gGGGGCGGGGCcgGgsssgggs
CC  transfac_consensus: 
CC  matrix.nb: 73
XX
//
AC  cluster_179
XX
ID  cluster_179
XX
DE  rrmtTGAytGAtk
P0       A     C     G     T
1        4     2     3     0
2        3     1     4     1
3        3     3     2     1
4        0     2     2     5
5        0     0     0     9
6        0     0     9     0
7        9     0     0     0
8        0     4     0     5
9        2     1     1     5
10       0     0     8     1
11       7     2     0     0
12       1     2     1     5
13       2     0     4     3
XX
CC  program: transfac
CC  matrix.nb: 89
CC  accession: cluster_179
CC  AC: cluster_179
CC  id: cluster_179
CC  name: cluster_179
CC  version: 
CC  name: cluster_179
CC  description: rrmtTGAytGAtk
CC  transfac_consensus: 
CC  matrix.nb: 89
XX
//
AC  cluster_180
XX
ID  cluster_180
XX
DE  tGCTGAGTCA
P0       A     C     G     T
1        2     3     0    11
2        0     0    16     0
3        0    16     0     0
4        0     3     0    13
5        1     0    14     1
6       15     0     0     1
7        0     3    13     0
8        0     0     0    16
9        0    16     0     0
10      16     0     0     0
XX
CC  program: transfac
CC  matrix.nb: 91
CC  accession: cluster_180
CC  AC: cluster_180
CC  id: cluster_180
CC  name: cluster_180
CC  version: 
CC  name: cluster_180
CC  description: tGCTGAGTCA
CC  transfac_consensus: 
CC  matrix.nb: 91
XX
//
AC  cluster_188
XX
ID  cluster_188
XX
DE  kgggkGGGGGaGGGG
P0       A     C     G     T
1        5     4    16     9
2        3     7    20     4
3        8     7    15     4
4        4     3    20     7
5        6     6     9    13
6        1     3    29     1
7        3     0    31     0
8        0     0    34     0
9        2     0    30     2
10       0     0    34     0
11      16     7     5     6
12       7     0    26     1
13       1     1    32     0
14       7     0    26     1
15       4     3    25     2
XX
CC  program: transfac
CC  matrix.nb: 99
CC  accession: cluster_188
CC  AC: cluster_188
CC  id: cluster_188
CC  name: cluster_188
CC  version: 
CC  name: cluster_188
CC  description: kgggkGGGGGaGGGG
CC  transfac_consensus: 
CC  matrix.nb: 99
XX
//
AC  cluster_189
XX
ID  cluster_189
XX
DE  kkkGGGGGGGyrgk
P0       A     C     G     T
1      541   315   912   643
2      331   186  1002   892
3      190   115  1228   878
4       44    27  2234   106
5       26     8  2320    57
6       29     5  2353    24
7       28    11  2325    47
8       21     4  2328    58
9       31     4  2016   360
10       1     1  2409     0
11     562   754    42  1053
12     777   207   991   436
13     518   305  1025   563
14     439   275  1038   659
XX
CC  program: transfac
CC  matrix.nb: 100
CC  accession: cluster_189
CC  AC: cluster_189
CC  id: cluster_189
CC  name: cluster_189
CC  version: 
CC  name: cluster_189
CC  description: kkkGGGGGGGyrgk
CC  transfac_consensus: 
CC  matrix.nb: 100
XX
//
AC  cluster_206
XX
ID  cluster_206
XX
DE  raAAGtGAAAGTga
P0       A     C     G     T
1      183    37   225    53
2      316    48    72    62
3      359    58    50    31
4      440     1    50     7
5       34     3   450    11
6       84    12    81   321
7       16     2   479     1
8      494     1     1     2
9      477     0    20     1
10     486     1     4     7
11      11     6   469    12
12       8    38    24   428
13     114    37   240   107
14     216   100    87    95
XX
CC  program: transfac
CC  matrix.nb: 120
CC  accession: cluster_206
CC  AC: cluster_206
CC  id: cluster_206
CC  name: cluster_206
CC  version: 
CC  name: cluster_206
CC  description: raAAGtGAAAGTga
CC  transfac_consensus: 
CC  matrix.nb: 120
XX
//
AC  cluster_227
XX
ID  cluster_227
XX
DE  GGGtGGGyAGGGGtGGGcdk
P0       A     C     G     T
1      541   368  3364   443
2      399   376  3451   490
3      317   375  3491   533
4      784  1094   865  1973
5      306   238  3926   246
6      237   186  3815   478
7      255   129  4104   228
8      908  1931   695  1182
9     3778   187   520   231
10      84    72  4473    87
11      18    28  4596    74
12      16    66  3554  1080
13      32    41  4599    44
14     351  1037   234  3094
15      25    40  4514   137
16    1005    20  3657    34
17     833    79  3642   162
18    1127  1894   821   874
19    1244   493  1610  1369
20     990   603  1766  1357
XX
CC  program: transfac
CC  matrix.nb: 143
CC  accession: cluster_227
CC  AC: cluster_227
CC  id: cluster_227
CC  name: cluster_227
CC  version: 
CC  name: cluster_227
CC  description: GGGtGGGyAGGGGtGGGcdk
CC  transfac_consensus: 
CC  matrix.nb: 143
XX
//
AC  cluster_230
XX
ID  cluster_230
XX
DE  dtyGCGTGm
P0       A     C     G     T
1       41    18    56    39
2       11    12    35    96
3       22    44    21    67
4        3     1   146     4
5        1   150     1     2
6        3     1   149     1
7        0     3     1   150
8        0     0   154     0
9       43    67    16    28
XX
CC  program: transfac
CC  matrix.nb: 147
CC  accession: cluster_230
CC  AC: cluster_230
CC  id: cluster_230
CC  name: cluster_230
CC  version: 
CC  name: cluster_230
CC  description: dtyGCGTGm
CC  transfac_consensus: 
CC  matrix.nb: 147
XX
//
AC  cluster_245
XX
ID  cluster_245
XX
DE  gkGGGyGkGksgGGgGGGg
P0       A     C     G     T
1      126   106   291   127
2       64    79   343   164
3       40    22   558    30
4        6    14   617    13
5       22     6   598    24
6       83   315     8   244
7       30     9   579    32
8       26    20   366   238
9       43     4   567    36
10     101    32   341   176
11      91   174   233   152
12     123    57   385    85
13      48    21   526    55
14      40    42   512    56
15      50    45   440   115
16      47    71   453    79
17      38    82   462    68
18      27    68   464    91
19      67    95   340   148
XX
CC  program: transfac
CC  matrix.nb: 163
CC  accession: cluster_245
CC  AC: cluster_245
CC  id: cluster_245
CC  name: cluster_245
CC  version: 
CC  name: cluster_245
CC  description: gkGGGyGkGksgGGgGGGg
CC  transfac_consensus: 
CC  matrix.nb: 163
XX
//
AC  cluster_266
XX
ID  cluster_266
XX
DE  raaagtGAAAGtGAAAsyra
P0       A     C     G     T
1       86    49   113    46
2      167    32    67    28
3      192    19    48    35
4      191    30    51    22
5       58    68   128    40
6       66    66    66    96
7       61     4   223     6
8      273     3    14     4
9      276     7     9     2
10     282     4     4     4
11      35    50   207     2
12      31    64     7   192
13       9     4   281     0
14     252     7    32     3
15     269     9    14     2
16     280     2     9     3
17       5   103   180     6
18      19    74    16   185
19     111    39   117    27
20     139    56    73    26
XX
CC  program: transfac
CC  matrix.nb: 186
CC  accession: cluster_266
CC  AC: cluster_266
CC  id: cluster_266
CC  name: cluster_266
CC  version: 
CC  name: cluster_266
CC  description: raaagtGAAAGtGAAAsyra
CC  transfac_consensus: 
CC  matrix.nb: 186
XX
//
AC  cluster_285
XX
ID  cluster_285
XX
DE  gsssssgsgsvrggGGCGGGrv
P0       A     C     G     T
1      100   122   242    36
2       95   142   214    49
3       60   130   259    51
4      117   130   210    43
5       89   152   227    32
6       86   174   199    41
7       78    94   245    83
8       81   150   228    41
9       87   113   242    58
10      96   136   224    44
11     131   148   188    33
12     136    70   249    45
13     101    56   282    61
14     112    72   269    47
15      11    50   433     6
16       2     7   486     5
17      29   453     1    17
18       3     7   476    14
19      10     6   475     9
20      64    21   405    10
21     298    20   168    14
22     197   163   131     9
XX
CC  program: transfac
CC  matrix.nb: 207
CC  accession: cluster_285
CC  AC: cluster_285
CC  id: cluster_285
CC  name: cluster_285
CC  version: 
CC  name: cluster_285
CC  description: gsssssgsgsvrggGGCGGGrv
CC  transfac_consensus: 
CC  matrix.nb: 207
XX
//
AC  cluster_30
XX
ID  cluster_30
XX
DE  ktyTGTGGTTwk
P0       A     C     G     T
1    28.3333    35 50.6667 52.6667
2    44.3333    83 83.3333 130.333
3       28 147.667    22 143.333
4    9.33333 4.66667 3.33333 323.667
5    0.333333 2.33333   336 2.33333
6    7.33333 59.6667     2   272
7    1.33333     0 339.333 0.333333
8    0.333333     0 338.333 2.33333
9    0.333333 19.6667 11.3333 309.667
10      10    73 21.3333 236.667
11   116.333     8    37 179.667
12   24.6667 42.6667    63    44
XX
CC  program: transfac
CC  matrix.nb: 224
CC  accession: cluster_30
CC  AC: cluster_30
CC  id: cluster_30
CC  name: cluster_30
CC  version: 
CC  name: cluster_30
CC  description: ktyTGTGGTTwk
CC  transfac_consensus: 
CC  matrix.nb: 224
XX
//
AC  cluster_312
XX
ID  cluster_312
XX
DE  kttgwkTGwTTwy
P0       A     C     G     T
1        2     3     4     4
2        0     2     2     9
3        3     0     1     9
4        2     1     9     1
5        6     0     0     7
6        0     2     4     7
7        0     1     0    12
8        2     1    10     0
9        4     1     1     7
10       0     0     1    12
11       0     1     2    10
12       4     2     1     6
13       1     6     1     5
XX
CC  program: transfac
CC  matrix.nb: 238
CC  accession: cluster_312
CC  AC: cluster_312
CC  id: cluster_312
CC  name: cluster_312
CC  version: 
CC  name: cluster_312
CC  description: kttgwkTGwTTwy
CC  transfac_consensus: 
CC  matrix.nb: 238
XX
//
AC  cluster_349
XX
ID  cluster_349
XX
DE  tTTTTATtytTTTTttt
P0       A     C     G     T
1      295   393   161   954
2        1     1     1  1800
3        9     7     1  1786
4        1     3     1  1798
5        0     3     0  1800
6     1797     1     4     1
7        4    98     6  1695
8      231   425    51  1096
9      328   452   132   891
10     167   403    87  1146
11     146   332    62  1263
12     118   245    47  1393
13     181   213    29  1380
14     327   179    21  1276
15     404   210    42  1147
16     433   304    63  1003
17     403   379   127   894
XX
CC  program: transfac
CC  matrix.nb: 278
CC  accession: cluster_349
CC  AC: cluster_349
CC  id: cluster_349
CC  name: cluster_349
CC  version: 
CC  name: cluster_349
CC  description: tTTTTATtytTTTTttt
CC  transfac_consensus: 
CC  matrix.nb: 278
XX
//
AC  cluster_377
XX
ID  cluster_377
XX
DE  kggkrrGRcGkskkggg
P0       A     C     G     T
1      116   125   333   199
2      125   118   353   177
3      102   113   383   175
4       61    89   291   332
5      416    24   200   133
6      380     8   370    15
7       21     4   734    14
8      203     4   542    24
9      140   429    41   163
10      28     2   696    47
11      19    63   360   331
12      13   343   386    31
13     119    35   213   406
14      64    57   229   423
15     140    97   467    69
16     120    64   477   112
17     166   103   357   147
XX
CC  program: transfac
CC  matrix.nb: 309
CC  accession: cluster_377
CC  AC: cluster_377
CC  id: cluster_377
CC  name: cluster_377
CC  version: 
CC  name: cluster_377
CC  description: kggkrrGRcGkskkggg
CC  transfac_consensus: 
CC  matrix.nb: 309
XX
//
AC  cluster_408
XX
ID  cluster_408
XX
DE  CaAkrGyGrTkGyGr
P0       A     C     G     T
1        1    25     6     0
2       19     5     6     2
3       28     2     2     0
4        2     2     8    20
5       20     0    12     0
6        6     0    26     0
7        0    21     0    11
8        0     0    32     0
9       12     1    12     7
10       0     1     7    24
11       6     4    12    10
12       1     0    26     5
13       0    11     6    15
14       3     1    28     0
15       9     5    13     5
XX
CC  program: transfac
CC  matrix.nb: 344
CC  accession: cluster_408
CC  AC: cluster_408
CC  id: cluster_408
CC  name: cluster_408
CC  version: 
CC  name: cluster_408
CC  description: CaAkrGyGrTkGyGr
CC  transfac_consensus: 
CC  matrix.nb: 344
XX
//
AC  cluster_411
XX
ID  cluster_411
XX
DE  sGccGCCatstyggctscGGGC
P0       A     C     G     T
1       69   222   161    46
2       31    70   369    28
3       31   336    93    38
4       75   291    89    43
5       33    40   366    59
6        8   456    27     7
7       18   452    21     7
8      298    91    70    39
9        9   113    88   288
10      30   261   150    57
11      20    91   104   283
12      12   180    25   281
13      51    93   251   103
14      81   118   264    35
15      92   196   114    96
16      56   123    88   231
17      59   170   170    99
18      87   249    56   106
19      12    64   341    81
20      12    74   394    18
21       9    97   371    21
22      36   391    34    37
XX
CC  program: transfac
CC  matrix.nb: 348
CC  accession: cluster_411
CC  AC: cluster_411
CC  id: cluster_411
CC  name: cluster_411
CC  version: 
CC  name: cluster_411
CC  description: sGccGCCatstyggctscGGGC
CC  transfac_consensus: 
CC  matrix.nb: 348
XX
//
AC  cluster_414
XX
ID  cluster_414
XX
DE  svsssscscAGGCCbsGsc
P0       A     C     G     T
1       65   177   179    59
2      145   140   146    49
3       48   121   270    41
4       77   150   214    39
5       53   134   246    47
6       57   179   163    81
7       71   225   113    71
8       43   189   163    85
9       35   279    96    70
10     476     2     1     1
11       2     2   476     0
12       0     1   479     0
13       0   477     3     0
14       1   425    42    12
15       9   157   128   186
16     108   197   144    31
17      97     4   373     6
18      67   165   238    10
19      61   252    59   108
XX
CC  program: transfac
CC  matrix.nb: 351
CC  accession: cluster_414
CC  AC: cluster_414
CC  id: cluster_414
CC  name: cluster_414
CC  version: 
CC  name: cluster_414
CC  description: svsssscscAGGCCbsGsc
CC  transfac_consensus: 
CC  matrix.nb: 351
XX
//
AC  cluster_419
XX
ID  cluster_419
XX
DE  rkskkkGkkGkTGkTTkGkGkt
P0       A     C     G     T
1        7     3    10     6
2        2     4    11     9
3        0     7    16     3
4        2     0    15     9
5        2     3    12     9
6        4     1    11    10
7        0     1    20     5
8        0     4    13     9
9        0     1    17     8
10       0     1    23     2
11       0     3    15     8
12       0     0     2    24
13       1     0    24     1
14       2     3    12     9
15       0     1     6    19
16       2     1     4    19
17       0     1     8    17
18       0     2    24     0
19       2     0    17     7
20       0     1    22     3
21       0     0    17     9
22       4     4     5    13
XX
CC  program: transfac
CC  matrix.nb: 356
CC  accession: cluster_419
CC  AC: cluster_419
CC  id: cluster_419
CC  name: cluster_419
CC  version: 
CC  name: cluster_419
CC  description: rkskkkGkkGkTGkTTkGkGkt
CC  transfac_consensus: 
CC  matrix.nb: 356
XX
//
AC  cluster_427
XX
ID  cluster_427
XX
DE  asTyyyrwwwTGAsTcAk
P0       A     C     G     T
1      246    48   119    86
2       58   179   148   114
3       21    11    67   400
4       70   172    28   229
5       48   164    31   256
6       88   237    43   131
7      139    55   209    96
8      178    55    95   171
9      147     9    51   292
10     309    49    12   129
11       0     0     0   499
12       1     8   448    42
13     482     0     5    12
14      38   260   190    11
15      17     7     2   473
16     106   308    66    19
17     363    10    48    78
18      81    89   198   131
XX
CC  program: transfac
CC  matrix.nb: 365
CC  accession: cluster_427
CC  AC: cluster_427
CC  id: cluster_427
CC  name: cluster_427
CC  version: 
CC  name: cluster_427
CC  description: asTyyyrwwwTGAsTcAk
CC  transfac_consensus: 
CC  matrix.nb: 365
XX
//
AC  cluster_5
XX
ID  cluster_5
XX
DE  grTGAGTCAycs
P0       A     C     G     T
1      127 103.556 214.778 107.222
2    540.556 279.111   318 93.1111
3    25.2222 7.66667 23.8889  1174
4    9.88889    24 1133.56 63.3333
5    1155.33 31.4444 10.7778 33.2222
6    87.4444   159 959.333    25
7    19.2222 11.1111 14.4444  1186
8       49 1168.67 9.44444 3.66667
9    1206.44 6.11111 5.88889 12.3333
10      58 377.222 240.111 555.444
11   241.222 350.111   148 244.111
12   41.4444 86.1111 57.7778 36.8889
XX
CC  program: transfac
CC  matrix.nb: 391
CC  accession: cluster_5
CC  AC: cluster_5
CC  id: cluster_5
CC  name: cluster_5
CC  version: 
CC  name: cluster_5
CC  description: grTGAGTCAycs
CC  transfac_consensus: 
CC  matrix.nb: 391
XX
//
AC  cluster_58
XX
ID  cluster_58
XX
DE  gTGAykyma
P0       A     C     G     T
1    239.5 226.5   457   140
2      7.5   4.5     6  1045
3        2    10  1048     3
4     1031   7.5   4.5    20
5        7   694  72.5 289.5
6       27    66   602   368
7       57   516  14.5 475.5
8      354 469.5   133 106.5
9    605.5   132 164.5   161
XX
CC  program: transfac
CC  matrix.nb: 400
CC  accession: cluster_58
CC  AC: cluster_58
CC  id: cluster_58
CC  name: cluster_58
CC  version: 
CC  name: cluster_58
CC  description: gTGAykyma
CC  transfac_consensus: 
CC  matrix.nb: 400
XX
//
AC  cluster_62
XX
ID  cluster_62
XX
DE  TGCTGAGTCAyssy
P0       A     C     G     T
1    81.6667    78 43.6667 471.667
2    49.3333    25   589 11.6667
3    21.6667 620.333 5.33333 27.6667
4    41.6667 54.3333    42   537
5       28 43.6667 555.667 47.6667
6    535.667    25 41.6667 72.6667
7       17 109.333 517.333 31.3333
8       21 15.3333    29 609.667
9    7.33333 623.667    39     5
10     592 45.3333 16.6667    21
11      25 213.333 63.6667   373
12   86.3333   183   247 158.667
13   67.6667 258.667   204   114
14   60.3333 132.667 93.6667   191
XX
CC  program: transfac
CC  matrix.nb: 405
CC  accession: cluster_62
CC  AC: cluster_62
CC  id: cluster_62
CC  name: cluster_62
CC  version: 
CC  name: cluster_62
CC  description: TGCTGAGTCAyssy
CC  transfac_consensus: 
CC  matrix.nb: 405
XX
//
AC  cluster_70
XX
ID  cluster_70
XX
DE  rraaagrGGAAsTGArAv
P0       A     C     G     T
1       61    26    39    25
2    140.333 44.3333 97.6667 35.3333
3    280.667 35.3333 105.333    63
4    274.667    41   112 56.6667
5      279    16 101.333    88
6    91.6667 48.6667 318.333 25.6667
7    287.333 25.6667 149.667 21.6667
8       65 3.66667 400.667    15
9    57.3333    13 411.333 2.66667
10   451.333 11.6667 14.3333     7
11     428 2.33333 33.6667 20.3333
12   65.6667 131.333 279.667 7.66667
13      60 38.3333 28.6667 357.333
14   24.3333 19.6667   427 13.3333
15     336    27 107.667 13.6667
16   224.333 54.6667   171 34.3333
17   330.333    45 53.3333 55.6667
18   88.3333 92.3333 120.667    32
XX
CC  program: transfac
CC  matrix.nb: 414
CC  accession: cluster_70
CC  AC: cluster_70
CC  id: cluster_70
CC  name: cluster_70
CC  version: 
CC  name: cluster_70
CC  description: rraaagrGGAAsTGArAv
CC  transfac_consensus: 
CC  matrix.nb: 414
XX
//
AC  cluster_8
XX
ID  cluster_8
XX
DE  aaawwTGCTGAsTcAGCAwwww
P0       A     C     G     T
1    209.667    46    88 66.3333
2      210 33.6667 73.6667 92.6667
3    257.333    23 32.6667    97
4    201.667    21 44.6667 142.667
5      127    81 81.3333 120.667
6       26 11.6667     9 363.333
7        1     2 404.667 2.33333
8    40.6667   364 1.66667 3.66667
9    10.3333 2.66667 0.666667 396.333
10   6.33333 0.666667 392.333 10.6667
11   397.333     2     5 5.66667
12   19.3333   203 156.333 31.3333
13   40.6667    27 16.6667 325.667
14      51 273.667    22 63.3333
15   301.667 3.33333 56.3333 48.6667
16   27.6667     7   324 51.3333
17   9.66667 340.667 47.6667    12
18     298    33 33.6667 45.3333
19   174.667    54    77 104.333
20   193.333 26.6667 40.3333 149.667
21     165 35.6667    37 172.333
22      80    57    23 128.667
XX
CC  program: transfac
CC  matrix.nb: 424
CC  accession: cluster_8
CC  AC: cluster_8
CC  id: cluster_8
CC  name: cluster_8
CC  version: 
CC  name: cluster_8
CC  description: aaawwTGCTGAsTcAGCAwwww
CC  transfac_consensus: 
CC  matrix.nb: 424
XX
//
AC  cluster_81
XX
ID  cluster_81
XX
DE  ggggrgGGGGcGGGGcsrgsgs
P0       A     C     G     T
1    198.667 226.333 614.333 104.667
2    273.667 242.667   605   161
3      210 214.667 730.333 127.333
4    311.333   279 565.333 126.667
5    436.333 241.333   429 175.667
6    261.333 167.667 755.667 97.6667
7    145.333    41 919.333 176.667
8    227.333    14 1031.33 9.66667
9    9.66667 6.66667 1254.33 11.6667
10   9.66667 17.3333  1252 3.33333
11   218.333 868.333 29.6667   166
12   28.3333    17  1192    45
13   17.6667    20 1184.33 60.3333
14   203.667 68.3333   958 52.3333
15     133   117 978.667 53.6667
16     177   714 256.333   135
17     117   558 338.667 268.667
18   332.333 203.667 490.667 255.667
19   251.667 282.667   623   125
20      53 79.6667 140.333 31.6667
21   50.3333 64.3333 160.333 29.6667
22   20.3333 54.3333 82.3333 9.33333
XX
CC  program: transfac
CC  matrix.nb: 426
CC  accession: cluster_81
CC  AC: cluster_81
CC  id: cluster_81
CC  name: cluster_81
CC  version: 
CC  name: cluster_81
CC  description: ggggrgGGGGcGGGGcsrgsgs
CC  transfac_consensus: 
CC  matrix.nb: 426
XX
//
AC  cluster_84
XX
ID  cluster_84
XX
DE  ggrrAAmyGAAAstga
P0       A     C     G     T
1        2     2 6.33333     2
2    169.333 95.6667   297 153.333
3    194.667 46.6667   458    16
4    448.667 10.3333   219 37.3333
5      689 6.66667 13.6667     6
6    665.667    13 21.6667    15
7    262.667 195.333 170.333    87
8    152.333   194   169   200
9    97.6667 12.6667   595    10
10     617 12.3333 79.6667 6.33333
11     644 44.6667    20 6.66667
12   626.333 7.66667    27 54.3333
13   16.3333 401.333 256.333 41.3333
14   92.3333 127.333    47 448.667
15     170   101 314.667 129.667
16     380   148 66.3333    81
XX
CC  program: transfac
CC  matrix.nb: 429
CC  accession: cluster_84
CC  AC: cluster_84
CC  id: cluster_84
CC  name: cluster_84
CC  version: 
CC  name: cluster_84
CC  description: ggrrAAmyGAAAstga
CC  transfac_consensus: 
CC  matrix.nb: 429
XX
//
AC  cluster_9
XX
ID  cluster_9
XX
DE  ATGACGTCAycr
P0       A     C     G     T
1    2431.33 162.5   332 54.3333
2    36.6667 36.1667    22 2972.67
3    24.1667  24.5  2382 636.833
4     2985  14.5 45.1667 22.8333
5    24.8333  3021 3.33333 18.3333
6    49.1667 21.8333 2987.83 8.66667
7        1 12.3333     1 3053.17
8    43.3333  3020     1 3.16667
9    3048.5 11.1667     3 4.83333
10    32.5  1137 40.3333 1857.67
11     659  1457 300.833 650.667
12     222 133.667 155.5    82
XX
CC  program: transfac
CC  matrix.nb: 435
CC  accession: cluster_9
CC  AC: cluster_9
CC  id: cluster_9
CC  name: cluster_9
CC  version: 
CC  name: cluster_9
CC  description: ATGACGTCAycr
CC  transfac_consensus: 
CC  matrix.nb: 435
XX
//
AC  cluster_97
XX
ID  cluster_97
XX
DE  gGGGGcGGGGc
P0       A     C     G     T
1      683   443  2076   287
2      300   104  2614   471
3      591    43  2845    10
4        7    45  3433     4
5       91    19  3374     5
6      809  2247    62   371
7      121    42  3242    84
8       14    36  3332   107
9      558    90  2695   146
10     264   255  2797   173
11     519  1875   655   440
XX
CC  program: transfac
CC  matrix.nb: 443
CC  accession: cluster_97
CC  AC: cluster_97
CC  id: cluster_97
CC  name: cluster_97
CC  version: 
CC  name: cluster_97
CC  description: gGGGGcGGGGc
CC  transfac_consensus: 
CC  matrix.nb: 443
XX
//
AC  cluster_98
XX
ID  cluster_98
XX
DE  GGsGGyrGGGm
P0       A     C     G     T
1        0     0     6     0
2        0     0     6     0
3        0     3     2     1
4        0     0     6     0
5        0     0     6     0
6        0     3     1     2
7        3     0     3     0
8        0     0     5     1
9        0     1     5     0
10       0     0     6     0
11       2     2     1     1
XX
CC  program: transfac
CC  matrix.nb: 444
CC  accession: cluster_98
CC  AC: cluster_98
CC  id: cluster_98
CC  name: cluster_98
CC  version: 
CC  name: cluster_98
CC  description: GGsGGyrGGGm
CC  transfac_consensus: 
CC  matrix.nb: 444
XX
//
";

$descr="<H4>Comment on the demonstration example : </H4><blockquote class ='demo'>

In this demonstration, we will scan a set of PSSMs on ChIP-seq peaks (length = 310 nt) of JunD on GM12878 cells, taken from the Encode project website. The program shows the binding profiles of each input PSSMs, some of them displaying a hill profile, others a valley profile and other a flat distribution.

Data set: wgEncodeAwgTfbsSydhGm12878JundUniPk.narrowPeak

The Null Hypothesis is that the distribution of TFBSs follows a flat distribution.<p/>";

## demo 1
print "<TD><B>";
print $query->hidden(-name=>'html_title',-default=>$demo_title);
print $query->hidden(-name=>'bg_method',-default=>'bginput');
print $query->hidden(-name=>'thresh_field',-default=>'pval');
print $query->hidden(-name=>'thresh_value',-default=>'1e-3');
print $query->hidden(-name=>'bgfile',-default=>'CHECKED');
print $query->hidden(-name=>'background',-default=>'');
print $query->hidden(-name=>'markov_order',-default=>'1');
print $query->hidden(-name=>'organism',-default=>'');
print $query->hidden(-name=>'return_field',-default=>'pval');
print $query->hidden(-name=>'matrix',-default=>$demo_matrix);
print $query->hidden(-name=>'matrix_format',-default=>'transfac');
print $query->hidden(-name=>'sequence',-default=>$demo_seq);
print $query->hidden(-name=>'sequence_format',-default=>$default{sequence_format});
print $query->hidden(-name=>'origin',-default=>'center');
print $query->submit(-label=>"DEMO");
print "</B></TD>";
print $query->end_form;


print "<TD><B><A HREF='help.matrix-scan.html'>MANUAL</A></B></TD>\n";
#print "<TD><B><A HREF='tutorials/tut_matrix-scan.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:Jacques.van-Helden\@univ-amu.fr'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

&ListParameters() if ($ENV{rsat_echo} >= 2);
&ListDefaultParameters() if ($ENV{rsat_echo} >= 2);

print $query->end_html;

exit(0);
