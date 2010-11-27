################################################################
#### Library for patser CGI interface


################################################################
### default values for filling the form
$default{matrix_format} = "consensus";
$default{matrix} = ""; ### [-m <Name of matrix file---default name is "matrix">]
$default{matrix_is_weight} = ""; ### [-w <Matrix is a weight matrix>]
$default{matrix_is_vertical} = ""; ### [-v <Vertical matrix---rows correspond to positions>]
$default{pseudo_counts} = 1; ### [-b <Correction added to the elements of the alignment matrix (default: 1)>]
$default{alphabet_file} = ""; ### [-a <Name of ascii alphabet file---default name is "alphabet">]
$default{alphabet} = "a:t 0.3 c:g 0.2"; ### [-A <Ascii alphabet information>]
$default{case} = "insensitive"; ### [-CS <Ascii alphabet is case sensitive (default: ascii alphabets are case insensitive)>]
$default{strands} = "both"; ### [-c <Score the complementary strand>]

$default{lthreshold_method} = "weight score";
$default{lthreshold} = "7"; ### [-ls <Lower-threshold score, inclusive (formerly the -l option)>]
$default{uthreshold} = "none"; ### [-u <Upper-threshold score, exclusive>]


$default{return} = "all matches";
$default{output_format} = "feature list";
$default{table_score} = "checked";
$default{table_start} = "";
$default{table_end} = "";
$default{table_strand} = "";
$default{table_position} = "";
$default{top_scores} = "3"; ### [-t <Print only the top scores>]
$default{positions} = "checked"; ### convert the result into a score table
$default{table} = ""; ### convert the result into a score table
$default{sort} = "checked"; ### [-ds <Print top scores in order of decreasing score (default: print in order of position)>]

$default{unrecognized} = "discontinuities (with warning)"; ### [-d1 <Treat unrecognized characters as discontinuities, but print warning (the default)>]

$default{vertically_print} = "checked"; ### [-p <Vertically print the weight matrix>]
$default{min_calc_P} = 0; # [-M <Set the minimum score for calculating the p-value of scores (default: 0)>]

#### additional options
$default{origin} = "end";
$default{flanking} = "4";

################################################################
#
# Display patser  options
#
sub DisplayPatserOptions {
    
    print $query->h4('Patser options');
    
    ################################################################
    #### Matrix specification
    print "<A HREF='help.patser.html#matrix'><B>\n";
    print "Matrix</B></A>\n";
    print "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\n";
    print "<B>Format</B>&nbsp;";
    
    #### matrix format
    print $query->popup_menu(-name=>'matrix_format',
			     -Values=>['tab', 
				       'consensus',
				       'gibbs',
				       'transfac'
				       ],
			     -default=>$default{matrix_format});
    
    #### weight matrix
    print "&nbsp;"x6;
    print $query->checkbox(-name=>'matrix_is_weight',
			   -label=>" contains weigths",
			   -checked=>$default{matrix_is_weight});

    #### vertical matrix
    print "&nbsp;"x6;
    print $query->checkbox(-name=>'matrix_is_vertical',
			   -label=>" vertical",
			   -checked=>$default{matrix_is_vertical});


    #### text area to enter the matrix
    print "<BR>\n";
    print $query->textarea(-name=>'matrix',
			   -default=>$default{matrix},
			   -rows=>4,
			   -columns=>60);

    
    ################################################################
    #### strands
    print "<BR>\n";
    print "<A HREF='help.patser.html#strands'><B>Search strands</B></A>&nbsp;\n";
    print $query->popup_menu(-name=>'strands',
			     -Values=>['single',
				       'both'],
			     -default=>$default{strands});

    ################################################################
    #### return value
    print "<BR>\n";
    print "<table border=0>\n";
    print "<tr><td>\n";
    print "<A HREF='help.patser.html#return'><B>Return</B></A>&nbsp;\n";
    print "</td><td>\n";
    print $query->radio_group(-name=>'return',
			      -Values=>[ 'all matches'],
			      -default=>$default{return});

    print "</td><td>\n";
    print $query->radio_group(-name=>'return',
			      -Values=>['top values for each sequence'],
			      -labels=>{'top values for each sequence'=>''},
			      -default=>$default{return});
    print $query->textfield(-name=>top_scores,
			    -default=>$default{top_scores},
			    -size=>3
			    );
    print "&nbsp;top value(s) for each sequence";


    print "</td></tr><tr><td>\n";
    print "<A HREF='help.patser.html#output_format'><B>Output format</B></A>&nbsp;\n";
    print "</td><td>\n";
    print $query->radio_group(-name=>'output_format',
			      -Values=>['feature list'],
			      -label=>'feature list',
			      -default=>$default{output_format});

    print "</td></tr><tr><td>\n";
    print "</td><td>\n";
    print $query->radio_group(-name=>'output_format',
			      -Values=>['table'],
			      -label=>'table',
			      -default=>$default{output_format});

    print "</td><td>\n";
    print $query->checkbox(-name=>'table_score',
			   -label=>' score',
			   -checked=>$default{table_score});
    print $query->checkbox(-name=>'table_start',
			   -label=>' start',
			   -checked=>$default{table_start});
    print $query->checkbox(-name=>'table_end',
			   -label=>' end',
			   -checked=>$default{table_end});
    print $query->checkbox(-name=>'table_strand',
			   -label=>' strand',
			   -checked=>$default{table_strand});
    print $query->checkbox(-name=>'table_position',
			   -label=>' position',
			   -checked=>$default{table_position});

    print "</td></tr>\n";
    print "</table>\n";
    

    ################################################################
    #### pseudo-counts
    print "<BR>\n";
    print "<B><A HREF='help.patser.html#pseudo_counts'>Pseudo-counts</A>\n";
    print $query->textfield(-name=>'pseudo_counts',
			    -default=>$default{pseudo_counts},
			    -size=>2);
    
    #### TEMPORARILY DISACTIVATED, BECAUSE INTERFERES WITH ADJUSTED INFO THRESHOLD
    ################################################################
    ### [-M <Set the minimum score for calculating the p-value of scores (default: 0)>]
    # print "<BR>\n";
    # print "<B><A HREF='help.patser.html#min_calc_P'>Minimum score for calculating the p-value</A>\n";
    # print $query->textfield(-name=>'min_calc_P',
    # 			-default=>$default{min_calc_P},
    # 			-size=>5);
    
    ################################################################
    #### Thresholds
    
    #### lower threshold
    print "<br>\n";
    print "<A HREF='help.patser.html#lthreshold'><B>Lower threshold estimation</B></A>";
    print $query->popup_menu(-name=>'lthreshold_method',
			     -Values=>['weight score', 
				       'maximum ln(p-value)', 
				       'adjusted information content (auto)', 
				       'none'],
			     -default=>$default{lthreshold_method});
    
    print $query->textfield(-name=>'lthreshold',
			    -default=>$default{lthreshold},
			    -size=>6);
    
    
    #### upper threshold
    print "<br>\n";
    print "<A HREF='help.patser.html#uthreshold'><B>Upper threshold</B></A>", "&nbsp"x6;
    print $query->textfield(-name=>'uthreshold',
			    -default=>$default{uthreshold},
			    -size=>6);
    
    ################################################################
    #### alphabet
    print "<BR>\n";
    print "<B><A HREF='help.patser.html#alphabet'>\n";
    print "Alphabet</A></b>\n";
    print $query->textfield(-name=>'alphabet',
			    -default=>$default{alphabet},
			    -size=>50);
    
    ################################################################
    #### case sensitivity
    print "<BR>\n";
    print "<B><A HREF='help.patser.html#case'>\n";
    print "Case</A></b>\n";
    print "&nbsp;"x6;
    print $query->popup_menu(-name=>'case',
			     -Values=>[ "sensitive", "insensitive","insensitive, but mark lowercases"],
			     -default=>$default{case});
    
    ################################################################
    #### unrecognized characters
    print "<BR>\n";
    print "<B><A HREF='help.patser.html#unrecognized'>\n";
    print "Treat unrecognized characters as</A></b>\n";
    print "&nbsp;"x2;
    print $query->popup_menu(-name=>'unrecognized',
			     -Values=>[ "errors", "discontinuities (with warning)","discontinuities (no warning)"],
			     -default=>$default{unrecognized});

    ################################################################
    #### vertically print the matrix
    print "<BR>\n";
    print "<a href=help.patser.html#vertically_print>";
    print $query->checkbox(-name=>'vertically_print',
			   -label=>' print the weight matrix',
			   -checked=>$default{vertically_print});
    print "</a>\n</br>";
}




################################################################
#
# read patser parameters
#
sub ReadPatserParameters {
    $matrix_format = lc($query->param('matrix_format'));
    $matrix_from_transfac_command = $SCRIPTS."/matrix-from-transfac";
    $matrix_from_gibbs_command = $SCRIPTS."/matrix-from-gibbs";
    
    $convert_matrix_command = $SCRIPTS."/convert-matrix -from ".$matrix_format." -to tab -return counts";

    #### alphabet
    my $alphabet = $query->param('alphabet') || " a:t c:g ";
    $patser_parameters = " -A $alphabet";


    ################################################################
    #### matrix specification
    unless ($query->param('matrix') =~ /\S/) { ### empty matrix
	&cgiError("You did not enter the matrix");
    }

    $matrix_file = "$TMP/$tmp_file_name.matrix";

    if ($matrix_format =~ /tab/i) {
	open MAT, "| $convert_matrix_command | grep -v '\/\/' | perl -pe 's|\t|   |g' > $matrix_file";
    } elsif ($matrix_format =~ /transfac/i) {
	open MAT, "| $matrix_from_transfac_command > $matrix_file";
    } elsif ($matrix_format =~ /gibbs/i) {
	open MAT, "| $matrix_from_gibbs_command > $matrix_file";
    } elsif ($matrix_format =~ /consensus/i) {
	open MAT, "> $matrix_file";
    } else {
	&cgiError("Invalid matrix format.");
    }
    print MAT $query->param('matrix');
    close MAT;
    &DelayedRemoval($matrix_file);
    $patser_parameters .= " -m $matrix_file";

    #### [-w <Matrix is a weight matrix>]
    if ($query->param('matrix_is_weight')) {
	$patser_parameters .= " -w";
    } else {
	#### pseudo-counts and weights are mutually exclusive
	if (&IsReal($query->param('pseudo_counts'))) {
	    $patser_parameters .= " -b ".$query->param('pseudo_counts');
	}

    }

    #### [-v <Vertical matrix---rows correspond to positions>]
    if ($query->param('matrix_is_vertical')) {
	$patser_parameters .= " -v";
    }
    ################################################################
    #### strands 
    if ($query->param('strands') =~ /both/i) {
	$patser_parameters .= " -c";
    }

    ################################################################
    #### return top values
    if ($query->param('return') =~ /top/i) {
	$patser_parameters .= " -t";
	$top_scores = $query->param('top_scores');
	if (&IsNatural($query->param('top_scores'))) {
	    if ($top_scores == 0) {
		&FatalError("number of top scores must be >= 1");
	    } else {
		$patser_parameters .= " $top_scores";
	    }
	} else {
	    &FatalError("Number of top scores must be a strictly positive integer");
	}
    }

    ################################################################
    #### case sensitivity
    if ($query->param('case') eq "sensitive") {
	$patser_parameters .= ' -CS'; #### [-CS <Ascii alphabet is case sensitive (default: ascii alphabets are case insensitive)>]
    } elsif ($query->param('case') =~ /mark/) {
	$patser_parameters .= ' -CM'; #### [-CM <Ascii alphabet is case insensitive, but mark the location of lowercase letters>]
    }

    ################################################################
    #### unrecognized characters
    if ($query->param('unrecognized') eq "errors") {
	$patser_parameters .= " -d0";
    } elsif ($query->param('unrecognized') =~ /no warning/) {
	$patser_parameters .= " -d2";
    } else {
	$patser_parameters .= " -d1";
    }

    ################################################################
    #### thresholds ###

    #### lower threshold on the weight
    if ($query->param('lthreshold_method') =~ /weight score/) {
	if (&IsReal($query->param('lthreshold'))) {
	    $patser_parameters .= " -ls ".$query->param('lthreshold');
	} elsif ($query->param('lthreshold') eq 'none') {
	    ### no lower threshold
	} else {
	    &Warning("Lower threshold ignored (not a real value)");
	}
	
	#### lower threshold on P-value
    } elsif  ($query->param('lthreshold_method') =~ /p\-value/) {
	### [-lp <Determine lower-threshold score from a maximum ln(p-value)>]
	if (&IsReal($query->param('lthreshold'))) {
	    $patser_parameters .= " -lp ".$query->param('lthreshold');
	} elsif ($query->param('lthreshold') eq 'none') {
	    ### no lower threshold
	} else {
	    &FatalError ("Lower threshold value must be a real number");
	}

	#### automatic threshold on the basis of adjusted information content
    } elsif  ($query->param('lthreshold_method') =~ /adjusted information content/) {
	if ($query->param('lthreshold') eq 'auto') {
	    ### [-li <Determine lower-threshold score from adjusted information content>]
	    $patser_parameters .= ' -li';
	} elsif ($query->param('lthreshold') eq 'none') {
	    ### no lower threshold
	} else {
	    &Warning ("Invalid value (".$query->param('lthreshold').") for adjusted information content (supported: 'auto' or 'none'). The value will be ignored");
	}

	#### no lower threshold
    } elsif  ($query->param('lthreshold_method') =~ /none/) {

    } else {
	&FatalError("Unknown method for estimating lower threshold");
    }

    #### upper threshold
    if (&IsReal($query->param('uthreshold'))) {
	$patser_parameters .= " -u ".$query->param('uthreshold');
    }


    #### TEMPORARILY DISACTIVATED, BECAUSE INTERFERES WITH ADJUSTED INFO THRESHOLD
    ################################################################
    #### minimum score for calculating the P-value
    # if (&IsReal($query->param('min_calc_P'))) {
    #     $patser_parameters .= " -M ".$query->param('min_calc_P');
    # } else {
    #     &FatalError("Minimum score for calculating P-value must e a real number");
    # }
    
}



################################################################
### parameters for the piping to the feature map ###
#$feature_file =  "$TMP/$tmp_file_name.ft";
sub ReadFeaturesFromPatserParams {

    #$feature_file =  "$TMP/$tmp_file_name.ft";
    $features_from_patser_cmd .= " -seq $sequence_file ";

    #### origin 
    if ($query->param('origin') =~ /end/i) {
	$features_from_patser_cmd .= " -origin -0";
    }

    #### flanking residues for the matching sequences
    if ($query->param('flanking') =~ /^\d+$/) {
	$features_from_patser_cmd .= " -N ".$query->param('flanking');
    }
    #$features_from_patser_cmd .= " -o $feature_file";
    
    &ReadPatserTableOutputFields();

}

################################################################
#### output fields for the table
sub ReadPatserTableOutputFields {
    ### return table with one row per sequence
    if ($query->param('output_format') eq 'table') {
	push @table_fields, 'score' if ($query->param('table_score'));
	push @table_fields, 'start' if ($query->param('table_start'));
	push @table_fields, 'end' if ($query->param('table_end'));
	push @table_fields, 'strand' if ($query->param('table_strand'));
	push @table_fields, 'position' if ($query->param('table_position'));
	$table_fields = join ",", @table_fields;
	$features_from_patser_cmd .= " -table $table_fields";
    }
}

################################################################
### if a matrix file is specified in the query,
### read matrix from this file
sub ReadMatrixFromFile {
    if (($matrix_file = $query->param("matrix_file")) &&
	(-e $matrix_file)) {
	open MATRIX, $matrix_file;
	while (<MATRIX>) {
	    $default{matrix} .= $_;
	}
	close MATRIX;
    }
}

return 1;

