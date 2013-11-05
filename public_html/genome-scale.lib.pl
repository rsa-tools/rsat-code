### default values for sequence retrieval
$default{sequence_format} = "fasta";
$default{seq_label} = "name";
$default{organism} = "Saccharomyces cerevisiae";
$default{from} = "default";
$default{to} = "default";
$default{feattype} = "CDS";
$default{sequence_type} = "upstream";


################################################################
#
# retrieve-seq options
#
sub DisplayRetrieveSeqOptions {
    print $query->h4("Sequence retrieval options");

    print $query->hidden(-name=>'genes',-default=>"all");
    print $query->hidden(-name=>'sequence_format',-default=>$default{sequence_format});


    &OrganismPopUp();

    #### feature type
    print "<B><A HREF='help.retrieve-seq.html#feattype'>Feature type</A></B>&nbsp;";
    print $query->radio_group(-name=>'feattype',
			      -values=>[@supported_feature_types],
			  -default=>$default{feattype});
    print "<BR>\n";

    ### sequence type
    print "<B><A HREF='help.retrieve-seq.html#sequence_type'>Sequence type</A></B>&nbsp;";
    print $query->popup_menu(-name=>'sequence_type',
			     -Values=>['upstream','downstream','ORF (unspliced)','whole chromosomes'],
			     -default=>$default{sequence_type});

    ### from to
    print "<B><A HREF='help.retrieve-seq.html#from_to'>From</A></B>&nbsp;\n";
    print $query->textfield(-name=>'from',
			    -default=>$default{from},
			    -size=>10);

    print "&nbsp;&nbsp;";
    print "<B><A HREF='help.retrieve-seq.html#from_to'>To</A></B>&nbsp;\n";
    print $query->textfield(-name=>'to',
			    -default=>$default{to},
			    -size=>10);
    print "<BR>\n";


    ### allow ORF overlap
    print $query->checkbox(-name=>'orf_overlap',
			 -checked=>'checked',
			 -label=>'');
    print "&nbsp;<A HREF='help.retrieve-seq.html#noorf'><B>allow overlap with upstream ORFs</B></A>";
    print "<BR>\n";

### temporarily inactivated because it does not work with all organisms
#    print $query->hidden(-name=>'orf_overlap',-default=>'on');

    ### sequence label
    print "<B><A HREF='help.retrieve-seq.html#seq_label'>Sequence label</A></B>&nbsp;";
    print $query->popup_menu(-name=>'seq_label',
			     -Values=>['gene identifier', 
				       'gene name',
				       'gene identifier + name',
				       'gene identifier + organism + gene name',
				       'full identifier'
				      ],
			     -default=>$default{seq_label});
    print "<BR>\n";

    ### sequence label
#    print $query->hidden(-name=>'seq_label',-default=>$default{seq_label});

    print "<BR>\n";
}

################################################################
#
# retrieve-seq parameters
#
sub ReadRetrieveSeqParams {
    $org = $query->param('organism');
    if ($query->param('sequence_type') =~ /chromosome/i) {
	$retrieve_seq_command = "$SCRIPTS/convert-seq";

	#### take whole genome as input file
	$retrieve_seq_parameters .= " -i $supported_organism{$org}->{'genome'}";
	$retrieve_seq_parameters .= " -from $supported_organism{$org}->{'seq_format'}";
	
	### output format ###
	if ($accepted_output_seq{$query->param('sequence_format')}) {
	    $retrieve_seq_parameters .= " -to ".$query->param('sequence_format');;
	}
	

    } else {
	$retrieve_seq_command = "$SCRIPTS/retrieve-seq";
	$retrieve_seq_parameters = " -all -nocomment";
	
	#### organism
	if (defined($supported_organism{$query->param('organism')})) {
	    $org = $query->param('organism');
	} else {
	    $org = "Saccharomyces_cerevisiae";
	}
	$retrieve_seq_parameters .= " -org ".$org;
	
	### feature type
	if ($query->param('feattype')) {
	    my ($feattype) = split " ", $query->param('feattype'); ### take the first word
	    $retrieve_seq_parameters .= " -feattype ".$feattype;
	}
	
	### sequence type
	if ($query->param('sequence_type')) {
	    $retrieve_seq_parameters .= " -type ".$query->param('sequence_type');
	}

	### output format ###
	if ($accepted_output_seq{$query->param('sequence_format')}) {
	    $retrieve_seq_parameters .= " -format ".$query->param('sequence_format');;
	}



	### sequence label
	my $seq_label = lc($query->param('seq_label'));
#	$retrieve_seq_parameters .= " -label ".$seq_label;
  	if (($seq_label =~ /name/) && 
	    ($seq_label =~ /organism/) && 
 	    ($seq_label =~ /identifier/)) {
 	    $retrieve_seq_parameters .= " -label id,organism_name,name";
	  } elsif (($seq_label =~ /name/) && 
		   ($seq_label =~ /identifier/)) {
 	    $retrieve_seq_parameters .= " -label id,name";
 	} elsif ($seq_label =~ /name/) {
 	    $retrieve_seq_parameters .= " -label name";
 	} elsif ($seq_label =~ /identifier/) {
 	    $retrieve_seq_parameters .= " -label id";
 	} elsif ($seq_label =~ /full/) {
 	    $retrieve_seq_parameters .= " -label full";
 	} else {
 	    &cgiError("Invalid option for sequence label '$seq_label'");
 	}

	### limits ###
	if (&IsInteger($query->param('from'))) {
	    $retrieve_seq_parameters .= " -from ".$query->param('from');
	}  
	if (&IsInteger($query->param('to'))) {
	    $retrieve_seq_parameters .= " -to ".$query->param('to');
	}

	### orf overlap ###
	unless (lc($query->param('orf_overlap')) eq "on") {
	    $retrieve_seq_parameters .= " -noorf ";
	}
    }
    

    ## return command and parameters
    return($retrieve_seq_command, $retrieve_seq_parameters);
}


return 1;

