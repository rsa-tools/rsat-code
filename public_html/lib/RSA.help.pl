################################################################
#
# Help messages for describing formats or options that appear in
# several programs
#
################################################################


sub help_message {
    my ($help_request) = @_;
    my $help_mesage = "";
    
    if ($help_request eq "family file") {
	$help_message = "FAMILY FILE
	The family file specifies the co;position of several gene families.
	It is a text file containing 2 tab-separated columns.
	col 1:	gene name (or ID)
	col 2: family name
	
	Lines starting with a semicolumn (;) are ignored, allowing to
	document the family files with comments..
	
	Example
		; genes responding to Phosphate stress
		pho5	PHO
		pho8	PHO
		; genes responding to nitrogen starvation
		DAL5    NIT
		GAP1    NIT
		...
";
    } else {
	&warn ("No help message about $help_request");
    }
    
    return $help_message;
}
    
