################################################################
#
# Help messages for describing formats or options that appear in
# several programs
#
################################################################

package main;

sub help_message {
    my ($help_request) = @_;
    my $help_mesage = "";
    
    if ($help_request eq "class file") {
	$help_message = "CLASS FILE
	The class file specifies the composition of several gene
	families.  It is a text file containing 2 columns separated by
	a tab character.

	    col 1:   class member
	    col 2:   class name

        Additional columns are ignored. 

	Lines starting with a semicolumn (;) are ignored, allowing to
	document the class files with comments.

        A given element (e.g. gene) can belong simultaneously to
        several families. In such a case, the element will appear on
        several rows (one per class),

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
	warn ("No help message about $help_request");
    }
    
    return $help_message;
}
    
return 1;
