#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
$pattern_assembly_command = "$SCRIPTS/pattern-assembly";
$tmp_file_name = sprintf "pattern-assembly.%s", &AlphaDate;
$result_file = "$TMP/$tmp_file_name.res";

### Read the CGI query
$query = new CGI;

#### update log file ####
&UpdateLogFile();

### Print the header
&RSA_header("pattern-assembly result");
&ListParameters if ($ECHO >=2);

#### read parameters ####
$parameters .= " -v 1";

################################################################
#### patterns
if ( $query->param('patterns') =~ /\S/) {
    open QUERY, ">$TMP/$tmp_file_name.pat";
    print QUERY $query->param('patterns');
    close QUERY;
    &DelayedRemoval("$TMP/$tmp_file_name.pat");
    $parameters .= " -i $TMP/$tmp_file_name.pat";
#} elsif ($query->param('pattern_file') =~ /\S/) {
} else {
    &cgiError("You should either enter query patterns in the box, or specify a pattern file\n");
}  

################################################################
####  max number of flanking residues
$maxfl = $query->param('maxfl');
if (&IsNatural($maxfl)) {
    if ($maxfl < 0) {
	&cgiError("Maximum flanking residues should be positive ($maxfl is invalid)");
    } else {
	$parameters .= " -maxfl $maxfl";
    }
} else {
    &cgiError("Maximum flanking residues should be a positive natural number ($maxfl is invalid)");
}

################################################################
####  max number of substitutions
$subst = $query->param('subst');
if (&IsNatural($subst)) {
    if ($subst < 0) {
	&cgiError("Maximum substitutions should be positive ($subst is invalid)");
    } else {
	$parameters .= " -subst $subst";
    }
} else {
    &cgiError("Maximum substitutions should be a positive natural number ($subst is invalid)");
}

################################################################
####  max number of patterns
$maxpat = $query->param('maxpat');
if (&IsNatural($maxpat)) {
    if ($maxpat < 0) {
	&cgiError("Maximum number of patterns should be positive ($maxpat is invalid)");
    } else {
	$parameters .= " -maxpat $maxpat";
    }
} else {
    &cgiError("Maximum number of patterns should be a positive natural number ($maxpat is invalid)");
}

################################################################
####  max number of patterns
$sc = $query->param('sc');
if (&IsNatural($sc)) {
    if ($sc < 2) {
	&cgiError("Score column should be >= 2 ($sc is invalid)");
    } else {
	$parameters .= " -sc $sc";
    }
}

################################################################
#### single or double strand assembly
if ($query->param('strand') =~ /single/) {
  $parameters .= " -1str";
} else {
  $parameters .= " -2str";
}


################################################################
#### run the command
print "<PRE>command: $pattern_assembly_command $parameters<P>\n</PRE>" if ($ECHO >=1);
if ($query->param('output') eq "display") {
#    &PipingWarning();

    open RESULT, "$pattern_assembly_command $parameters |";

    print '<H2>Result</H2>';
    &PrintHtmlTable(RESULT, $result_file, true);
    close(RESULT);

    &PipingForm();

    print "<HR SIZE = 3>";
} else { 
    &EmailTheResult("$pattern_assembly_command $parameters", $query->param('user_email'));
}
print $query->end_html();

exit(0);


################################################################
#
# Pipe the result to other commands
#
sub PipingForm {
    my $genes = `cat $result_file`;
    ### prepare data for piping
    print <<End_of_form;
<HR SIZE = 3>
<TABLE>
<TR>
<TD>
<H3>Next step</H3>
</TD>
<TD>
<FORM METHOD="POST" ACTION="retrieve-seq_form.cgi">
<INPUT type="hidden" NAME="organism" VALUE="$organism">
<INPUT type="hidden" NAME="genes" VALUE="selection">
<INPUT type="hidden" NAME="gene_selection" VALUE="$genes">
<INPUT type="submit" value="retrieve sequences">
</FORM>
</TD>
</TR>
</TABLE>
End_of_form

}

