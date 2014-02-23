################################################################
##
## A class for message handling
##
package RSAT::message;

use RSAT::GenericObject;
@ISA = qw( RSAT::GenericObject );


=pod

=head1 NAME

    RSAT::message

=head1 DESCRIPTION

Message handling for RSAT.

=cut



=pod

=item new()

Create a new message.

=cut
sub new {
    my ($class, %args) = @_;
    my $message = bless {
	}, $class;
    return $message;
}


=pod

=item Warning(@warninng_message)

Print a warning message (on the STDERR). The message can be provided
as a list. Each element of the list is then printed on a separate
line.

=cut
sub Warning {
    my @warning_message = @_;
    my $message = join "\t", @warning_message;
    if ((defined($ENV{RSA_OUTPUT_CONTEXT})) &&
	($ENV{RSA_OUTPUT_CONTEXT} eq "cgi")) {
	$message =~ s/\n/<br>\n/g;
	&cgiWarning($message);
    } else {
	warn("; WARNING\t", $message, "\n");
    }
}



=pod

=item cgiWarning

Print a warning message in HTML format (STDOUT)

=cut
sub cgiWarning {
    my $warning_message = join "<P>", @_;
    my $warning_color = "#FFAA00";
    &cgiMessage($warning_message, "Warning", $warning_color);
}


=pod

=item Info

Information message

=cut
sub Info {
    my (@info_message) = @_;
    my $message = join "\t", @info_message;
    if (defined($ENV{RSA_OUTPUT_CONTEXT}) 
	&& ($ENV{RSA_OUTPUT_CONTEXT} eq "cgi")) {
	$message =~ s/\n/<br>\n/g;
	&cgiMessage($message);
    } else {
	warn("; INFO\t", $message, "\n");
    }
}


=pod

=item Debug

Debug message

=cut
sub Debug {
    my (@debug_message) = @_;
    my $message = join "\t", @debug_message;
    if ((defined($ENV{RSA_OUTPUT_CONTEXT})) && ($ENV{RSA_OUTPUT_CONTEXT} eq "cgi")) {
	$message =~ s/\n/<br>\n/g;
	&cgiMessage($message);
    } else {
	warn("; DEBUG\t", $message, "\n");
    }
}


=pod

=item cgiMessage

Print a message in HTML format (STDOUT)

=cut
sub cgiMessage {
  my ($message, $message_type, $color) = @_;
  $color = "#006600" unless ($color);
  $message_type = "Information" unless ($message_type);
  print  ("<blockquote class='",lc($message_type),"'>",
	  "\n",
	  "<font color='".$color."'><b>",$message_type,": </b>",
	  $message,
	  "</font>\n",
	  "</blockquote>",
	  "<br><hr size=3>\n");
}




=pod

=item TimeWarn

Warning with time 

=cut
sub TimeWarn {
    my $message = "";
    $message .= join( "\t", &RSAT::util::AlphaDate(), @_);
    chomp($message);
    if (defined($ENV{RSA_OUTPUT_CONTEXT}) 
	&& ($ENV{RSA_OUTPUT_CONTEXT} eq "cgi")) {
	$message =~ s/\n/<br>\n/g;
	&cgiMessage($message, "TimeWarn");
    } else {
	$message .= "\n";
	warn ("; ".$message);
    }
}


=pod

=item psWarn

Warning with details about hte current Perl process (memory usage, cpu usage,
...). This is useful to track problems during the execution of perl scripts.

=cut
sub psWarn {
    my (@message) = @_;
    my $message = join "\t", @message;

    ## Get information on the current process
    my $pid = $$;
    my $ps_cmd = "ps -p $pid -O '%mem %cpu size rss vsize'";
    my $ps = `$ps_cmd`;
    $message .= "\n";
    $message .= $ps;
    $message .= "\n";
    &TimeWarn($message);
#    warn ($message);
}


=pod

=item MessageToAdmin

Report an error by sending an email to RSAT administrator

=cut
sub MessageToAdmin {
    my ($message) = @_;
    my $script_name = &RSAT::util::ShortFileName($0);
    $mail_command = "mail -s \'RSA-tools - $script_name - $date\'";
    system "echo \"$message\" | $mail_command $ENV{SERVER_ADMIN} &"; 
}

return 1;


__END__


=pod

=back

