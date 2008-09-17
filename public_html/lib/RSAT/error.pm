###############################################################
#
# A class for error handling
#
package RSAT::error;

use RSAT::GenericObject;
use RSAT::message;
@ISA = qw( RSAT::GenericObject );

=pod

=head1 NAME

    RSAT::error

=head1 DESCRIPTION

Error handling for RSAT.

=cut



################################################################
=pod

=item FatalError

Send an error message and die

=cut
sub FatalError {
    my @error_message = @_;
    my $message = join "\t", @error_message;
    $ENV{RSA_ERROR} = "1";
    my $context = $ENV{RSA_OUTPUT_CONTEXT} || "screen";
    if ($context eq "cgi") {
	$message =~ s/\n\t/ /g;
	$message =~ s/\t/ /g;
	&cgiError($message);
    } else {
      die("Error\n\t", $message, "\n");
    }
}

################################################################
=pod

=item cgiError

Print a HTML message and die.

=cut

sub cgiError {
  my $error_message = join "", @_;
  my $error_color = "#FF0000";
  $error_message .= "\n<p>(error occurred on host ".`hostname`.")";
  RSAT::message::cgiMessage($error_message, "Error", $error_color);
  print "</body></html>\n";
  exit(1);
}

return 1;


__END__

=pod

=back

