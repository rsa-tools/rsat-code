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
    $ENV{RSA_ERROR} = "1";
    my $context = $ENV{RSA_OUTPUT_CONTEXT} || "screen";
    if ($context eq "cgi") {
	&cgiError(@error_message);
    } else {
      my $message = join "\t", @error_message;
#      $message =~ s/\n\t/\n/g;
#      $message =~ s/\t+/ /g;
      die("Error\n\t", $message, "\n");
    }
}

################################################################
=pod

=item cgiError

Print a HTML message and die.

=cut

sub cgiError {
  my $error_message = join "<br>\n", @_;
  my $error_color = "#DD0000";
  my $hostname = `hostname`;
  chomp($hostname);
  $error_message .= "<br>\nError occurred on host ".$hostname;
  &RSAT::message::cgiMessage($error_message, "Error", $error_color);
  print "</body></html>\n";
  exit(1);
}

return 1;


__END__

=pod

=back

