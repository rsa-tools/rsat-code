###############################################################
#
# Class TaskManager
#
package RSAT::TaskManager;

use RSAT::GenericObject;
use RSAT::error;
@ISA = qw( RSAT::GenericObject );

### class attributes

=pod

=head1 NAME

    RSAT::TaskManager

=head1 DESCRIPTION

TaskManager class. 

=cut

################################################################
## NOTE: THIS FUNCTION IS IN RSA.cgi
## It is doing specific treatment of the Web site, it is thus better placed in RSA.cgi than in the TaskManager
# ################################################################
# ## Send e-mail with the result of a command
# sub EmailTheResult {
#     my ($command, $email_address, $tmp_file_name, %args) = @_;
 
#    my $delay = "72 hours";
#     my $mail_title = $args{title};
#     unless ($mail_title) {
# 	my $script_name = $script_name || $0;
# 	$mail_title = join " ", "[RSAT]", $script_name, &AlphaDate();
#     }

#     ## check the email address
#     &CheckEmailAddress($email_address);

#     ## temporary file for storing the result
#     unless ($tmp_file_name) {
#       $tmp_file_name = join("", "result_", &AlphaDate(), ".txt");
#     }
# #    my $result_URL = "$ENV{rsat_www}/tmp/$tmp_file_name";
#     my $result_URL = "$tmp_file_name";
# #    $result_URL =~ s|$ENV{RSAT}/public_html/tmp||;
# #    $result_URL = $ENV{rsat_www}.$result_URL;
# #    $result_URL = join(";", "server config keys", keys %RSAT::server::config);
# #    $result_URL = join(";", "ENV", keys(%ENV));

#     ## debugging: report the command
#     print "<pre>$command > $TMP/$tmp_file_name</pre>" if ($ECHO >= 1);

#     ## Indicate the URL of the future result file
#     my $message = "The server is now processing your request.\n"; 
#     $message .= "Once it will be finished, the result will become available at the following URL\n";
#     $message .= "\t${result_URL}\n";
#     $message .= "When the result will be ready, you will be warned at your email address ($email_address).\n";
#     $message .= "The result file will remain on the server for $delay.\n";

# #    &SendMail($message, $email_address, $title);

#     ## Display a message announcing the result
#     $html_message = $message;
#     $html_message =~ s|(http://\S+)|<a href=$1>$1</a>|gm;
#     &RSAT::message::Info($html_message);

#     ## Prepare the email command
#     my $command =  "$command | perl -pe 's|$ENV{RSAT}/(public_html/)*||g' > $TMP/$tmp_file_name; ";

#     ## Send completion mail
#     my $done_message = "Your result is available at the following URL:\n\t${result_URL}";
#     $done_message .= "\nThe result file will remain there for $delay.";
#     $mail_command = "mail -s \'".$mail_title."\'";
#     $command .= "echo \"$done_message\" | $mail_command $email_address &"; 
#     print "<PRE>$command</PRE>" if ($ECHO >= 1);
#     system $command;

#     ## Prepare removal of the temporary file
#     &RSAT::server::DelayedRemoval("$TMP/$tmp_file_name", $delay) if (-e $TMP/$tmp_file_name);
# }


return 1;


__END__
