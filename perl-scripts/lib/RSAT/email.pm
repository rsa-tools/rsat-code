package RSAT::email;

use Email::Simple;
use Email::Simple::Creator;
# check dispatch modules are actually installed; they are not in conda!
our $SEND_SIMPLE_AVAIL = eval "use Email::Sender::Simple; 1" ? 1 : 0;
our $SEND_SMTP_AVAIL = eval "use Email::Sender::Transport::SMTP; 1" ? 1 : 0;

=pod

=item B<CheckEmailAddress($email_address)>

check email address format

=cut
sub CheckEmailAddress {
    my ($email_address) = @_;
    if ($email_address eq "") {
	&RSAT::error::FatalError ("You did not enter your e-mail address");
    } if ($email_address =~ /http:\/\//) {
	&RSAT::error::FatalError ("Invalid email address: $email_address<br>", "Submitting http links in Web form is a hacking practice. This attempt will be reported.");

    } elsif ($email_address !~ /(\S+\@\S+)/) {
	&RSAT::error::FatalError ("Invalid email address: $email_address<br>");
    }
}



=pod

=item MessageToAdmin

Report an error by sending an email to RSAT administrator

=cut
sub MessageToAdmin {
    my ($message) = @_;

    ## Check if server admin has been specified
    unless (defined( $ENV{SERVER_ADMIN})) {
      &RSAT::message::Warning("Cannot send mail to server admin. Variable SERVER_ADMIN should be defined in RSAT_config.props");
      return();
    }

    ## Define title based on script name
    my $script_name = &RSAT::util::ShortFileName($0);
    my $title = join(" - " , "RSAT", $script_name, $date);
#    $mail_command = "mail -s \'".$title."\'";
#    $mail_command = "mail -s \'RSAT - $script_name - $date\'";
#    system "echo \"$message\" | $mail_command $ENV{SERVER_ADMIN} &"; 
    &send_mail($message, $ENV{SERVER_ADMIN}, $title);
}



=pod

=item B<Send an email from a TLS server, uses $ENV{mail_supported} and $SEND_SMTP_AVAIL>

=cut

sub send_mail_STARTTLS {
    my ($message, $recipient, $subject, $smtp_server, $user, $pass, $from, $warn_message) = @_;

    if ($ENV{mail_supported} eq 'no' || $SEND_SMTP_AVAIL == 0) {
	&RSAT::message::Warning("This RSAT Web site does not support email sending (TLS). ", $subject);
	return();
    } 

    ## Load the resquired module I cannot make the use() conditional
    ## so I load a 1-line perl file with require().
    require($ENV{RSAT}."/perl-scripts/lib/use_smtp_lib.pl");
	
    ## Check if recipient argument contains a valid email address
    &CheckEmailAddress($recipient);
    
    ## Set a subject if not specificed in arguents
    unless ($subject) {
	$script_name = $0;
	$subject = join " ", "[RSAT]", $script_name, &RSAT::util::AlphaDate();
    }
    
    ## Set the STARTTLS connection
    if (!$smtp_server) { $smtp_server = $ENV{starttls} }
    if (!$user) { $user = $ENV{starttls_user} }
    if (!$pass) { $pass = $ENV{starttls_pass} }  
    if (!$from) { $from = $ENV{smtp_sender} }
    
    ## Issue a warning to indicate that mail will be sent
    if (($ENV{rsat_echo} >= 1) || ($main::verbose >= 2)) {
	my $mail_warn = "Sending mail";
	$mail_warn .= " to \"".$recipient."\"" if ($recipient);
	$mail_warn .= " ; Subject: \"".$subject."\"";
	$mail_warn .= " SMTP STARTTLS server: ".$smtp_server if (($ENV{rsat_echo} >= 2) || ($main::verbose >= 2));
	&RSAT::message::TimeWarn($mail_warn) if ($warn_message);
    }
    
    ## Compose message
    my $email = Email::Simple->create(
	header => [
	    To      => $recipient,
        From    => $from,
	    Subject => $subject,
	],
	body => $message,
	);
    
    ## Try to send the email only if required module is available
    my $transport = Email::Sender::Transport::SMTP->new(
        host => $smtp_server,
        ssl  => 'starttls',
        sasl_username => $user,
        sasl_password => $pass,
        debug => 0, # or 1
    );
    
    eval { Email::Sender::Simple->send($email, {transport => $transport}) };
    &RSAT::error::FatalError ("Error sending email ", $@) if $@;

    return();
}


=pod

=item B<Send an email message, uses $ENV{mail_supported} and $SEND_SIMPLE_AVAIL>

=cut
sub send_mail {
  my ($message, $recipient, $subject, $warn_message) = @_;

  if ((defined($ENV{RSA_OUTPUT_CONTEXT})) &&
      (($ENV{RSA_OUTPUT_CONTEXT}eq "cgi") || ($ENV{RSA_OUTPUT_CONTEXT} eq "RSATWS"))) {
    if (($ENV{starttls} ne "") &&
	($ENV{starttls_user} ne "") &&
	($ENV{starttls_pass} ne "")) {
      &send_mail_STARTTLS($message, 
			  $recipient, 
			  $subject);
    }
  }

  if ($ENV{mail_supported} eq 'no' || $SEND_SIMPLE_AVAIL == 0) {
    &RSAT::message::Warning("This RSAT Web site does not support email sending (SIMPLE). ", $subject);
  } else {

    ## Check if recipient argument contains a valid email address
    &CheckEmailAddress($recipient);

    ## Set a subject if not specificed in arguents
    unless ($subject) {
      $script_name = $0;
      $subject = join " ", "[RSAT]", $script_name, &RSAT::util::AlphaDate();
    }

    ## Define the SMTP server
    my $smtp_server = "localhost"; ## Default is send by local machine
    if (($ENV{smtp}) && ($ENV{smtp} !~ /smtp.at.your.site/)) {
      $smtp_server = $ENV{smtp};
    }
    my $smtp_port = 25;
    if ($ENV{smtp_port}) {
	$smtp_port = $ENV{smtp_port};
      &RSAT::message::Debug("smtp_port", $smtp_port) if ($main::verbose >= 5);
    }
    &RSAT::message::Debug("smtp_server", $smtp_server) if ($main::verbose >= 5);

    ## Define the "from" email (can be defined in RSAT_config.props or
    ## as environment variable smtp_sender)
    my $from = "";
    if ($ENV{smtp_sender}) {
      $from = $ENV{smtp_sender};
    }

    ## Issue a warning to indicate that mail will be sent
    if (($ENV{rsat_echo} >= 1) || ($main::verbose >= 2)) {
      my $mail_warn = "Sending mail";
      $mail_warn .= " from \"".$from."\"" if ($from);
      $mail_warn .= " to \"".$recipient."\"" if ($recipient);
      $mail_warn .= " ; Subject: \"".$subject."\"";
      $mail_warn .= "SMTP server: ".$smtp_server if (($ENV{rsat_echo} >= 2) || ($main::verbose >= 2));
      &RSAT::message::TimeWarn($mail_warn) if ($warn_message);
    }

    # ## Send the message using MIME::Lite    
    # my $msg = MIME::Lite->new(
    # 	From    => $from,
    # 	To      => $recipient,
    # 	Subject => $subject." [MIME::Lite]",
    # 	Type    => 'text/plain',
    # 	Data    => $message,
    # 	);
    # $msg->send('smtp', $smtp_server);


    ## Sent message using Email::Sender
    my $email = Email::Simple->create(
      header => [
	To      => $recipient,
	From    => $from,
	Subject => $subject,
      ],
      body => $message,
	);
    &RSAT::message::Debug( "email", $email) if ($main::verbose >= 3);

    print "# $SEND_SIMPLE_AVAIL $SEND_SMTP_AVAIL\n";

#    &RSAT::message::Debug("INC", join (";", @INC)) if ($main::verbose >= 3);
    #    eval  {use Email::Sender::Transport::SMTP} ;  die $@ if $@;  ## Load the Perl module only if required
    &RSAT::message::Debug("smtp_server", $smtp_server) if ($main::verbose >= 5);
    
    my $transport;
    if ($smtp_port) {
	$transport = Email::Sender::Transport::SMTP->new({
	host => $smtp_server,
	port => $smtp_port});
    } else {
	$transport = Email::Sender::Transport::SMTP->new(host => $smtp_server);
    }
    Email::Sender::Simple->send($email, {transport => $transport});
  }
}

1;
