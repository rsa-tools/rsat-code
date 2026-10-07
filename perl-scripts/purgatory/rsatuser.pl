#!/usr/bin/perl -w

use strict;
use Getopt::Std;
use MIME::Lite;
use Email::Simple::Creator;
use Email::Sender::Transport::SMTP;
use Email::Sender::Simple qw(sendmail);

my $EMAILADD   = '';
my $MAILSERVER = '';
my $SMTPUSER   = ''; 
my $APACHEFILE = '';
my $LOGFILE    = '';
my $BOILERTEXT = "";
my $BOILERTEXT2= "";
my $HTPASSEXE  = '';

my ($user_email_file, $admin_pass, $test) = ('', '', 0);
my ($email, $pass, @emails, %opts);

getopts('htu:f:p:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)) {
  print "\nThis script registers new users to RSAT Web app\n";	
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-u user email address                            (required unless -f is set)\n";
  print "-f file with user addressess, 1 per line         (optional)\n";
  print "-p 'password' of SMTPUSER\@MAILSERVER            (required, see script source)\n";
  print "-t test mode, email sent to EMAILADD             (EMAILADD=$EMAILADD)\n\n\n";
  print "Example:\n\n  sudo perl rsatuser.pl -f emails.txt -p 'F...'\n";
  exit(0);
}

if(defined($opts{'f'})) {
  $user_email_file = $opts{'f'};	
  open(LIST,"<",$user_email_file) ||
    die "# ERROR: cannot read $user_email_file\n";
  while(<LIST>) {
    if(/^([^\@]+\@\S+)/) {
      $email = $1;
      push(@emails, $email);
    }
  }
  close(LIST);  

} elsif(defined($opts{'u'})) {
  push(@emails, $opts{'u'});
} 

if(defined($opts{'p'})) {
  $admin_pass = $opts{'p'};

} else {
  die "# ERROR: need a valid -p password to send notification email\n";
}	

if(defined($opts{'t'})) {
  $test = 1
} else {
  open(LOG,">>",$LOGFILE) ||
    die "# ERROR: cannot append to $LOGFILE\n";
}

# assign password to emails, register users and send email notification
foreach $email (@emails) {
  print "# $email\n";

  # create new password, https://stackoverflow.com/a/801361
  $pass = join('', map +(0..9,'a'..'z','A'..'Z')[rand(10+26*2)], 1..30);

  if($test == 1) {

    # don't register, test only
    print "$email ".localtime()."\n";
    $email = $EMAILADD

  } else {
    # register user in Apache file	  
    system("$HTPASSEXE -b $APACHEFILE $email $pass"); 
    if ( $? != 0 ) {
      die "# ERROR: failed running $HTPASSEXE $APACHEFILE $email $pass\n";
    } else {  	    
      print LOG "$email ".localtime()."\n";
    } 
  }

  # send notification
  send_email($MAILSERVER,$SMTPUSER,$admin_pass,
    $EMAILADD,$email,'RSAT registration',
    "\n$BOILERTEXT\n\nuser: $email\npassword: $pass\n\n$BOILERTEXT2");
}

if($test == 1) {
  exit(0);

} else {
  close(LOG);
  exit(0);
}




sub send_email
{
  my ($host,$huser,$hpass,$sender,$address,$subject,$text) = @_;

  my $transport = Email::Sender::Transport::SMTP->new(
    host => $host,
    ssl  => 'starttls',
    sasl_username => $huser,
    sasl_password => $hpass,
   debug => 0, # or 1
  );

  my $email = Email::Simple->create(
    header => [
      From    => $sender,
      To      => $address,
      Subject => $subject,
    ],
    body => $text,
  );   

  sendmail($email, {transport => $transport});
}

