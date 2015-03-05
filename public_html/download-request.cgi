#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
  push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
#### redirect error log to a file
BEGIN {
  $ERR_LOG = "/dev/null";
#    $ERR_LOG = &RSAT::util::get_pub_temp()."/RSA_ERROR_LOG.txt";
  use CGI::Carp qw(carpout);
  open (LOG, ">> $ERR_LOG")
      || die "Unable to redirect log\n";
    carpout(*LOG);
}
require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

## Initialize parameters
my @mandatory_fields  = qw ( first_name
                             last_name
                             email_address
                             institution
                             city
                             country
                             );

my @all_fields = @mandatory_fields;

### Read the CGI query
$query = new CGI;

### Print the header
&RSA_header("Download", "results");

## Check security issues
&CheckWebInput($query);

## Check mandatory fields
foreach my $field (@mandatory_fields) {
  if ($query->param($field) eq "") {
    &RSAT::error::FatalError($field." field cannot be empty.");
  }
}

## Check email
my $email_address = "";
if ($query->param('email_address') =~ /(\S+\@\S+\.\S+)/) {
  $email_address = $1;
} else {
  &RSAT::error::FatalError("Invalid email address", $query->param('email_address'));
}

## Send mail
my $message = "RSAT download request\n\n";
foreach my $field (@all_fields) {
  $message .= sprintf("%-22s\t%s", $field, $query->param($field));
  $message .= "\n";
}

print "<pre>";
print $message;
print "</pre>";
my $recipient = 'Jacques.van-Helden@univ-amu.fr'; ## All download requests should be sent to JvH
my $subject = join(" ", 
		   "RSAT download request from",
		   $query->param("first_name"), 
		   $query->param("last_name"), 
    );
&RSAT::server::send_mail($message, $recipient, $subject);

## Indicate download URL
print "<h2>Download URL</h2>";

print "<p>To download the Regulatory Sequence Analysis Tools, please follow this link</p>\n";

my $download_url = "http://download.rsat.eu/";
print "<p>", "<b><a href='",$download_url, "'>",$download_url, "</a></b>\n"; 

print $query->end_html();

exit(0);


