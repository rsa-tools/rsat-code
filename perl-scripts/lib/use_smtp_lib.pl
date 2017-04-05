## Load the dependency to Email::Sender::Transport::SMTP on ly if
## required since this module has only to be installed on some
## servers.

use Email::Sender::Transport::SMTP;
