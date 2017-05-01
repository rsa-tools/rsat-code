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

### Read the CGI query
$query = new CGI;

### print the header
&RSA_header("Supported organisms", "results");

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

&ListParameters() if ($ENV{rsat_echo} >= 2);


##$tmp_file_name = sprintf "supported-organisms.%s", &AlphaDate();
#$prefix = "supported-organisms";
#$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); ($tmp_file_dir, $tmp_file_name) = &SplitFileName($tmp_file_path);

$font{variable} = 1;
#$command = "$SCRIPTS/supported-organisms -v 1";

my @return_fields = qw (ID taxid nb source last_update up_from taxonomy data);
my $return_fields = join",", @return_fields;
#$parameters = " -return ".$return_fields;

################################################################
## treat taxon specificity of the server if required
my $group = "";
if ($ENV{group_specificity}) {
  $group = $ENV{group_specificity};
}

## Export the table with header and absolute paths
my $organism_table = &RSAT::OrganismManager::supported_organism_table(1, 0, $source, $taxon, $group, $depth, @return_fields); 
my @organism_rows = split("\n", $organism_table);
my $nb_organisms = scalar(@organism_rows);

## Check if at least one organism is supported
if ($nb_organisms == 0) {
  &RSAT::error::FatalError("Not a single organism supported on this server.");
}


################################################################
## Print general information about this RSAT instance
print "<h2>RSAT instance: ", $ENV{rsat_site}, "</h2>\n";

print "<p><b>Organisms supported: </b>", $nb_organisms, "</p>\n";

if ($group) {
  print "<p><b>Group specificity: </b>", $group, "</p>\n";
}

################################################################
## Print the table with supported organisms
print "<table class='sortable' width='600'>";
foreach my $row (@organism_rows) {
  my @fields = split("\t", $row);
  print "<tr>\n";
  if ($fields[0] =~'^#') {
    $fields[0]=~s/#//;
    foreach my $field (@fields) {
      print "<th>", $field, "</th>\n";
    }
  } else {
    foreach my $field (@fields) {
      ## Convert full RSAT path to relative
      $field = &RSAT::util::hide_RSAT_path($field); 

      ## Replace data directory by a link
      if ($field =~ /public_html\/(data\/\S*)/) {
	$data_dir = $1;
	$field = join("", "<a href='", $data_dir, "'>data</a>");
      }
      print "<td>", $field, "</td>\n";
    }
  }
  print "</tr>\n";
}
print "</table>";

#$command .= " ".$parameters;
#&ReportWebCommand($command) if ($ENV{rsat_echo} >= 1);

## Add link to the data folder of each organism
#$command .= " | awk '{print \$0\"\t<a href=$ENV{rsat_www}/data/genomes/\"\$1\"/>data</a>\"}'";
#$command .= " | perl -pe 's|<a href=.+/data/genomes/;/>data</a>||'";

# open RESULT, "$command|";

# ### Print result on the web page
# print "<CENTER>\n";
# &PrintHtmlTable(RESULT, $tmp_file_path, 0, 10000);
# print "</CENTER>\n";

# close(RESULT);

print '<hr size=3>';


print $query->end_html;

exit(0);

