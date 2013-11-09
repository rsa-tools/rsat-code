#################################################################
#
# Quick test for parsing transduction signal. let us start with references
#

### add the program's directory to the lib path
BEGIN {
    if ($0 =~ /([^(\/)]+)$/) {
	push (@INC, "$`"); ### add the program's directory to the lib path
    }
}
require "config.pl";
require "lib/load_classes.pl";
require "lib/parsing_util.pl";


################################################################
#### bibliographic references
package SIGTRANS::Reference;
{
    @ISA = qw ( classes::DatabaseObject );
    %_attribute_cardinality = (
			       id=>"SCALAR",
			       description=>"SCALAR",
			       datastarID=>"SCALAR",
			       pubmedID=>"SCALAR",
			       comments=>"ARRAY"
			       );
}


################################################################
#### biochemical entities
package SIGTRANS::Entity;
{
    @ISA = qw ( classes::DatabaseObject );
    %_attribute_cardinality = (
			       id=>"SCALAR",
			       description=>"SCALAR",
			       names=>"ARRAY"
			       );
}


################################################################
#### transformation


################################################################
#### control


################################################################
#### main package
package main;
{
    #### initialisation
    $schema = "sigtrans";
    $user = "sigtrans";
    $password = "sigtrans";

    #### bibliographic references
    $dir{input} = "/win/amaze/amaze_team/sandra/exported_tables/export_20021223/";

    $export_subdir = "sigtrans";
    $dir{output} = "$parsed_data/$export_subdir/$delivery_date";


    &ReadArguments();

    #### output directory
    &CheckOutputDir();

    $in_file{references} = $dir{input}."/reference_description.txt";


    ################################################################
    #### define class holders 

    #### bibliographic references
    $references = classes::ClassFactory->new_class(object_type=>"SIGTRANS::Reference",
						prefix=>"ref_");

    $references->set_out_fields( qw (id pubmedID datastarID description) );

    #### biochemical entities
    $entities =  classes::ClassFactory->new_class(object_type=>"SIGTRANS::Entity",
					       prefix=>"ent_");

    
    #### some verbosity 
    &DefaultVerbose() if ($verbose >= 1);

    &ParseReferences($references, $in_file{references});

    foreach my $holder ($references,
			$entities) {
	$suffix="";
	$holder->dump_tables($suffix, 0, $dir{output});
	$holder->generate_sql(schema=>$main::schema, 
			      user=>$main::user,
			      password=>$main::password,
			      dir=>"$dir{output}/sql_scripts",
			      host=>$main::host,
			      dbms=>$main::dbms
			      );
    }

    exit(0);
}

################################################################
#### parse bibliographic references
sub ParseReferences {
    my ($class_holder, $inputfile) = @_;
    
    unless (-e $inputfile) {
	die "Error: the file $inputfile does not exist\n";
    }

    open REF, $inputfile;
    while (<REF>) {
	s/\r//;
	chomp;
	my @fields = split "\t";
	$new_ref = $class_holder->new_object();
	$new_ref->new_attribute_value("description", $fields[0]);
	$new_ref->new_attribute_value("datastarID", $fields[1]);
	$new_ref->new_attribute_value("pubmedID", $fields[2]);
	$new_ref->new_attribute_value("comments", $fields[3]);
    }
    close REF;

}



################################################################
### read arguments from the command line
sub ReadArguments {
    for my $a (0..$#ARGV) {

	&ReadGenericOptions($a);
	    
	#### input directory
	if ($ARGV[$a] =~ /^-indir/) {
	    $main::dir{input} = $ARGV[$a+1];


	}
    }
}


### print the help message 
### when the program is called with the -h or -help option
sub PrintHelp {
  open HELP, "| more";
  print <<EndHelp;
NAME
        parse_signal_transduction.pl

DESCRIPTION

	Parse signal transduction data, exported in tab-delimited
	format from the Access database.

AUTHOR
	Jacques van Helden (jvanheld\@ucmb.ulb.ac.be)  

OPTIONS	
$generic_option_message
	-indir	input directory
		This directory should contain the .txt files resulting
		from the export of the Access database.
EndHelp
  close HELP;
}

  
