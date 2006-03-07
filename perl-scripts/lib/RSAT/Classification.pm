###############################################################
#
# Class Classification
#
package RSAT::Classification;

use RSAT::GenericObject;
use RSAT::util;
use RSAT::Family;
use RSAT::error;
@ISA = qw( RSAT::GenericObject );

### class attributes

=pod

=head1 NAME

    RSAT::Classification

=head1 DESCRIPTION

This class stores a classification, i.e. a set of classes, each
containing a set of objects.

The classification can be a partition (i.e. each objects belogs to one
and only one class) or contain multiple assignations.

The class is equipped with methods for importing and exporting a
classification from various formats.

Hierarchical classifications can be loaded, but the inter-class
hierarchy is not treated in the current version.

=cut


################################################################
=pod

=item new()

Create a new Classification.

=cut
sub new {
    my ($class, %args) = @_;
    my $object = bless {
	}, $class;
    return $object;
}


################################################################
=pod

=item B<read_from_file>

Read the classification from a text file. 


 Title    : read_from_file
 Usage    : $classificaion->->read_from_file($input_file, $input_format)
 Function : Read the Classification from a text file. 
 Supported formats: mcl

=cut

sub read_from_file {
    my ($self, $input_file, $input_format) = @_;
    if ($input_format eq "mcl") {
	$self->read_mcl($input_file);
    } elsif ($input_format eq "tab") {
	$self->read_tab($input_file);
    } else {
	&RSAT::error::FatalError(join ("\t", "Classification::read_from_file", $input_format, "is not a supported input format"));
    }
}


################################################################
=pod

=item B<read_mcl>

Read the classification from a mcl text file.


 Title    : read_mcl
 Usage    : $classificaion->->read_mcl($input_file)
 Function : Read the Classification from a mcl text file.

=cut

sub read_mcl {
    my ($self, $input_file) = @_;
    ($main::in) = &RSAT::util::OpenInputFile($input_file);
    
    ## Load the classification
    my $class_number = 0;
    while (<$main::in>) {
	next unless (/\S/);
	chomp;
	$class_number++;
	my @labels = split(/\s+/);
	my $class = new RSAT::Family(name=>$class_number);
	&RSAT::message::Debug("Reading", $class, $class_number);
	foreach my $label (@labels) {
	    $class->new_member($label);
	}
	$self->push_attribute("classes", $class);
    }
}


################################################################
=pod

=item B<read_tab>

Read the classification from a tab-delimited text file.  Each row
describes one class membership.  By default, the first column
indicates the member, the second column the class. 

 Title    : read_tab
 Usage    : $classificaion->->read_tab($input_file)
 Function : Read the Classification from a tab-delimited text file.

The columns can be specified by specifying additional arguments. For
instance, if the input is a tab-delimited file with members in the
third and classes in the seventh column :

 Usage :  $classificaion->->read_tab($input_file,member_column=>3,class_column=>7)

A score column can optionally be specified with the argument score_column. 

 Usage :  $classificaion->->read_tab($input_file,score_column=3)

=cut

sub read_tab {
    my ($self, $input_file, %args) = @_;
    my ($in) = &RSAT::util::OpenInputFile($input_file);
    my %class = (); ## classes indexed by name
    $member_column = $args{member_column} || 1;
    $class_column = $args{class_column} || 2;
    $score_column = $args{score_column} || 0;

    ## Verbosity
    if ($main::verbose >= 0) {
	&RSAT::message::TimeWarn(join("\t", 
				  "Reading classification from tab file", $input_file) ) if ($main::verbose >= 2);
	&RSAT::message::Info(join("\t", 
				  "member column: ".$member_column,
				  "class column: ".$class_column,
				  "score column: ".$score_column,
				 )) if ($main::verbose >= 3);
    }

    ## Load the classification
    my $line = 0;
    while (<$in>) {
	$line++;
	next if (/^;/);
	next if (/^\#/);
	next unless (/\S/);
	chomp();
	s/\r$//;
	my @fields = split /\t/;
	
	### class member
	$member_name = &RSAT::util::trim(uc($fields[$member_column -1]));
	unless ($member_name =~ /\S/) {
	    &RSAT::message::Warning(join ("\t", "Error class file", 
					  $class_file,  "line", 
					  $line, "member not specified")) if ($main::verbose >= 1);
	    next;
	}
	
	### class name
	$class_name = &RSAT::util::trim($fields[$class_column - 1]);
	unless ($class_name) {
	    &RSAT::message::Warning(join("\t", "Error class file", 
					 $class_file,  "line", 
					 $line, "class not specified")) if ($main::verbose >= 1);
	    next;
	}
	
	#### create a new class if required
	unless ($class{$class_name}) {
	    $class{$class_name} = new RSAT::Family(name=>$class_name);
	    $self->push_attribute("classes", $class{$class_name});
	}
	
	## Add the member to the class
	if ($score_column) {
	    ### read score
	    local $score = &RSAT::util::trim($fields[$score_column - 1]);
	    unless (&IsReal($score)) {
		&RSAT::error::FatalError(join("\t", $score, "Invalid score (must be a Real number).", 
					      "class file", $class_file,  
					      "line", $line,
					      "member", $member_name,
					      "class", $class_name,
					     ));
	    }
	    $class{$class_name}->new_member($member_name, 0, score=>$score);
#		&RSAT::message::Debug("member score", $class_name, $member_name, $score) if ($main::verbose >= 0);
	} else {
	    $class{$class_name}->new_member($member_name);
	}
	&RSAT::message::Warning( join ("\t",  ";", $class_name,
				       $member_name,
				       $name{$id}) ) if ($main::verbose >= 5);
#	    &RSAT::message::Debug("Scores", $class_name, $class{$class_name}->get_attribute("scores")) if ($main::verbose >= 0);
	
    }
    close $in if ($input_file);

}


################################################################
=pod

=item B<to_text>

Convert the classification into a string for export to a text file.


 Title    : to_text
 Usage    : $classificaion->->to_text($output_format)
 Supported formats: tab

=cut

sub to_text {
    my ($self, $output_format) = @_;
    if ($output_format eq "tab") {
	$self->to_tab();
    } elsif ($output_format eq "mcl") {
	$self->to_mcl();
    } else {
	&RSAT::error::FatalError(join ("\t", "Classification::to_text", $output_format, "is not a supported output format"));
    }
}

################################################################
=pod

=item B<to_tab>

Convert the classification to a tab-delimited string. 
One row per element-class association. 
First column contains the element.
Second column contains the class.

 Title    : to_tab
 Usage    : $classificaion->->to_tab()

=cut

sub to_tab {
    my ($self) = @_;
    my @classes = $self->get_attribute("classes");
    my $string = "";
    $string .= join ("\t", "#element", "class");
    foreach my $class (@classes) {
	my $class_name = $class->get_attribute("name");
	&RSAT::message::Debug("Printing", $class, $class_name);
	foreach my $member ($class->get_members()) {
	    $string .= join ("\t", $member, $class_name);
	    $string .= "\n";
	}
    }
    return $string;
}

################################################################
=pod

=item B<to_mcl>

Convert the classification to the mcl cluster format. MCL is a
graph-based clustering algorithm developed by Stijn Van Dongen
(http://micans.org/mcl/).


One row per cluster. Each row provides a list of elements, separated
by tab characters. 

 Title    : to_mcl
 Usage    : $classificaion->->to_mcl()

=cut

sub to_mcl {
    my ($self) = @_;
    my @classes = $self->get_attribute("classes");
    my $string = "";
    foreach my $class (@classes) {
	my $class_name = $class->get_attribute("name");
	&RSAT::message::Debug("Printing", $class, $class_name) if ($main::verbose >= 5);
	my @members = $class->get_members();
	$string .= join ("\t", @members);
	$string .= "\n";
    }
    return $string;
}

return 1;

__END__

