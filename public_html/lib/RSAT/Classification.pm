###############################################################
#
# Class Classification
#
package RSAT::Classification;

use RSAT::GenericObject;
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
 Usage    : $classificaion->->read_from_file($inputfile, $input_format)
 Function : Read the Classification from a text file. 
 Supported formats: mcl

=cut

sub read_from_file {
    my ($self, $inputfile, $input_format) = @_;
    if ($input_format eq "mcl") {
	$self->read_mcl($inputfile);
    } else {
	&RSAT::error::FatalError(join ("\t", "Classification::read_from_file", $input_format, "is not a supported input format"));
    }
}


################################################################
=pod

=item B<read_mcl>

Read the classification from a mcl text file.


 Title    : read_mcl
 Usage    : $classificaion->->read_mcl($inputfile)
 Function : Read the Classification from a mcl text file.

=cut

sub read_mcl {
    my ($self, $inputfile) = @_;
    ($main::in) = &RSAT::util::OpenInputFile($inputfile);
    
    ## Load the classification
    my $family_number = 0;
    while (<$main::in>) {
	next unless (/\S/);
	chomp;
	$family_number++;
	my @labels = split(/\s+/);
	my $family = new RSAT::Family(name=>$family_number);
	&RSAT::message::Debug("Reading", $family, $family_number);
	foreach my $label (@labels) {
	    $family->new_member($label);
	}
	$self->push_attribute("families", $family);
    }
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
    my @families = $self->get_attribute("families");
    my $string = "";
    $string .= join ("\t", "#element", "class");
    foreach my $family (@families) {
	my $family_name = $family->get_attribute("name");
	&RSAT::message::Debug("Printing", $family, $family_name);
	foreach my $member ($family->get_members()) {
	    $string .= join ("\t", $member, $family_name);
	    $string .= "\n";
	}
    }
    return $string;
}

return 1;

__END__

