package MyTypes::PeakMotifsRequest;
use strict;
use warnings;


__PACKAGE__->_set_element_form_qualified(0);

sub get_xmlns { 'urn:RSATWS' };

our $XML_ATTRIBUTE_CLASS;
undef $XML_ATTRIBUTE_CLASS;

sub __get_attr_class {
    return $XML_ATTRIBUTE_CLASS;
}

use Class::Std::Fast::Storable constructor => 'none';
use base qw(SOAP::WSDL::XSD::Typelib::ComplexType);

Class::Std::initialize();

{ # BLOCK to scope variables

my %output_of :ATTR(:get<output>);
my %verbosity_of :ATTR(:get<verbosity>);
my %test_of :ATTR(:get<test>);
my %tmp_test_infile_of :ATTR(:get<tmp_test_infile>);
my %control_of :ATTR(:get<control>);
my %tmp_control_infile_of :ATTR(:get<tmp_control_infile>);
my %max_seq_length_of :ATTR(:get<max_seq_length>);
my %max_motif_number_of :ATTR(:get<max_motif_number>);
my %ref_motif_of :ATTR(:get<ref_motif>);
my %top_peaks_of :ATTR(:get<top_peaks>);
my %min_length_of :ATTR(:get<min_length>);
my %max_length_of :ATTR(:get<max_length>);
my %markov_of :ATTR(:get<markov>);
my %min_markov_of :ATTR(:get<min_markov>);
my %max_markov_of :ATTR(:get<max_markov>);
my %noov_of :ATTR(:get<noov>);
my %class_int_of :ATTR(:get<class_int>);
my %str_of :ATTR(:get<str>);
my %graph_title_of :ATTR(:get<graph_title>);
my %image_format_of :ATTR(:get<image_format>);
my %disco_of :ATTR(:get<disco>);
my %source_of :ATTR(:get<source>);
my %task_of :ATTR(:get<task>);

__PACKAGE__->_factory(
    [ qw(        output
        verbosity
        test
        tmp_test_infile
        control
        tmp_control_infile
        max_seq_length
        max_motif_number
        ref_motif
        top_peaks
        min_length
        max_length
        markov
        min_markov
        max_markov
        noov
        class_int
        str
        graph_title
        image_format
        disco
        source
        task

    ) ],
    {
        'output' => \%output_of,
        'verbosity' => \%verbosity_of,
        'test' => \%test_of,
        'tmp_test_infile' => \%tmp_test_infile_of,
        'control' => \%control_of,
        'tmp_control_infile' => \%tmp_control_infile_of,
        'max_seq_length' => \%max_seq_length_of,
        'max_motif_number' => \%max_motif_number_of,
        'ref_motif' => \%ref_motif_of,
        'top_peaks' => \%top_peaks_of,
        'min_length' => \%min_length_of,
        'max_length' => \%max_length_of,
        'markov' => \%markov_of,
        'min_markov' => \%min_markov_of,
        'max_markov' => \%max_markov_of,
        'noov' => \%noov_of,
        'class_int' => \%class_int_of,
        'str' => \%str_of,
        'graph_title' => \%graph_title_of,
        'image_format' => \%image_format_of,
        'disco' => \%disco_of,
        'source' => \%source_of,
        'task' => \%task_of,
    },
    {
        'output' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'verbosity' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'test' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'tmp_test_infile' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'control' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'tmp_control_infile' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'max_seq_length' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'max_motif_number' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'ref_motif' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'top_peaks' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'min_length' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'max_length' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'markov' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'min_markov' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'max_markov' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'noov' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'class_int' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'str' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'graph_title' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'image_format' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'disco' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'source' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'task' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
    },
    {

        'output' => 'output',
        'verbosity' => 'verbosity',
        'test' => 'test',
        'tmp_test_infile' => 'tmp_test_infile',
        'control' => 'control',
        'tmp_control_infile' => 'tmp_control_infile',
        'max_seq_length' => 'max_seq_length',
        'max_motif_number' => 'max_motif_number',
        'ref_motif' => 'ref_motif',
        'top_peaks' => 'top_peaks',
        'min_length' => 'min_length',
        'max_length' => 'max_length',
        'markov' => 'markov',
        'min_markov' => 'min_markov',
        'max_markov' => 'max_markov',
        'noov' => 'noov',
        'class_int' => 'class_int',
        'str' => 'str',
        'graph_title' => 'graph_title',
        'image_format' => 'image_format',
        'disco' => 'disco',
        'source' => 'source',
        'task' => 'task',
    }
);

} # end BLOCK







1;


=pod

=head1 NAME

MyTypes::PeakMotifsRequest

=head1 DESCRIPTION

Perl data type class for the XML Schema defined complexType
PeakMotifsRequest from the namespace urn:RSATWS.

Parameters for the operation chip_motifs.




=head2 PROPERTIES

The following properties may be accessed using get_PROPERTY / set_PROPERTY
methods:

=over

=item * output


=item * verbosity


=item * test


=item * tmp_test_infile


=item * control


=item * tmp_control_infile


=item * max_seq_length


=item * max_motif_number


=item * ref_motif


=item * top_peaks


=item * min_length


=item * max_length


=item * markov


=item * min_markov


=item * max_markov


=item * noov


=item * class_int


=item * str


=item * graph_title


=item * image_format


=item * disco


=item * source


=item * task




=back


=head1 METHODS

=head2 new

Constructor. The following data structure may be passed to new():

 { # MyTypes::PeakMotifsRequest
   output =>  $some_value, # string
   verbosity =>  $some_value, # int
   test =>  $some_value, # string
   tmp_test_infile =>  $some_value, # string
   control =>  $some_value, # string
   tmp_control_infile =>  $some_value, # string
   max_seq_length =>  $some_value, # int
   max_motif_number =>  $some_value, # int
   ref_motif =>  $some_value, # string
   top_peaks =>  $some_value, # int
   min_length =>  $some_value, # int
   max_length =>  $some_value, # int
   markov =>  $some_value, # int
   min_markov =>  $some_value, # int
   max_markov =>  $some_value, # int
   noov =>  $some_value, # int
   class_int =>  $some_value, # int
   str =>  $some_value, # int
   graph_title =>  $some_value, # string
   image_format =>  $some_value, # string
   disco =>  $some_value, # string
   source =>  $some_value, # string
   task =>  $some_value, # string
 },




=head1 AUTHOR

Generated by SOAP::WSDL

=cut

