package MyTypes::MatrixScanRequest;
use strict;
use warnings;


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
my %sequence_file_of :ATTR(:get<sequence_file>);
my %matrix_file_of :ATTR(:get<matrix_file>);
my %matrix_format_of :ATTR(:get<matrix_format>);
my %matrix_list_of :ATTR(:get<matrix_list>);
my %top_matrices_of :ATTR(:get<top_matrices>);
my %background_of :ATTR(:get<background>);
my %background_input_of :ATTR(:get<background_input>);
my %background_window_of :ATTR(:get<background_window>);
my %markov_of :ATTR(:get<markov>);
my %background_pseudo_of :ATTR(:get<background_pseudo>);
my %return_fields_of :ATTR(:get<return_fields>);
my %lth_of :ATTR(:get<lth>);
my %uth_of :ATTR(:get<uth>);
my %both_strand_of :ATTR(:get<both_strand>);
my %single_strand_of :ATTR(:get<single_strand>);
my %verbosity_of :ATTR(:get<verbosity>);
my %origin_of :ATTR(:get<origin>);
my %decimals_of :ATTR(:get<decimals>);
my %crer_ids_of :ATTR(:get<crer_ids>);

__PACKAGE__->_factory(
    [ qw(
        output
        sequence_file
        matrix_file
        matrix_format
        matrix_list
        top_matrices
        background
        background_input
        background_window
        markov
        background_pseudo
        return_fields
        lth
        uth
        both_strand
        single_strand
        verbosity
        origin
        decimals
        crer_ids
    ) ],
    {
        'output' => \%output_of,
        'sequence_file' => \%sequence_file_of,
        'matrix_file' => \%matrix_file_of,
        'matrix_format' => \%matrix_format_of,
        'matrix_list' => \%matrix_list_of,
        'top_matrices' => \%top_matrices_of,
        'background' => \%background_of,
        'background_input' => \%background_input_of,
        'background_window' => \%background_window_of,
        'markov' => \%markov_of,
        'background_pseudo' => \%background_pseudo_of,
        'return_fields' => \%return_fields_of,
        'lth' => \%lth_of,
        'uth' => \%uth_of,
        'both_strand' => \%both_strand_of,
        'single_strand' => \%single_strand_of,
        'verbosity' => \%verbosity_of,
        'origin' => \%origin_of,
        'decimals' => \%decimals_of,
        'crer_ids' => \%crer_ids_of,
    },
    {
        'output' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'sequence_file' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'matrix_file' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'matrix_format' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'matrix_list' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'top_matrices' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'background' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'background_input' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'background_window' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'markov' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'background_pseudo' => 'SOAP::WSDL::XSD::Typelib::Builtin::float',
        'return_fields' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'lth' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'uth' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'both_strand' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'single_strand' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'verbosity' => 'SOAP::WSDL::XSD::Typelib::Builtin::float',
        'origin' => 'SOAP::WSDL::XSD::Typelib::Builtin::float',
        'decimals' => 'SOAP::WSDL::XSD::Typelib::Builtin::float',
        'crer_ids' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
    }
);

} # end BLOCK







1;


=pod

=head1 NAME

MyTypes::MatrixScanRequest

=head1 DESCRIPTION

Perl data type class for the XML Schema defined complexType
MatrixScanRequest from the namespace urn:RSATWS.

Parameters for the operation matrix scan




=head2 PROPERTIES

The following properties may be accessed using get_PROPERTY / set_PROPERTY
methods:

=over

=item * output

=item * sequence_file

=item * matrix_file

=item * matrix_format

=item * matrix_list

=item * top_matrices

=item * background

=item * background_input

=item * background_window

=item * markov

=item * background_pseudo

=item * return_fields

=item * lth

=item * uth

=item * both_strand

=item * single_strand

=item * verbosity

=item * origin

=item * decimals

=item * crer_ids



=back


=head1 METHODS

=head2 new

Constructor. The following data structure may be passed to new():

 { # MyTypes::MatrixScanRequest
   output =>  $some_value, # string
   sequence_file =>  $some_value, # string
   matrix_file =>  $some_value, # string
   matrix_format =>  $some_value, # string
   matrix_list =>  $some_value, # string
   top_matrices =>  $some_value, # int
   background =>  $some_value, # string
   background_input =>  $some_value, # int
   background_window =>  $some_value, # int
   markov =>  $some_value, # int
   background_pseudo =>  $some_value, # float
   return_fields =>  $some_value, # string
   lth =>  $some_value, # string
   uth =>  $some_value, # string
   both_strand =>  $some_value, # int
   single_strand =>  $some_value, # int
   verbosity =>  $some_value, # float
   origin =>  $some_value, # float
   decimals =>  $some_value, # float
   crer_ids =>  $some_value, # int
 },




=head1 AUTHOR

Generated by SOAP::WSDL

=cut

