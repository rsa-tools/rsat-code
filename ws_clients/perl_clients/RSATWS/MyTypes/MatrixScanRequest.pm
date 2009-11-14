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
my %sequence_of :ATTR(:get<sequence>);
my %tmp_sequence_infile_of :ATTR(:get<tmp_sequence_infile>);
my %matrix_of :ATTR(:get<matrix>);
my %tmp_matrix_infile_of :ATTR(:get<tmp_matrix_infile>);
my %sequence_format_of :ATTR(:get<sequence_format>);
my %matrix_format_of :ATTR(:get<matrix_format>);
my %n_treatment_of :ATTR(:get<n_treatment>);
my %consensus_name_of :ATTR(:get<consensus_name>);
my %pseudo_of :ATTR(:get<pseudo>);
my %equi_pseudo_of :ATTR(:get<equi_pseudo>);
my %top_matrices_of :ATTR(:get<top_matrices>);
my %background_model_of :ATTR(:get<background_model>);
my %tmp_background_infile_of :ATTR(:get<tmp_background_infile>);
my %organism_of :ATTR(:get<organism>);
my %background_of :ATTR(:get<background>);
my %background_input_of :ATTR(:get<background_input>);
my %background_window_of :ATTR(:get<background_window>);
my %markov_of :ATTR(:get<markov>);
my %background_pseudo_of :ATTR(:get<background_pseudo>);
my %return_fields_of :ATTR(:get<return_fields>);
my %sort_distrib_of :ATTR(:get<sort_distrib>);
my %lth_of :ATTR(:get<lth>);
my %uth_of :ATTR(:get<uth>);
my %str_of :ATTR(:get<str>);
my %verbosity_of :ATTR(:get<verbosity>);
my %origin_of :ATTR(:get<origin>);
my %decimals_of :ATTR(:get<decimals>);
my %crer_ids_of :ATTR(:get<crer_ids>);

__PACKAGE__->_factory(
    [ qw(        output
        sequence
        tmp_sequence_infile
        matrix
        tmp_matrix_infile
        sequence_format
        matrix_format
        n_treatment
        consensus_name
        pseudo
        equi_pseudo
        top_matrices
        background_model
        tmp_background_infile
        organism
        background
        background_input
        background_window
        markov
        background_pseudo
        return_fields
        sort_distrib
        lth
        uth
        str
        verbosity
        origin
        decimals
        crer_ids

    ) ],
    {
        'output' => \%output_of,
        'sequence' => \%sequence_of,
        'tmp_sequence_infile' => \%tmp_sequence_infile_of,
        'matrix' => \%matrix_of,
        'tmp_matrix_infile' => \%tmp_matrix_infile_of,
        'sequence_format' => \%sequence_format_of,
        'matrix_format' => \%matrix_format_of,
        'n_treatment' => \%n_treatment_of,
        'consensus_name' => \%consensus_name_of,
        'pseudo' => \%pseudo_of,
        'equi_pseudo' => \%equi_pseudo_of,
        'top_matrices' => \%top_matrices_of,
        'background_model' => \%background_model_of,
        'tmp_background_infile' => \%tmp_background_infile_of,
        'organism' => \%organism_of,
        'background' => \%background_of,
        'background_input' => \%background_input_of,
        'background_window' => \%background_window_of,
        'markov' => \%markov_of,
        'background_pseudo' => \%background_pseudo_of,
        'return_fields' => \%return_fields_of,
        'sort_distrib' => \%sort_distrib_of,
        'lth' => \%lth_of,
        'uth' => \%uth_of,
        'str' => \%str_of,
        'verbosity' => \%verbosity_of,
        'origin' => \%origin_of,
        'decimals' => \%decimals_of,
        'crer_ids' => \%crer_ids_of,
    },
    {
        'output' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'sequence' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'tmp_sequence_infile' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'matrix' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'tmp_matrix_infile' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'sequence_format' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'matrix_format' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'n_treatment' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'consensus_name' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'pseudo' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'equi_pseudo' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'top_matrices' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'background_model' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'tmp_background_infile' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'organism' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'background' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'background_input' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'background_window' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'markov' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'background_pseudo' => 'SOAP::WSDL::XSD::Typelib::Builtin::float',
        'return_fields' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'sort_distrib' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'lth' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'uth' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'str' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'verbosity' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'origin' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'decimals' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'crer_ids' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
    },
    {

        'output' => 'output',
        'sequence' => 'sequence',
        'tmp_sequence_infile' => 'tmp_sequence_infile',
        'matrix' => 'matrix',
        'tmp_matrix_infile' => 'tmp_matrix_infile',
        'sequence_format' => 'sequence_format',
        'matrix_format' => 'matrix_format',
        'n_treatment' => 'n_treatment',
        'consensus_name' => 'consensus_name',
        'pseudo' => 'pseudo',
        'equi_pseudo' => 'equi_pseudo',
        'top_matrices' => 'top_matrices',
        'background_model' => 'background_model',
        'tmp_background_infile' => 'tmp_background_infile',
        'organism' => 'organism',
        'background' => 'background',
        'background_input' => 'background_input',
        'background_window' => 'background_window',
        'markov' => 'markov',
        'background_pseudo' => 'background_pseudo',
        'return_fields' => 'return_fields',
        'sort_distrib' => 'sort_distrib',
        'lth' => 'lth',
        'uth' => 'uth',
        'str' => 'str',
        'verbosity' => 'verbosity',
        'origin' => 'origin',
        'decimals' => 'decimals',
        'crer_ids' => 'crer_ids',
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

Parameters for the operation matrix_scan




=head2 PROPERTIES

The following properties may be accessed using get_PROPERTY / set_PROPERTY
methods:

=over

=item * output


=item * sequence


=item * tmp_sequence_infile


=item * matrix


=item * tmp_matrix_infile


=item * sequence_format


=item * matrix_format


=item * n_treatment


=item * consensus_name


=item * pseudo


=item * equi_pseudo


=item * top_matrices


=item * background_model


=item * tmp_background_infile


=item * organism


=item * background


=item * background_input


=item * background_window


=item * markov


=item * background_pseudo


=item * return_fields


=item * sort_distrib


=item * lth


=item * uth


=item * str


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
   sequence =>  $some_value, # string
   tmp_sequence_infile =>  $some_value, # string
   matrix =>  $some_value, # string
   tmp_matrix_infile =>  $some_value, # string
   sequence_format =>  $some_value, # string
   matrix_format =>  $some_value, # string
   n_treatment =>  $some_value, # string
   consensus_name =>  $some_value, # string
   pseudo =>  $some_value, # int
   equi_pseudo =>  $some_value, # int
   top_matrices =>  $some_value, # int
   background_model =>  $some_value, # string
   tmp_background_infile =>  $some_value, # string
   organism =>  $some_value, # string
   background =>  $some_value, # string
   background_input =>  $some_value, # int
   background_window =>  $some_value, # int
   markov =>  $some_value, # int
   background_pseudo =>  $some_value, # float
   return_fields =>  $some_value, # string
   sort_distrib =>  $some_value, # int
   lth =>  $some_value, # string
   uth =>  $some_value, # string
   str =>  $some_value, # int
   verbosity =>  $some_value, # int
   origin =>  $some_value, # string
   decimals =>  $some_value, # int
   crer_ids =>  $some_value, # int
 },




=head1 AUTHOR

Generated by SOAP::WSDL

=cut

