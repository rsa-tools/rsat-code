package MyTypes::CompareMatricesRequest;
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
my %matrix_1_of :ATTR(:get<matrix_1>);
my %matrix_2_of :ATTR(:get<matrix_2>);
my %matrix_of :ATTR(:get<matrix>);
my %tmp_matrix1_infile_of :ATTR(:get<tmp_matrix1_infile>);
my %tmp_matrix2_infile_of :ATTR(:get<tmp_matrix2_infile>);
my %tmp_matrix_infile_of :ATTR(:get<tmp_matrix_infile>);
my %format1_of :ATTR(:get<format1>);
my %format2_of :ATTR(:get<format2>);
my %format_of :ATTR(:get<format>);
my %background_model_of :ATTR(:get<background_model>);
my %tmp_background_infile_of :ATTR(:get<tmp_background_infile>);
my %background_format_of :ATTR(:get<background_format>);
my %top1_of :ATTR(:get<top1>);
my %top2_of :ATTR(:get<top2>);
my %output_prefix_of :ATTR(:get<output_prefix>);
my %mode_of :ATTR(:get<mode>);
my %distinct_of :ATTR(:get<distinct>);
my %strand_of :ATTR(:get<strand>);
my %matrix_id_of :ATTR(:get<matrix_id>);
my %return_of :ATTR(:get<return>);
my %sort_of :ATTR(:get<sort>);
my %lth_of :ATTR(:get<lth>);
my %uth_of :ATTR(:get<uth>);

__PACKAGE__->_factory(
    [ qw(        output
        matrix_1
        matrix_2
        matrix
        tmp_matrix1_infile
        tmp_matrix2_infile
        tmp_matrix_infile
        format1
        format2
        format
        background_model
        tmp_background_infile
        background_format
        top1
        top2
        output_prefix
        mode
        distinct
        strand
        matrix_id
        return
        sort
        lth
        uth

    ) ],
    {
        'output' => \%output_of,
        'matrix_1' => \%matrix_1_of,
        'matrix_2' => \%matrix_2_of,
        'matrix' => \%matrix_of,
        'tmp_matrix1_infile' => \%tmp_matrix1_infile_of,
        'tmp_matrix2_infile' => \%tmp_matrix2_infile_of,
        'tmp_matrix_infile' => \%tmp_matrix_infile_of,
        'format1' => \%format1_of,
        'format2' => \%format2_of,
        'format' => \%format_of,
        'background_model' => \%background_model_of,
        'tmp_background_infile' => \%tmp_background_infile_of,
        'background_format' => \%background_format_of,
        'top1' => \%top1_of,
        'top2' => \%top2_of,
        'output_prefix' => \%output_prefix_of,
        'mode' => \%mode_of,
        'distinct' => \%distinct_of,
        'strand' => \%strand_of,
        'matrix_id' => \%matrix_id_of,
        'return' => \%return_of,
        'sort' => \%sort_of,
        'lth' => \%lth_of,
        'uth' => \%uth_of,
    },
    {
        'output' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'matrix_1' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'matrix_2' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'matrix' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'tmp_matrix1_infile' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'tmp_matrix2_infile' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'tmp_matrix_infile' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'format1' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'format2' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'format' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'background_model' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'tmp_background_infile' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'background_format' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'top1' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'top2' => 'SOAP::WSDL::XSD::Typelib::Builtin::int',
        'output_prefix' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'mode' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'distinct' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'strand' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'matrix_id' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'return' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'sort' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'lth' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
        'uth' => 'SOAP::WSDL::XSD::Typelib::Builtin::string',
    },
    {

        'output' => 'output',
        'matrix_1' => 'matrix_1',
        'matrix_2' => 'matrix_2',
        'matrix' => 'matrix',
        'tmp_matrix1_infile' => 'tmp_matrix1_infile',
        'tmp_matrix2_infile' => 'tmp_matrix2_infile',
        'tmp_matrix_infile' => 'tmp_matrix_infile',
        'format1' => 'format1',
        'format2' => 'format2',
        'format' => 'format',
        'background_model' => 'background_model',
        'tmp_background_infile' => 'tmp_background_infile',
        'background_format' => 'background_format',
        'top1' => 'top1',
        'top2' => 'top2',
        'output_prefix' => 'output_prefix',
        'mode' => 'mode',
        'distinct' => 'distinct',
        'strand' => 'strand',
        'matrix_id' => 'matrix_id',
        'return' => 'return',
        'sort' => 'sort',
        'lth' => 'lth',
        'uth' => 'uth',
    }
);

} # end BLOCK







1;


=pod

=head1 NAME

MyTypes::CompareMatricesRequest

=head1 DESCRIPTION

Perl data type class for the XML Schema defined complexType
CompareMatricesRequest from the namespace urn:RSATWS.

Parameters for the operation compare_matrices




=head2 PROPERTIES

The following properties may be accessed using get_PROPERTY / set_PROPERTY
methods:

=over

=item * output


=item * matrix_1


=item * matrix_2


=item * matrix


=item * tmp_matrix1_infile


=item * tmp_matrix2_infile


=item * tmp_matrix_infile


=item * format1


=item * format2


=item * format


=item * background_model


=item * tmp_background_infile


=item * background_format


=item * top1


=item * top2


=item * output_prefix


=item * mode


=item * distinct


=item * strand


=item * matrix_id


=item * return


=item * sort


=item * lth


=item * uth




=back


=head1 METHODS

=head2 new

Constructor. The following data structure may be passed to new():

 { # MyTypes::CompareMatricesRequest
   output =>  $some_value, # string
   matrix_1 =>  $some_value, # string
   matrix_2 =>  $some_value, # string
   matrix =>  $some_value, # string
   tmp_matrix1_infile =>  $some_value, # string
   tmp_matrix2_infile =>  $some_value, # string
   tmp_matrix_infile =>  $some_value, # string
   format1 =>  $some_value, # string
   format2 =>  $some_value, # string
   format =>  $some_value, # string
   background_model =>  $some_value, # string
   tmp_background_infile =>  $some_value, # string
   background_format =>  $some_value, # string
   top1 =>  $some_value, # int
   top2 =>  $some_value, # int
   output_prefix =>  $some_value, # string
   mode =>  $some_value, # string
   distinct =>  $some_value, # string
   strand =>  $some_value, # string
   matrix_id =>  $some_value, # string
   return =>  $some_value, # string
   sort =>  $some_value, # string
   lth =>  $some_value, # string
   uth =>  $some_value, # string
 },




=head1 AUTHOR

Generated by SOAP::WSDL

=cut

